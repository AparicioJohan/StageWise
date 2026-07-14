#' Prepare data for BLUP 
#' 
#' Prepare data for BLUP
#' 
#' The \code{method} argument can be used to control how the linear system is solved. "MME" leads to inversion of the MME coefficient matrix, while "Vinv" leads to inversion of the overall var-cov matrix for the response vector. If NULL, the software uses whichever method involves inverting the smaller matrix. If the number of random effects (m) is less than the number of BLUEs (n), "MME" is used.
#' 
#' For the multi-location model, if all of the environments for a location are masked, the average of the other locations is used when computing average fixed effects.
#' 
#' @param data data frame of BLUEs from Stage 1 
#' @param vcov list of variance-covariance matrices for the BLUEs
#' @param geno object of \code{\link{class_geno}} from \code{\link{read_geno}}
#' @param vars object of \code{\link{class_var}} from \code{\link{Stage2}}
#' @param mask (optional) data frame with possible columns "id","env","loc","trait"
#' @param method (optional) "MME", "Vinv", NULL (default). see Details
#' 
#' @return Object of \code{\link{class_prep}}
#' 
#' @import methods
#' @import Matrix
#' @importFrom stats model.matrix
#' @importFrom MASS ginv
#' @export

blup_prep <- function(data,vcov=NULL,geno=NULL,vars,mask=NULL,method=NULL) {
  
  stopifnot(is(data,"data.frame"))
  stopifnot(is(vars,"class_var"))
  data$id <- as.character(data$id)
  data$env <- as.character(data$env)
  if (!is.null(method)) {
    method <- toupper(method)
    stopifnot(method %in% c("MME","VINV"))
  }

  if (vars@model > 0)
    stopifnot(inherits(geno,"class_geno"))
  if (vars@model==3L)
    stopifnot(is(geno,"class_genoD"))
  
  n.trait <- 1
  if ("trait" %in% colnames(data)) {
    traits <- rownames(vars@resid)
    n.trait <- length(traits)
    data$trait <- as.character(data$trait)
    stopifnot(n.trait > 1)
  } 
  
  #if (length(vars@diagG)>0) {
  #  stopifnot(!is.null(geno))
  #}
  
  data <- data[!is.na(data$BLUE),]
  
  if (!is.null(mask)) {
    tmp <- intersect(colnames(mask),c("id","env","trait","loc"))
    stopifnot(length(tmp)>0)
    if (length(tmp) > 1) {
      tmp2 <- apply(mask[,tmp],1,paste,collapse=":")
      ix <- which(apply(data[,tmp],1,paste,collapse=":") %in% tmp2)
    } else {
      tmp2 <- mask[,tmp]
      ix <- which(data[,tmp] %in% tmp2)
    }
    if (length(ix) > 0) {
      data <- data[-ix,]
    }
  }
  
  if (!is.null(geno)) {
    stopifnot(inherits(geno,"class_geno"))
    #stopifnot(length(vars@diagG)>0)
    stopifnot(vars@model > 0)
    id <- intersect(data$id,rownames(geno@G))
    data <- data[data$id %in% id,]
    id <- rownames(geno@G)
    ploidy <- geno@ploidy
  } else {
    id <- unique(data$id)
    ploidy <- 0L
  }
  
  if (n.trait > 1) {
    data <- data[order(data$env,data$id,data$trait),]
    tmp <- paste(data$id,data$trait,sep=":")
  } else {
    data <- data[order(data$env,data$id),]
    tmp <- data$id
  }
  if (!is.null(vcov)) {
    #stopifnot(nrow(vars@meanOmega) > 0)
    data$env <- factor(data$env,levels=names(vcov))
    omega.list <- mapply(FUN=function(Q,ix){
                        ix2 <- match(ix,rownames(Q))
                        as(Q[ix2,ix2,drop=FALSE],"dpoMatrix")
                      },Q=vcov,ix=split(tmp,data$env))
  } else {
    data$env <- factor(data$env)
    omega.list <- lapply(split(tmp,data$env),function(rnames){
                                n <- length(rnames)
                                Q <- Matrix(0,nrow=n,ncol=n,dimnames=list(rnames,rnames))
                                return(Q)})
  }
  
  #redo envs because some may have been dropped
  data$env <- as.character(data$env)
  envs <- unique(data$env)
  n.env <- length(envs)
  omega.list <- omega.list[envs]
  data$env <- factor(data$env,levels=envs)
  n.obs <- sapply(omega.list,nrow)
  
  n.loc <- 1
  if ("loc" %in% colnames(data)) {
    locations <- rownames(vars@resid)
    n.loc <- length(locations)
    data$loc <- as.character(data$loc)
    stopifnot(n.loc > 1)
    missing.loc <- setdiff(locations,as.character(data$loc))
  } else {
    missing.loc <- character(0)
  }
  
  #Rlist
  if (n.trait==1) {
    if (n.loc > 1) {
      tmp <- vars@resid[as.character(data$loc)[match(envs,as.character(data$env))],1]
    } else {
      tmp <- rep(vars@resid[1,1],n.env)
    }
    Rlist <- mapply(FUN=function(n,v){Diagonal(n=n)*v},n=as.list(n.obs),v=as.list(tmp))
  } else {
    #multi-trait
    
    tmp <- split(data$id,data$env)
    tmp2 <- lapply(tmp,function(id) {
      I.mat <- Diagonal(n=length(id))
      dimnames(I.mat) <- list(id,id)
      kron2(A=I.mat, A.inv=I.mat, B=vars@resid)$mat
    })
    Rlist <- mapply(function(Q,rnames){
                    Q[rnames,rnames]},
                    Q=tmp2,rnames=lapply(omega.list,rownames))
  }

  n.id <- length(id)
  data$id <- factor(data$id,levels=id)
  
  if (n.trait > 1)
    data$trait <- factor(data$trait,levels=traits)
  if (n.loc > 1)
    data$loc <- factor(data$loc,levels=locations)
  
  n.mark <- length(vars@fix.eff.marker)
  if (n.mark > 0) {
    dat2 <- data.frame(id=rownames(geno@coeff),as.matrix(geno@coeff[,vars@fix.eff.marker]))
    colnames(dat2) <- c("id",vars@fix.eff.marker)
    data <- merge(data,dat2,by="id")
  }
  
  I.mat <- Diagonal(n=n.id)
  dimnames(I.mat) <- list(id,id)

  #Gmat$mat is the var-cov matrix of the random effects, Gmat$inv its inverse
  Gvar <- function(B, which=c("I","G","D")) {
    switch(match.arg(which),
           I = kron2(A=I.mat, A.inv=I.mat, B=B),
           G = kron2(A=geno@G, A.inv=eigen_inverse(geno@eigen.G), B=B),
           D = kron2(A=geno@D, A.inv=eigen_inverse(geno@eigen.D), B=B))
  }

  if (n.trait==1) {
    if (n.env > 1) {
      X <- sparse.model.matrix(~env-1,data,sep = "__")
      colnames(X) <- sub("env__","",colnames(X),fixed=T)
    } else {
      X <- sparse.model.matrix(~1,data)
      colnames(X) <- envs
    }
    
    if (n.loc > 1) {
      Z <- sparse.model.matrix(~id:loc-1,data,sep="__")
      colnames(Z) <- sub("loc__","",colnames(Z),fixed=T)
      colnames(Z) <- sub("id__","",colnames(Z),fixed=T)
      loc.id <- data.frame(as.matrix(expand.grid(locations,id))[,c(2,1)])
      colnames(loc.id) <- c("id","loc")
      Znames <- apply(loc.id,1,paste,collapse=":")
      Z <- Z[,Znames] #loc within id
      
      if (is.null(geno)) {
        Gmat <- Gvar(B=vars@geno1, which="I")
      } else {
        if (vars@model==1L) {
          Gmat <- Gvar(B=vars@geno1, which="G")
        } else {
          Gmat1 <- Gvar(B=vars@geno1, which="G")
          if (vars@model==3L) {
            Gmat2 <- Gvar(B=vars@geno2, which="D")
          
            #add to X matrix
            dom.covariate <- Z %*% kronecker((geno@coeff.D/(geno@scale*(ploidy-1))) %*% 
                                               matrix(1,nrow=ncol(geno@coeff.D),ncol=1),diag(n.loc))
            colnames(dom.covariate) <- paste("heterosis",locations,sep=":")
            X <- cbind(X,dom.covariate)
          } else {
            Gmat2 <- Gvar(B=vars@geno2, which="I")
          }
          Gmat <- list(mat=bdiag2(Gmat1$mat,Gmat2$mat),inv=bdiag2(Gmat1$inv,Gmat2$inv))
        }
      } 
    } else {
      
      Z <- sparse.model.matrix(~id-1,data,sep="__")
      colnames(Z) <- sub("id__","",colnames(Z),fixed=T)

      if (is.null(geno)) {
        Gmat <- Gvar(B=vars@geno1, which="I")
      } else {
        if (vars@model==1L) {
          Gmat <- Gvar(B=vars@geno1, which="G")
        } else {
          Gmat1 <- Gvar(B=vars@geno1, which="G")
          if (vars@model==3L) {
            Gmat2 <- Gvar(B=vars@geno2, which="D")
            dom.covariate <- as.numeric(Z %*% (geno@coeff.D/(geno@scale*(ploidy-1))) %*% 
                                          matrix(1,nrow=ncol(geno@coeff.D),ncol=1))
            X <- cbind(X,heterosis=dom.covariate)
          } else {
            Gmat2 <- Gvar(B=vars@geno2, which="I")
          }
          Gmat <- list(mat=bdiag2(Gmat1$mat,Gmat2$mat),inv=bdiag2(Gmat1$inv,Gmat2$inv))
        }
      }
    }
  } else {
    #multi-trait 
    Z <- sparse.model.matrix(~id:trait-1,data,sep="__")
    colnames(Z) <- sub("trait__","",colnames(Z),fixed=T)
    colnames(Z) <- sub("id__","",colnames(Z),fixed=T)
    
    if (n.env > 1) {
      X <- sparse.model.matrix(~env:trait-1,data,sep="__")
      colnames(X) <- sub("trait__","",colnames(X),fixed=T)
      colnames(X) <- sub("env__","",colnames(X),fixed=T)
    } else {
      X <- sparse.model.matrix(~trait-1,data,sep="__")
      colnames(X) <- sub("trait__","",colnames(X),fixed=T)
      colnames(X) <- paste(envs,colnames(X),sep=":")
    }
    
    trait.id <- data.frame(as.matrix(expand.grid(traits,id))[,c(2,1)])
    colnames(trait.id) <- c("id","trait")
    Znames <- apply(trait.id,1,paste,collapse=":")
    Z <- Z[,Znames] #trait within id

    if (is.null(geno)) {
      Gmat <- Gvar(B=vars@geno1, which="I")
    } else {
      if (vars@model==1L) {
        Gmat <- Gvar(B=vars@geno1, which="G")
      } else {
        Gmat1 <- Gvar(B=vars@geno1, which="G")
        if (vars@model==3L) {
          Gmat2 <- Gvar(B=vars@geno2, which="D")
        
          #add to X matrix
          dom.covariate <- Z %*% kronecker((geno@coeff.D/(geno@scale*(ploidy-1))) %*%
                                             matrix(1,nrow=ncol(geno@coeff.D),ncol=1),diag(n.trait))
          colnames(dom.covariate) <- paste("heterosis",traits,sep=":")
          X <- cbind(X,dom.covariate)
        } else {
          Gmat2 <- Gvar(B=vars@geno2, which="I")
        }
        Gmat <- list(mat=bdiag2(Gmat1$mat,Gmat2$mat),inv=bdiag2(Gmat1$inv,Gmat2$inv))
      }
    } 
  }
  
  if (vars@model > 1L) {
    Z <- cbind(Z,Z)
  }

  if (n.mark > 0) {
    if (n.trait > 1) {
      q <- paste(vars@fix.eff.marker,"trait",sep=":")
    } else {
      if (n.loc > 1) {
        q <- paste(vars@fix.eff.marker,"loc",sep=":")
      } else {
        q <- vars@fix.eff.marker
      }
    }
    if (n.mark > 1) {
      q <- paste(q,collapse="+")
    } 
    q <- paste0("~",q)
    q <- paste0(q,"-1")
    tmp <- sparse.model.matrix(formula(q),data,sep="__")
    if (n.trait > 1) {
      colnames(tmp) <- sub("trait__","",colnames(tmp),fixed=T)
      tmp2 <- expand.grid(factor(vars@fix.eff.marker,levels=vars@fix.eff.marker,ordered=T),traits)
      tmp2 <- tmp2[order(tmp2$Var1),]
      marker.trait <- apply(tmp2,1,paste,collapse=":")
      colnames(tmp) <- marker.trait
    }
    if (n.loc > 1) {
      colnames(tmp) <- sub("loc__","",colnames(tmp),fixed=T)
      tmp2 <- expand.grid(factor(vars@fix.eff.marker,levels=vars@fix.eff.marker,ordered=T),locations)
      tmp2 <- tmp2[order(tmp2$Var1),]
      marker.loc <- apply(tmp2,1,paste,collapse=":")
      colnames(tmp) <- marker.loc
    }
    X <- cbind(X,tmp)
  }
  
  covariates <- unlist(attributes(vars@fix.eff.marker))
  n.covar <- length(covariates)
  if (n.covar > 0) {
    q <- paste(covariates,collapse="+")
    q <- paste0("~",q)
    q <- paste0(q,"-1")
    tmp <- sparse.model.matrix(formula(q),data,sep="__")
    X <- cbind(X,tmp)
  } 
  n.fix <- ncol(X)
  m <- ncol(Z)
  n <- nrow(Z)
  var.u <- Gmat$mat

  if (is.null(method)) {
    if (m < n) {
      method <- "MME"
    } else {
      method <- "VINV"
    }
  }
  
  if (method=="MME") {
    #Construct MME coefficient matrix.
    #Only crossprod(Rinv) = solve(omega+R) is needed, so factor once and invert the
    #triangle, rather than forming the inverse and factoring it again.
    tmp <- mapply(function(a,b){
                    U <- chol(as.matrix(a+b))
                    t(backsolve(U, diag(nrow(U))))
                  },omega.list,Rlist,SIMPLIFY=FALSE)
    Rinv <- bdiag(tmp)

    RZ <- Rinv %*% Z
    RX <- Rinv %*% X
    Q <- cbind(RX, RZ)

    random.ix <- (n.fix+1):(n.fix+m)

    #When geno is supplied, Gmat$inv is the dense G-inverse and the random-effect block
    #of MME is dense, so sparse storage/factorization is all overhead: use LAPACK.
    #Without markers Gmat$inv really is sparse, and the sparse path stays faster.
    if (nnzero(Gmat$inv)/prod(dim(Gmat$inv)) > 0.2) {
      MME <- as.matrix(crossprod(Q))
      MME[random.ix,random.ix] <- MME[random.ix,random.ix] + as.matrix(Gmat$inv)
      MME.inv <- chol2inv(chol(MME))
    } else {
      MME <- as(crossprod(Q) + bdiag(Matrix(0,nrow=n.fix,ncol=n.fix),Gmat$inv),"symmetricMatrix")
      MME.inv <- as(solve(MME),"symmetricMatrix")
    }

    soln <- as.numeric(MME.inv %*% as.numeric(crossprod(Q,Rinv%*%data$BLUE)))
    fixed <- soln[1:n.fix]
    names(fixed) <- colnames(X)
    random <- soln[random.ix]
    var.uhat <- var.u - MME.inv[random.ix,random.ix]
    var.bhat <- forceSymmetric(Matrix(MME.inv[(1:n.fix),(1:n.fix),drop=FALSE]))
    dimnames(var.bhat) <- list(names(fixed),names(fixed))
    cov.buhat <- Matrix(-MME.inv[(1:n.fix),random.ix,drop=FALSE])

  } else {
    #invert V
    tmp <- mapply(function(a,b){as(a+b,"dpoMatrix")},omega.list,Rlist)
    Rmat <- bdiag(tmp)

    Zvar.u <- Z %*% var.u
    Vinv <- as(solve(as(tcrossprod(Zvar.u, Z) + Rmat,"symmetricMatrix")),"symmetricMatrix")

    #P = Vinv - VX %*% var.bhat %*% t(VX) is the projection matrix. Forming it from
    #VX (n x n.fix) avoids the O(n^3) Cholesky factorizations of V and of I - HX.
    VX <- Vinv %*% X
    tmp <- crossprod(X, VX)
    tmp2 <- try(solve(tmp),silent=TRUE)
    if (is(tmp2,"try-error"))
      tmp2 <- MASS::ginv(as.matrix(tmp))
    var.bhat <- forceSymmetric(tmp2)
    fixed <- as.numeric(var.bhat %*% crossprod(VX, data$BLUE))
    names(fixed) <- colnames(X)
    dimnames(var.bhat) <- list(names(fixed),names(fixed))

    VXtZG <- crossprod(VX, Zvar.u)                          #n.fix x m
    PZvar.u <- Vinv %*% Zvar.u - VX %*% (var.bhat %*% VXtZG)
    random <- as.numeric(crossprod(PZvar.u, data$BLUE))
    var.uhat <- forceSymmetric(crossprod(Zvar.u, PZvar.u))
    cov.buhat <- var.bhat %*% VXtZG
  }
  
  fixed.marker <- numeric(0)
  heterosis <- numeric(0)
  
  if (n.loc > 1) {
    loc.env <- unique(data[,c("loc","env")])
    loc.env <- loc.env[order(loc.env$loc,loc.env$env),]
    ix <- match(names(fixed)[1:n.env],as.character(loc.env$env))
    avg.env <- tapply(fixed[1:n.env],loc.env$loc[ix],mean)
    tmp <- as.numeric(avg.env)
    names(tmp) <- names(avg.env)
    avg.env <- tmp
  } else {
    if (n.trait > 1) {
      env.trait <- expand.grid(env=envs,trait=traits)
      env.trait$et <- apply(env.trait[,1:2],1,paste,collapse=":")
      ix <- match(env.trait$et,names(fixed))
      avg.env <- tapply(fixed[ix],env.trait$trait,mean)
      tmp <- as.numeric(avg.env)
      names(tmp) <- names(avg.env)
      avg.env <- tmp
      n.env <- n.env*n.trait
    } else {
      avg.env <- mean(fixed[1:n.env])
    }
  }
  
  if (n.covar > 0) {
    #evaluate covariates at mean value
    avg.env <- avg.env + sum(apply(data[,covariates,drop=FALSE],2,mean)*fixed[(n.fix-n.covar+1):n.fix])
  }
  
  nlt <- max(n.loc,n.trait)
  if (vars@model==3L) {
    heterosis <- fixed[n.env+1:nlt]
    if (length(missing.loc)>0)
      heterosis[paste("heterosis",missing.loc,sep=":")] <- mean(heterosis,na.rm=T)
    if (n.mark > 0) 
      fixed.marker <- fixed[n.env+nlt+1:(n.mark*nlt)]
  } else {
    if (n.mark > 0) 
      fixed.marker <- fixed[n.env+1:(n.mark*nlt)]
  }
  if (length(missing.loc)>0) {
    avg.env[missing.loc] <- mean(avg.env,na.rm=T)
  }
  
  new(Class="class_prep",id=id,ploidy=ploidy,var.u=var.u,var.uhat=var.uhat,
      var.bhat=var.bhat, cov.buhat=cov.buhat,
      avg.env=avg.env,heterosis=heterosis,fixed.marker=fixed.marker,B=vars@B,
      random=random, geno1.var=vars@geno1, geno2.var=vars@geno2, 
      model=vars@model)
  
}    
