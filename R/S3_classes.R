#' Print an object of class \code{Stage1}
#'
#' @description Prints information about \code{Stage1} function.
#'
#' @aliases print.Stage1
#' @usage \method{print}{Stage1}(x, ...)
#' @param x An object fitted with the function \code{Stage1()}.
#' @param ... Options used by the tibble package to format the output.
#' @method print Stage1
#' @return an object inheriting from class \code{Stage1}.
#' @importFrom utils head
#' @export
#' @examples
#' \dontrun{
#' model <- Stage1(
#'   filename = "pheno.csv",
#'   traits = "total.yield",
#'   solver = "spats", 
#'   spline = c("row", "range")
#' )
#' print(model)
#' }
print.Stage1 <- function(x, ...) {
  cat("\nObject of class 'Stage1'\n")
  # Components
  cat("\nComponents:", paste(names(x), collapse = ", "), "\n")
  dt_blues <- x$blues
  # Summary
  n_ind <- length(unique(dt_blues$id))
  n_env <- length(unique(dt_blues$env))
  n_loc <- if ("loc" %in% names(dt_blues)) {
    length(unique(dt_blues$loc))
  } else {
    1L
  }
  n_traits <- if ("trait" %in% names(dt_blues)) {
    length(unique(dt_blues$trait))
  } else {
    1L
  }
  cat("\nSummary:\n")
  cat("  Environments :", n_env, "\n")
  cat("  Locations    :", n_loc, "\n")
  cat("  Individuals  :", n_ind, "\n")
  cat("  Traits       :", n_traits, "\n")
  # BLUEs
  cat("\nBLUEs (first 3 rows):\n")
  print(
    as.data.frame(utils::head(dt_blues, 3)),
    digits = 3,
    row.names = FALSE,
    ...
  )
  # Environmental connectivity
  conn <- table(unique(dt_blues[c("id", "env")]))
  connectivity <- crossprod(conn)
  names(dimnames(connectivity)) <- NULL
  if (n_env <= 10) {
    cat("\nConnectivity among environments:\n")
    print(connectivity)
  }
  cat("\n")
  invisible(x)
}

#' Show an object of class \code{class_genoD}
#'
#' @param object An object of class \code{class_genoD}.
#'
#' @return Invisibly returns \code{object}.
#'
#' @method show class_genoD
#' @export
setMethod("show", "class_genoD", function(object) {
  n_ind <- nrow(object@coeff)
  n_markers <- ncol(object@coeff)
  cat("<class_genoD>\n\n")
  cat("Genotypic data:\n")
  cat("  Individuals :", n_ind, "\n")
  cat("  Markers     :", n_markers, "\n")
  cat("  Ploidy      :", object@ploidy, "\n")
  cat("\nRelationship matrices:\n")
  cat("  Additive    :", paste(dim(object@G), collapse = " x "), "\n")
  cat("  Dominance   :", paste(dim(object@D), collapse = " x "), "\n")
  if (!is.null(object@map)) {
    cat("\nMap: available\n\n")
  } else {
    cat("\nMap: not available\n\n")
  }
  invisible(object)
})

#' Show an object of class \code{class_geno}
#'
#' @param object An object of class \code{class_geno}.
#'
#' @return Invisibly returns \code{object}.
#'
#' @method show class_geno
#' @export
setMethod("show", "class_geno", function(object) {
  n_ind <- nrow(object@coeff)
  n_markers <- ncol(object@coeff)
  cat("<class_geno>\n\n")
  cat("Genotypic data:\n")
  cat("  Individuals :", n_ind, "\n")
  cat("  Markers     :", n_markers, "\n")
  cat("  Ploidy      :", object@ploidy, "\n")
  cat("\nRelationship matrix:\n")
  cat("  Additive    :", paste(dim(object@G), collapse = " x "), "\n")
  if (!is.null(object@map)) {
    cat("\nMap: available\n\n")
  } else {
    cat("\nMap: not available\n\n")
  }
  invisible(object)
})
