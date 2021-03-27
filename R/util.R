#' Helper function to check that a GInteractions object has been binned
#'
#' @param x GInteractions object.
#' @param res Integer - resolution in base-pairs.
#'
#' @rdname check_binned
#'

check_binned <- function(x, res) {

  ## Check that each anchor is binned correctly
  a1 <- length(which(width(first(x)) != res+1)) == 0
  a2 <- length(which(width(second(x)) != res+1)) == 0

  return(a1 & a2)

}
