#' Filters out short BEDPE interactions
#'
#' Interactions from a BEDPE file that would result in APA plots intersecting the diagonal
#' (given a resolution and buffer size) are removed from the BEDPE file.
#'
#' @param bedpe data.frame or data.table in BEDPE format
#' @param res integer resolution of bedpe bins
#' @param buffer integer number of res-length bins from the center pixel
#'
#' @return Returns a filtered BEDPE
#'
#' @export
#'
filterBedpe <- function(bedpe, res, buffer) {

  ## TODO Generalize this function for inter & intra pairs

  ## Only valid for intrachromosomal bedpe files
  if(all(bedpe[[1]] != bedpe[[4]])) {
    stop("bb_filterBedpe is only valid for intrachromosomal pairs")
  }

  ## Distance from diagonal to center bedpe pixel
  d <- (sqrt((abs(bedpe[[5]] - bedpe[[2]])^2 * 2)) / 2)

  ## Distance from center bedpe pixel to corner of apa
  y <- sqrt((res*(buffer*2 + 1))^2 * 2) / 2

  ## Return loops that don't intersect the diagonal
  return(bedpe[d > y])

}
