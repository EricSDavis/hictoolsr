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

  ## Input checking and processing -------------------------------------------------------

  ## Convert to GInteractions object
  if (!"GInteractions" %in% class(bedpe)) {
    warning('class(bedpe) != "GInteractions". Using as_ginteractions() to convert.')
    bedpe <- as_ginteractions(bedpe)
  }

  ## Check that each anchor is binned correctly
  binned <- check_binned(bedpe, res)

  ## Bin if not correctly binned and give warning
  if (!binned) {
    warning(strwrap("bedpe is not binned correctly.
        It is a good idea to bin bedpe before filtering.
        Binning each anchor at center position.
        For more options use the binBedpe() function."),
            immediate. = TRUE,
            call. = FALSE)

    bedpe <- binBedpe(bedpe, res, a1Pos = 'center', a2Pos = 'center')
  }

  ## Convert to data.table in bedpe format
  bedpe <-
    as.data.table(bedpe)[,c(1:3, 6:8)] %>%
    `colnames<-`(c("chr1", "start1", "end1", "chr2", "start2", "end2"))


  ## Calculate and filter by distances ---------------------------------------------------

  ## Distance from center bedpe pixel to corner of apa
  y <- sqrt((res*(buffer*2 + 1))^2 * 2) / 2

  ## For interchromosomal reads:
  ## Distance from diagonal to center bedpe pixel
  bedpe[chr1 == chr2, d := (sqrt((abs(start2 - start1)^2 * 2)) / 2)]

  ## Return interactions that don't intersect the diagonal
  bedpe <- bedpe[is.na(d) | d > y, -c('d')]

  ## Return GInteractions object
  return(as_ginteractions(bedpe))

}
