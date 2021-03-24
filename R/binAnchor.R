#' Define helper function for binning an anchor
#'
#' @param a GRanges object
#' @param p Position within anchors to resize the bin. Can be a character or integer vector of length 1 or length(a) designating the position for each element in bedpe. Character options are "start", "end" and "center". Integers are referenced from the start position for '+' and '*' strands and from the end position for the '-' strand.
#' @param res Integer - resolution in which to bin the anchor
#'
#' @examples
#'
#' ## Create example GRanges
#' gr1 <- GRanges(seqnames = "chr1",
#'                ranges = IRanges(start = rep(5000,3), end = rep(6000,3)),
#'                strand = c('+', '-', '*'))
#'
#' gr2 <- gr1 %>% promoters(upstream = 2000, downstream = 200)
#'
#' ## Binning the results
#' binAnchor(gr1, 'start', 1000)
#' binAnchor(gr1, 'end', 1000)
#' binAnchor(gr1, 'center', 1000)
#'
#' ## Bin after shifting back to TSS
#' binAnchor(gr2, 2000, 1000)
#'
#'
binAnchor <- function(a, p, res) {

  ## Shift, bin, and trim anchors
  a %<>%
    shiftAnchor(p) %>%
    mutate(start = floor(start/res)*res,
           end = floor(start/res)*res + res) %>%
    trim() %>%
    suppressWarnings()

  return(a)

}
