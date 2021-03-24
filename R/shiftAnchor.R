#' Define helper function for flexibly shifting an anchor
#'
#' @param a GRanges object
#' @param p Position within anchors to resize the bin. Can be a character or integer vector of length 1 or length(a) designating the position for each element in bedpe. Character options are "start", "end" and "center". Integers are referenced from the start position for '+' and '*' strands and from the end position for the '-' strand.
#'
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
#' ## Shifting anchors by keyword
#' shiftAnchor(gr1, 'start')
#' shiftAnchor(gr1, 'end')
#' shiftAnchor(gr1, 'center')
#' # shiftAnchor(gr1, 'blah') error
#'
#' ## Shifting anchors by position
#' shiftAnchor(gr1, 100)
#' shiftAnchor(gr1, c(100, 200, 300))
#' # shiftAnchor(gr1, c(100, 200, 300, 400)) error
#' # shiftAnchor(gr1, c(100, 200)) error
#'
#' ## Shifting back to TSS
#' shiftAnchor(gr2, 2000)
#'
#'
shiftAnchor <- function(a, p) {

  ## Shift anchors appropriately
  if (class(p) == "character") {

    if (p %in% c('start', 'end', 'center')) {

      a %<>% resize(width = 1, fix = p)

    } else {

      stop('p character must be one of "start", "end" or "center"', call. = T)

    }

  } else if (class(p) %in% c('numeric', 'integer')) {

    ## Pull out strand of anchor (as vector)
    sa <- as.vector(strand(a))

    ## Subset for strand if length(p) > 1
    if (length(p) > 1) {
      stopifnot(length(p) == length(a))
      pp <- p[which(sa %in% c('+', '*'))]
      pm <- p[which(sa == '-')] * -1
    } else {
      pp <- p
      pm <- p * -1
    }

    ## Shift '+' or '*' strand
    a[sa %in% c('+', '*')] %<>%
      resize(width = 1, fix = 'start') %>%
      shift(shift = pp)

    ## Shift '-' strand
    a[sa == '-'] %<>%
      resize(width = 1, fix = 'start') %>%
      shift(shift = pm)

  } else {

    stop('class(p) must be either character or numeric.', call. = T)

  }

  ## Return anchor
  return(a)

}
