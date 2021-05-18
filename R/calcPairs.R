#' Define function to calculate all paired interactions within a windowSize
#'
#' @description
#' calcPairs will calculate pairs of interactions within a rolling windowSize.
#' By default, pairs of interactions are returned as StrictGInteractions where
#' anchor1 <= anchor2 for all interactions (mode = 'strict'). For more information
#' see (?InteractionSet::`GInteractions-class`). Self interactions (e.g. 1, 1 or 2, 2)
#' and order (e.g. 1, 2 is the same as 2, 1) are ignored.
#'
#' @param gr GRanges object
#' @param windowSize integer defining the window
#' @param mode If mode="strict", a StrictGInterctions object is returned with anchor
#'  indices swapped such that anchor1 <= anchor2 for all interactions. If mode="reverse",
#'  a ReverseStrictGInterctions object is returned with anchor indices swapped such that
#'  anchor1 >= anchor2.
#'
#' @return Returns a StrictGInterations object
#'
#' @examples
#' ## Load TxDb
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' ## Generate 50Kb bins across hg19
#' library(GenomicRanges)
#' bins <- tileGenome(seqinfo(txdb), tilewidth = 10e3, cut.last.tile.in.chrom = TRUE)
#'
#' ## Calculate all 50Kb bin pairs within 1Mb
#' calcPairs(gr = bins, windowSize = 1e6)
#'
#' @export
#'
calcPairs <- function(gr, windowSize, mode = 'strict') {

  ## Check arguments ---------------------------------------------------------------------

  stopifnot(isClass('GRanges', gr))
  stopifnot(isClass('integer', windowSize))
  stopifnot(length(windowSize) == 1L)
  mode <- match.arg(mode, choices = c('strict', 'reverse'))

  ## Begin processing --------------------------------------------------------------------

  ## Sort gr
  gr <- gr %>% sort()

  ## Define windows by windowSize
  windows <-
    gr %>%
    resize(width = windowSize, fix = 'start') %>%
    suppressWarnings() %>%
    trim()

  ## Group gr by windows
  ov <-
    findOverlaps(gr, windows, type = 'within') %>%
    as.data.table() %>%
    `colnames<-`(c('gr', 'windows'))

  ## Set up progress bar
  mx <- round(uniqueN(ov$windows)*1.01, 0)
  pb <- txtProgressBar(min = 1, max = mx, initial = NA, style = 3)

  ## Iterate over combinations
  ov <- ov[, {setTxtProgressBar(pb, .GRP); .(gr[1], gr[-1])}, by = windows]

  ## Remove NA's (last value for each chromosome)
  ov <- na.omit(ov)

  ## Convert coordinates to GInteractions object
  gi <- GInteractions(anchor1 = ov$V1,
                      anchor2 = ov$V2,
                      regions = gr,
                      mode = mode)

  ## Finish progress bar
  setTxtProgressBar(pb, value = mx)
  close(pb)
  cat('\n')

  return(gi)
}
