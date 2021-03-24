#' Define function to flexibly bin bedpe data by hic resolution
#'
#' @param bedpe GInteractions or data.table object with paired interactions
#' @param res Integer - resolution in which to bin bedpe anchors
#' @param a1Pos,a2Pos Position within anchors to resize the bin. Can be a character or integer vector of length 1 or length(bedpe) designating the position for each element in bedpe. Character options are "start", "end" and "center". Integers are referenced from the start position for '+' and '*' strands and from the end position for the '-' strand.
#'
#'
#' @export
#'
binBedpe <- function(bedpe, res, a1Pos, a2Pos) {

  ## Extract anchors
  a1 <- anchors(bedpe, type = "first")
  a2 <- anchors(bedpe, type = "second")

  ## Bin anchors
  a1 <- binAnchor(a = a1, p = a1Pos, res = res)
  a2 <- binAnchor(a = a2, p = a2Pos, res = res)

  ## Binned GInteractions object
  gi <- GInteractions(a1, a2)

  ## Add back metadata
  mcols(gi) <- mcols(bedpe)

  ## Return binned bedpe
  return(gi)

}
