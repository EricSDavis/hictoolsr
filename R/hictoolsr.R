#' hictoolsr: a set of tools for Hi-C data analysis
#'
#'
#'
#' @rawNamespace import(data.table, except = c(between, shift, first, second))
#' @rawNamespace import(grid)
#' @rawNamespace import(BentoBox)
#' @rawNamespace import(InteractionSet)
#' @rawNamespace import(plyranges)
#' @importFrom strawr straw
#' @importFrom tidyr complete full_seq
#' @importFrom reshape2 acast
#' @importFrom magrittr %<>%
#' @importFrom GenomicRanges shift resize trim makeGRangesFromDataFrame
#' @importFrom S4Vectors mcols `mcols<-` first second DataFrame
#'
"_PACKAGE"
