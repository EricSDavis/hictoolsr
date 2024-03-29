#' Convert DataFrames to GInteraction objects
#'
#' makeGInteractionsFromDataFrame takes a paired-interaction (i.e. BEDPE) formatted data-frame-like object and converts it to a GInteractions object. For convenience, \code{as_ginteractions()} can be used as an alias.
#'
#' @param df A data.table, data.frame, or DataFrame object. Assumes that the first 6   colummns are in the format chr1, start1, end1 and chr2, start2, end2, representing each pair of interactions.
#' @param keep.extra.columns TRUE or FALSE (the default). If TRUE, the columns in df that are not used to form the genomic ranges of the returned GRanges object are then returned as metadata columns on the object. Otherwise, they are ignored. If df has a width column, then it's always ignored.
#' @param starts.in.df.are.0based TRUE or FALSE (the default). If TRUE, then the start positions of the genomic ranges in df are considered to be 0-based and are converted to 1-based in the returned GRanges object. This feature is intended to make it more convenient to handle input that contains data obtained from resources using the "0-based start" convention. A notorious example of such resource is the UCSC Table Browser (http://genome.ucsc.edu/cgi-bin/hgTables).
#'
#' @return GInteraction object
#'
#' @rdname makeGInteractionsFromDataFrame
#' @aliases as_ginteractions
#'
#' @examples
#' ## data.frame
#' df <- data.frame(chr1 = "chr1", x1 = 10000, x2 = 20000,
#'                  chr2 = "chr1", y1 = 30000, y2 = 40000)
#' makeGInteractionsFromDataFrame(df)
#'
#' ## data.table
#' df <- data.table(chr1 = "chr1", x1 = 10000, x2 = 20000,
#'                  chr2 = "chr1", y1 = 30000, y2 = 40000)
#' makeGInteractionsFromDataFrame(df)
#'
#' ## DataFrame
#' df <- DataFrame(chr1 = "chr1", x1 = 10000, x2 = 20000,
#'                 chr2 = "chr1", y1 = 30000, y2 = 40000)
#' makeGInteractionsFromDataFrame(df)
#'
#' ## Alias
#' df <- data.frame(chr1 = "chr1", x1 = 10000, x2 = 20000,
#'                  chr2 = "chr1", y1 = 30000, y2 = 40000,
#'                  pval = 0.05, dist = 10000)
#' as_ginteractions(df)
#'
#' ## Additional metadata
#' df <- data.frame(chr1 = "chr1", x1 = 10000, x2 = 20000,
#'                  chr2 = "chr1", y1 = 30000, y2 = 40000,
#'                  pval = 0.05, dist = 10000)
#' as_ginteractions(df)
#'
#' ## Remove additional metadata
#' as_ginteractions(df, keep.extra.columns = F)
#'
#' ## Add 1 to starts (for 0-based programs)
#' as_ginteractions(df, starts.in.df.are.0based = T)
#'
#' ## Bad usage
#' df <- data.frame(chr1 = "chr1", x1 = 10000, x2 = 20000,
#'                  chr2 = "chr1", y1 = 30000)
#' makeGInteractionsFromDataFrame(df)
#'
#'
#' @export
#'
#'
makeGInteractionsFromDataFrame <- function(df,
                                           keep.extra.columns = TRUE,
                                           starts.in.df.are.0based = FALSE) {

  ## Convert data.table/data.frame to DataFrame
  if ("data.frame" %in% class(df)) {
    df <- DataFrame(df)
  } else if ("DFrame" %in% class(df)) {
    df <- df
  } else {
    stop("class(df) must be either 'data.frame', 'data.table', or 'DFrame'.")
  }

  ## Handle improper dimensions
  if(ncol(df) < 6) {
    stop("ncol(df) must be >= 6 and start with paired interactions (i.e. chr1, start1, end1 and chr2, start2, end2).")
  }

  ## Split into anchors
  a1 <- df[1:3] %>% `colnames<-`(c('seqnames', 'start', 'end'))
  a2 <- df[4:6] %>% `colnames<-`(c('seqnames', 'start', 'end'))

  ## Convert anchors to GRanges
  a1 <- makeGRangesFromDataFrame(a1, starts.in.df.are.0based = starts.in.df.are.0based)
  a2 <- makeGRangesFromDataFrame(a2, starts.in.df.are.0based = starts.in.df.are.0based)

  ## Create GInteractions object
  gi <- GInteractions(a1, a2)

  ## Add in metadata columns
  if (keep.extra.columns & ncol(df) > 6) {
    mcols(gi) <- df[7:ncol(df)]
  }

  ## Return
  return(gi)

}


#' @rdname makeGInteractionsFromDataFrame
#' @export
as_ginteractions <- makeGInteractionsFromDataFrame
