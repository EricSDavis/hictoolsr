#' Extract BEDPE interactions from a Hi-C matrix
#'
#' @param bedpe data.frame or data.table in BEDPE format.
#' @param hic string path to .hic file.
#' @param norm string hic normalization <NONE/VC/VC_SQRT/KR>.
#' @param res integer resolution of bedpe bins.
#' @param buffer integer number of res-length bins from the center pixel.
#' @param filter TRUE or FALSE (default). If TRUE, \code{filterBedpe()} will be used to remove short intrachromosomal interactions that would otherwise cross the diagonal.
#'
#' @return Returns a list APA matricies that can be compiled and plotted.
#'
#'
#' @export
#'
calcApa <- function(bedpe, hic, norm = 'NONE', res = 10000, buffer = 5, filter = FALSE, out = 'observed') {

  ## Set scipen to 999 for straw extraction
  oo <- options()
  options(scipen = 999)
  on.exit(options(oo))

  ## Input checking and processing -------------------------------------------------------

  ## Convert to GInteractions object
  if ("GInteractions" != class(bedpe)) {
    warning('class(bedpe) != "GInteractions". Using as_ginteractions() to convert.')
    bedpe <- as_ginteractions(bedpe)
  }

  ## Check that each anchor is binned correctly
  binned <- check_binned(bedpe, res)

  ## Bin if not correctly binned and give warning
  if (!binned) {
    warning(strwrap("bedpe is not binned correctly.
        Binning each anchor at center position.
        For more options use the binBedpe() function."),
            immediate. = TRUE,
            call. = FALSE)

    bedpe <- binBedpe(bedpe, res, a1Pos = 'center', a2Pos = 'center')
  }

  ## Filter bedpe
  if (filter) {
    bedpe <- filterBedpe(bedpe, res, buffer)
  }

  ## Convert to data.table in reduced format
  bedpe <-
    as.data.table(bedpe)[,c(1:2,6:7)] %>%
    `colnames<-`(c("chr1", "start1", "chr2", "start2"))

  ## Remove 'chr' for straw (also converts to character vector)
  bedpe$chr1 <- gsub('chr', '', bedpe$chr1)
  bedpe$chr2 <- gsub('chr', '', bedpe$chr2)

  ## Put data in correct order for straw extraction --------------------------------------

  ## Intrachromosomal: Flip column order where start1 > start2
  bedpe[chr1 == chr2 & start1 > start2, `:=`(start1=start2, start2=start1)]


  ## Interchromosomal: Flip column order where chr1 > chr2

  ## Replace X, Y with 23, 24
  bedpe[chr1 == 'X', chr1 := 23]
  bedpe[chr1 == 'Y', chr1 := 24]
  bedpe[chr2 == 'X', chr2 := 23]
  bedpe[chr2 == 'Y', chr2 := 24]

  ## Coerce chromsome strings to numeric
  bedpe[,chr1 := as.numeric(chr1)]
  bedpe[,chr2 := as.numeric(chr2)]

  ## Flip column order where chr1 > chr2
  bedpe[chr1 > chr2, `:=`(chr1=chr2, start1=start2, chr2=chr1, start2=start1)]

  ## Set order
  setorderv(bedpe, c('chr1', 'chr2', 'start1', 'start2'))

  ## Coerce chromosome numeric to strings
  bedpe[,chr1 := as.character(chr1)]
  bedpe[,chr2 := as.character(chr2)]

  ## Replace 23, 24 with X, Y
  bedpe[chr1 == '23', chr1 := 'X']
  bedpe[chr1 == '24', chr1 := 'Y']
  bedpe[chr2 == '23', chr2 := 'X']
  bedpe[chr2 == '24', chr2 := 'Y']

  ## Prepare groups and locs for extraction ----------------------------------------------

  ## Get min and max ranges to extract for each chr pair
  regs <- bedpe[, .(min = min(start1, start2) - res*buffer,
                    max = max(start1, start2) + res*buffer),
                by = .(chr1, chr2)]

  ## Paste ranges to get chr1 and chr2 locs
  locs <- regs[, .(chr1loc = paste(chr1, min, max, sep = ':'),
                   chr2loc = paste(chr2, min, max, sep = ':'))]

  ## Add condition to check and correct "chr" designation
  ## in the .hic file
  hicChroms <- strawr::readHicChroms(hic)$name
  if (length(grep("chr", hicChroms)) != 0) {
    locs$chr1loc <- paste0("chr", locs$chr1loc)
    locs$chr2loc <- paste0("chr", locs$chr2loc)
  }

  ## Get group information from chr combinations
  bedpe[, GRP := .GRP, by = .(chr1, chr2)]


  ## Extract slices from hic file --------------------------------------------------------

  ## Start progress bar
  pb <- progress::progress_bar$new(
    format = "  :step [:bar] :percent elapsed: :elapsedfull",
    clear = F, total = nrow(locs)*3+1)
  pb$tick(0)

  ## Initialize list for storing lists of wide matricies
  m <- list()

  ## Begin extraction and searching for counts
  for (i in 1:nrow(locs)) {

    ## Update progress
    pb$tick(tokens=list(step=sprintf('Extracting chr%s-chr%s',
                                     regs$chr1[i], regs$chr2[i])))

    ## Extract sparse matrix from group
    sparseMat <- as.data.table(strawr::straw(norm = norm,
                                             fname = hic,
                                             chr1loc = locs$chr1loc[i],
                                             chr2loc = locs$chr2loc[i],
                                             unit = "BP",
                                             binsize = res,
                                             matrix = out))

    ## Update progress
    pb$tick(tokens=list(step=sprintf('Matching chr%s-chr%s',
                                     regs$chr1[i], regs$chr2[i])))

    ## Subset interactions by group
    g <- bedpe[GRP == i,]

    ## Generate bins for rows and columns
    bins <- g[,.(x = seq(from = start1 - res*buffer, to = start1 + res*buffer, by = res),
                 y = seq(from = start2 - res*buffer, to = start2 + res*buffer, by = res),
                 counts = 0), by = .(groupRow = 1:nrow(g))]

    ## Expand bins to long format (fast cross-join)
    longMat <- bins[, do.call(CJ, c(.SD, sorted = F)),
                    .SDcols = c('x', 'y'), by = groupRow]
    longMat$counts <- 0

    ## Rename columns and rearrange
    longMat <- longMat[,.(x, y, counts, groupRow)]

    ## Set keys
    setkeyv(sparseMat, c('x', 'y'))

    ## Get counts by key
    longMat$counts <- sparseMat[longMat]$counts

    ## Set unmatched counts (NA) to 0
    longMat[is.na(counts), counts := 0]

    ## Update progress
    pb$tick(tokens=list(step=sprintf('Aggregating chr%s-chr%s',
                                     regs$chr1[i], regs$chr2[i])))

    ## Generate centered bins
    cbins <- seq(from = 0 - res*buffer, to = 0 + res*buffer, by = res)

    ## Expand centered bins (fast cross-join)
    expBins <- CJ(Var1=cbins, Var2=cbins)

    ## Replace x and y with the expanded, centered bins
    longMat[,`:=`(x = expBins$Var1,
                  y = expBins$Var2,
                  counts = counts), by = groupRow]

    ## Aggregate matrix by x and y position
    agg <- longMat[, .(counts = sum(counts)), by = .(x, y)]

    ## Cast to wide as a list of a single matrix
    wideMat <- list(reshape2::acast(agg, x ~ y, value.var = 'counts'))

    ## Append matrix list
    m <- append(m, wideMat)

  }

  ## Aggregate chr-chr matricies
  aggM <- Reduce('+', m)

  ## Close progress bar
  pb$tick(tokens = list(step = "Done!"))
  if(pb$finished) pb$terminate()

  ## Return aggregated matrix
  return(aggM)

}
