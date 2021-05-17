#' Extract counts from .hic file(s)
#'
#' @param bedpe GInteraction object that has been binned to the correct resolution
#' @param hic Path to .hic file
#' @param chroms Character vector of chromosomes to extract
#' @param res Resolution of bedpe bins
#' @param norm string hic normalization <NONE/VC/VC_SQRT/KR>.
#'
#' @export
#'
extractCounts <- function(bedpe, hic, chroms = c(1:22, 'X', 'Y'),
                          res = 10000, norm = 'NONE', matrix = 'observed') {

  ## TODO
  ## Write helper function to bin bedpe interactions by resolution
  ##    Make it flexible enough to handle data.tables/frames and GInteraction objects
  ## Apply helper functino in extract Counts

  ## Start progress bar
  pb <- progress::progress_bar$new(
    format = "  :step [:bar] :percent elapsed: :elapsedfull",
    clear = F, total = length(hic)*length(chroms)*3+1)
  pb$tick(0)

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

  ## Convert chroms to character
  chroms <- as.character(chroms)

  ## Split binned binnedEP by chromosome
  chrBedpe <- lapply(chroms, function(chr) {
    bedpe[seqnames(first(bedpe)) == paste0('chr',chr)]
  })

  ## Add names
  names(chrBedpe) <- chroms

  ## Loop through each hic file and chromosome
  for(i in 1:length(hic)) {
    for(j in 1:length(chroms)) {

      ## Update progress
      pb$tick(tokens=list(step=sprintf('Extracting chr%s from %s',
                                       chroms[j], basename(hic[i]))))

      ## Construct straw format from binnedEP
      ep <- data.table(
        x = start(first(chrBedpe[[chroms[j]]])),
        y = start(second(chrBedpe[[chroms[j]]])),
        counts = 0
      )

      ## Swap row assignment for out of order x and y
      ep[y < x, `:=`(x=y,y=x)]

      ## Pull out the entire chromosome sparse matrix
      sparseMat <- as.data.table(strawr::straw(norm = norm, fname = hic[i],
                                               chr1loc = chroms[j],
                                               chr2loc = chroms[j],
                                               unit = "BP",
                                               binsize = res,
                                               matrix = matrix))

      ## Update progress
      pb$tick(tokens=list(step=sprintf('Extracting chr%s from %s',
                                       chroms[j], basename(hic[i]))))

      ## Set keys
      setkeyv(sparseMat, c('x', 'y'))

      ## Get counts by key
      ep$counts <- sparseMat[ep]$counts

      ## Set unmatched counts (NA) to 0
      ep[is.na(counts), counts := 0]

      ## Add counts back to chrBedpe
      mcols(chrBedpe[[chroms[j]]])[basename(hic[i])] <- ep$counts

      ## Update progress
      pb$tick(tokens=list(step=sprintf('Extracting chr%s from %s',
                                       chroms[j], basename(hic[i]))))

    }
  }

  ## Combine results
  chrBedpe <- unname(chrBedpe)
  combinedResults <- do.call(c, chrBedpe)

  ## Close progress bar
  pb$tick(tokens = list(step = "Done!"))
  if(pb$finished) pb$terminate()

  ## Return combined results
  return(combinedResults)
}


# ## Extract loop counts for each biorep hic file
# result <- extractCounts(bedpe = binnedEP,
#                         hic = c("../hic/CI_THP1_O_1_inter_30.hic",
#                                 "../hic/CI_THP1_O_2_inter_30.hic",
#                                 "../hic/CI_THP1_A_1_inter_30.hic",
#                                 "../hic/CI_THP1_A_2_inter_30.hic"),
#                         chroms = c(1:22,'X','Y'),
#                         res = 10000)
#
# result
