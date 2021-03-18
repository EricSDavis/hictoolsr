#' Read-in and merge BEDPE Interactions by Genomic distance with DBSCAN
#'
#' @description
#' Given a list of BEDPE files, a resolution, and a selector column, this function will
#' use dbscan to merge loops by genomic distance.
#'
#' @param bedpeFiles Character vector of BEDPE file paths to be merged.
#' @param res Integer - resolution in base-pairs of all BEDPE interactions. Used to
#'  determine the epsilon value for \code{dbscan()} (i.e. \code{eps = res*2}).
#' @param selectCol Integer or character string denoting the column to be used to select
#'  clustered interactions.
#' @param dist_method Character string - distance measure used in the dist function. For
#'  available methods see \code{?dist()}.
#' @param minPts Integer - minimum number of interactions to form a DBSCAN cluster.
#'
#' @return Returns a data table of selected BEDPE interactions
#'
#'
#' @export
#'
mergeBedpe <- function(bedpeFiles,
                       res=5000,
                       selectCol=12,
                       dist_method="manhattan",
                       minPts=2) {

  ## Read in list of bedpe
  bedpe <- lapply(bedpeFiles, fread)

  ## Add new column for source
  bedpe <- Map(cbind, bedpe, source = basename(bedpeFiles))

  ## Concatentate bedpe
  bedpe <- do.call(rbind, bedpe)

  ## Confirm that all comparisons are intrachromosomal
  if(!all(bedpe[,1] == bedpe[,4])) stop("All interactions must be intrachromosomal")

  ## Bin BEDPE interactions to the specified resolution
  bedpe[,c(2,5)] <- floor(bedpe[,c(2,5)]/res)*res
  bedpe[,c(3,6)] <- bedpe[,c(2,5)] + res

  ## Split by chromosome
  chrBedpe <- split(x = bedpe, f = bedpe[,1])

  ## Calculate distance between x1 and y1 coordinates for each chromosome
  d <- lapply(chrBedpe, function(x) dist(x[,c(2,5)], method = dist_method))

  ## Cluster with dbscan for each chromosome
  clst <- lapply(d, function(x) dbscan(x, eps = res*2, minPts = minPts)$cluster)

  ## Add cluster info to chromosome bedpe list
  chrBedpe <- Map(cbind, chrBedpe, cluster = clst)

  ## Split by groups within each chromosome
  chrBedpeClust <- lapply(chrBedpe, function(x) split(x, x$cluster))

  ## Separate single interactions from clustered interactions per chromosome
  singleBedpe <- lapply(chrBedpeClust, `[[`, 1)
  clustBedpe  <- lapply(chrBedpeClust, `[`, -1)

  ## Select the "strongest interaction" (12) for each cluster (y) and chromosome (x)
  selectedBedpe <- lapply(clustBedpe, function(x) {
    do.call(rbind, lapply(x, function(y) y[which.max(y[[selectCol]]),]))
  })

  ## Join single and selected interactions by chromosome
  mergedChrBedpe <- Map(rbind, singleBedpe, selectedBedpe)

  ## Sort result by start position
  mergedChrBedpe <- lapply(mergedChrBedpe, function(x) x[order(x[,2]),])

  ## Combine into a single data table
  mergedBedpe <- do.call(rbind, mergedChrBedpe)

  ## Return final BEDPE data table
  return(mergedBedpe)

}
