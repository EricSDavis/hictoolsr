#' Calculates all potential paired interactions from a BED file
#'
#' @param bed data.frame in narrowpeak format
#' @param res bin resolution for output BEDPE
#' @param x integer top x numer of strongest or weakest peaks
#' @param decreasing logical order by decreasing
#' @param interactions string filter for interchromosomal or intrachromosomal interactions
#'
#' @return Returns a data.table of BEDPE interactions
#'
#' @export
#'
calcBedpe <- function(bed, res = 10000, x = 100, decreasing = TRUE,
                         interactions = "interchromosomal") {

  ## Format bed if it is a granges object
  if ("GRanges" %in% class(bed)) {
    bed <- as.data.table(bed)
    bed[, width:=NULL]
  }

  ## Ensure bed is in narrowPeak format (i.e. 10 columns)
  stopifnot(ncol(bed) == 10)

  ## Define narrowPeak column names
  colnames(bed) <- c("chrom", "chromStart", "chromEnd",
                     "name", "score", "strand", "signalValue",
                     "pValue", "qValue", "peak")

  ## Take the top x strongest or weakest peaks
  bed <- bed[order(bed$signalValue, decreasing = decreasing)][1:x,]

  ## Bin the bed file based on the resolution
  binned <- data.table(
    chrom = bed$chrom,
    binStart = floor((bed$chromStart + bed$peak)/res)*res,
    binEnd = floor((bed$chromStart + bed$peak)/res)*res + res
  )

  ## Remove duplicate bins (multiple peaks that fall into the same hic bin)
  binned <- unique(binned)

  ## Define number of rows
  n <- nrow(binned)

  ## Define empty matrix of the dimensions of the data
  m <- matrix(nrow = n, ncol = n)

  ## Fill in matrix with all i (rows) and j (column) comparisons
  for(i in 1:n){
    for(j in 1:n){
      m[i,j] <- paste0(c(i,j), collapse = "_")
    }
  }

  ## Take the upper triangular, split rows and combine into coordinates
  coords <- do.call(rbind, strsplit(m[upper.tri(m)], "_"))

  ## Coerce cooridnates to numeric
  coords <- apply(coords, 2, as.numeric)

  ## Extract coordinates from binned file
  bedpe <- cbind(binned[coords[,1],], binned[coords[,2],])
  colnames(bedpe) <- c("chr1", "x1", "x2", "chr2", "y1", "y2")

  ## Filter for inter or intrachromosomal interactions
  if (interactions == "interchromosomal")
  {
    ## Remove intra-chromosomal pairs
    bedpe <- bedpe[chr1 != chr2]
  }
  else if (interactions == "intrachromosomal")
  {
    ## Remove inter-chromosomal pairs
    bedpe <- bedpe[chr1 == chr2]
  }
  else {stop("interactions must be either 'interchromosomal' or 'intrachromosomal'.")}

  ## Remove chr prefix & non-standard chroms
  bedpe$chr1 <- gsub("chr", "", bedpe$chr1)
  bedpe$chr2 <- gsub("chr", "", bedpe$chr2)
  bedpe <- bedpe[bedpe$chr1 %in% c(1:22, "X", "Y") & bedpe$chr2 %in% c(1:22, "X", "Y"),]

  ## Return bedpe
  return(bedpe)

}
