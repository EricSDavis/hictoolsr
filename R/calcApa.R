#' Extract BEDPE interactions from a Hi-C matrix
#'
#' @param bedpe data.frame or data.table in BEDPE format.
#' @param hicFile string path to .hic file.
#' @param norm string hic normalization <NONE/VC/VC_SQRT/KR>.
#' @param res integer resolution of bedpe bins.
#' @param buffer integer number of res-length bins from the center pixel.
#' @param equalLengthOut logical can be set to false to get an APA list that is not the
#'   same length as the number of rows of the input bedpe file. Adds a slight efficieny
#'   boost.
#'
#' @return Returns a list APA matricies that can be compiled and plotted.
#'
#' @import data.table strawr
#' @importFrom tidyr complete full_seq
#'
#' @export
#'
calcApa <- function(bedpe, hicFile, norm = "NONE", res = 10000, buffer = 5, equalLengthOut = T) {

  ## Put data in correct order for straw extraction
  if (all(bedpe[[1]] == bedpe[[4]])) {

    ## Intrachromosomal: Flip columns where x1 > y1 if chr are equal
    ## Convert chromosome factors to characters
    bedpe[[1]] <- as.character(bedpe[[1]])
    bedpe[[4]] <- as.character(bedpe[[4]])

    ## Flip columns where x1 > y1 if chr are equal
    bedpe[bedpe$y1 - bedpe$x1 < 0,] <- bedpe[bedpe$y1 - bedpe$x1 < 0, c(4:6, 1:3)]

  } else {

    ## Interchromosomal: Flip columns where chr1 > chr2
    ## Replace chromosome X, Y with 23, 24
    bedpe[[1]][bedpe[[1]] == "X"] <- 23
    bedpe[[4]][bedpe[[4]] == "X"] <- 23
    bedpe[[1]][bedpe[[1]] == "Y"] <- 23
    bedpe[[4]][bedpe[[4]] == "Y"] <- 23

    ## Coerce chromosome strings to numeric
    bedpe[[1]] <- as.numeric(bedpe[[1]])
    bedpe[[4]] <- as.numeric(bedpe[[4]])

    ## Flip anchor order where chr1 > chr2
    bedpe[bedpe[[1]] > bedpe[[4]],] <- bedpe[bedpe[[1]] > bedpe[[4]], c(4:6, 1:3)]

    ## Coerce chromosome numerics to strings
    bedpe[[1]] <- as.numeric(bedpe[[1]])
    bedpe[[4]] <- as.numeric(bedpe[[4]])

    ## Replace 23, 24 with chromosome X, Y
    bedpe[[1]][bedpe[[1]] == "23"] <- "X"
    bedpe[[4]][bedpe[[4]] == "23"] <- "X"
    bedpe[[1]][bedpe[[1]] == "24"] <- "Y"
    bedpe[[4]][bedpe[[4]] == "24"] <- "Y"

  }

  ## Expand expand pixel to range
  bedpe2    <- bedpe
  bedpe2$x1 <- bedpe$x1 - (res*buffer)
  bedpe2$x2 <- bedpe$x1 + (res*buffer)
  bedpe2$y1 <- bedpe$y1 - (res*buffer)
  bedpe2$y2 <- bedpe$y1 + (res*buffer)

  ## Find all counts within slice of data
  counts <- lapply(1:nrow(bedpe2), function(i) {
    as.data.table(straw(norm = norm, fname = hicFile,
                        chr1loc = paste0(bedpe2[i,1:3], collapse = ":"),
                        chr2loc = paste0(bedpe2[i,4:6], collapse = ":"),
                        unit = "BP", binsize = res, matrix = "observed"))
  })

  if (!equalLengthOut) {
    ## Filter out empty counts
    bedpe2 <- bedpe2[lapply(counts, nrow) > 0]
    counts <- counts[lapply(counts, nrow) > 0]
  }

  ## Create vectorized version of seq function
  seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"), SIMPLIFY = F)

  ## Get rows and columns for matrix
  rows <- seq2(from = bedpe2$y1, to = bedpe2$y2, by = res)
  cols <- seq2(from = bedpe2$x1, to = bedpe2$x2, by = res)

  ## Fill missing combinations
  ct <- lapply(1:length(counts), function(i) {
    as.data.table(complete(counts[[i]],
                           x = full_seq(cols[[i]], res),
                           y = full_seq(rows[[i]], res),
                           fill = list(counts = 0)))
  })

  ## Cast to wide format
  dc <- lapply(ct, function(i) {
    dcast.data.table(i, y ~ x, value.var = "counts", fill = 0, drop = F)
  })

  ## Convert to matrix
  dm <- lapply(dc, function(i) {
    as.matrix(i[,-1], rownames = i$y)
  })

  ## Return list of apa matricies
  return(dm)
}
