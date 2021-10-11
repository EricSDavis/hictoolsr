## Script to generate loop counts for differential
## loop calling vignette

## Load packages
library(hictoolsr)
library(dbscan)

## Define WT and FS loop file paths
wt_loops <- system.file("extdata/WT_5kbLoops.txt", package = "hictoolsr")
fs_loops <- system.file("extdata/FS_5kbLoops.txt", package = "hictoolsr")

## Merge loops and convert to GInteractions
loops <-
  mergeBedpe(bedpeFiles = c(wt_loops, fs_loops), res = 10e3) |>
  as_ginteractions()

## Download biorep hic data from GEO
## (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143465)
hicFiles <-
  list.files(path = "../../largeData/hic/LEUK/replicate",
             full.names = TRUE) |>
  normalizePath()

# ## Extract Hi-C counts between loop pixels
# loopCounts <- extractCounts(bedpe = loops,
#                             hic = hicFiles,
#                             chroms = c(1:22, "X"),
#                             res = 10e3,
#                             norm = 'NONE',
#                             matrix = 'observed')

## Extract Hi-C counts between loop pixels (in parallel)
system.time({
  loopCounts <-
    BiocParallel::bplapply(X = c(1:22, 'X'), FUN = \(x) {
      message("chr",x)
      extractCounts(bedpe = loops,
                    hic = hicFiles,
                    chroms = x,
                    res = 10e3,
                    norm = 'NONE',
                    matrix = 'observed')
    })

  loopCounts <- do.call(c, loopCounts)
})

## Save loop counts
save(loopCounts, file = "data/loopCounts.rda")

