## Load package
library(InteractionSet)
library(DESeq2)
library(hictoolsr)

## Call differential loops ---------------------------------------------------------------

## Load loop count data
data("loopCounts")

## Simplify column names
colnames(mcols(loopCounts)) <-
  gsub(pattern = "GSM.*_IDR_(WT|FS)_A9_(1|2)_(1|2)_.*",
       replacement = "\\1_\\2_\\3",
       x = colnames(mcols(loopCounts)))

## Isolate count matrix
cnts <-
  mcols(loopCounts)[grep("WT|FS", colnames(mcols(loopCounts)))] |>
  as.matrix()

## Create colData from column names
colData <-
  do.call(rbind, strsplit(x = colnames(cnts), split = "_")) |>
  as.data.frame(stringsAsFactors = TRUE) |>
  `colnames<-`(value = c("condition", "biorep", "techrep"))

## Build DESeq data set
dds <-
  DESeqDataSetFromMatrix(countData = cnts,
                         colData = colData,
                         design = ~techrep + biorep + condition)

## Run DEseq analysis
res <-
  DESeq(dds) |>
  lfcShrink(coef = "condition_WT_vs_FS", type="apeglm")

## Attach DESeq2 results
mcols(loopCounts) <- cbind(mcols(loopCounts), res)

## Separate WT/FS-specific loops
wtLoops <- loopCounts[loopCounts$padj <= 0.01 &
                        loopCounts$log2FoldChange > 0]

fsLoops <- loopCounts[loopCounts$padj <= 0.01 &
                        loopCounts$log2FoldChange < 0]

## Make loop APAs ------------------------------------------------------------------------

## List loop categories
loopList <-
  list(allLoops = loopCounts,
       wtLoops = wtLoops,
       fsLoops = fsLoops)

## Define resolution and buffer (pixels from center)
res <- 10e3
buffer <- 10

## Filter out short interactionss
filteredLoops <-
  lapply(X = loopList,
         FUN = filterBedpe,
         res = res,
         buffer = buffer) |>
  `names<-`(value = names(loopList))

## Download Hi-C data from GEO
## (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143465)
wtHicPath <- "../../largeData/hic/LEUK/condition/HEK_HiC_NUP_IDR_WT_A9_megaMap_inter_30.hic"
fsHicPath <- "../../largeData/hic/LEUK/condition/HEK_HiC_NUP_IDR_FS_A9_megaMap_inter_30.hic"

## Extract interactions from WT and FS hic maps
loopApaWtHic <-
  lapply(X = filteredLoops,
         FUN = calcApa,
         hic = wtHicPath,
         norm = "KR",
         buffer = buffer)

loopApaFsHic <-
  lapply(X = filteredLoops,
         FUN = calcApa,
         hic = fsHicPath,
         norm = "KR",
         buffer = buffer)

## Save APA results
save(loopApaWtHic,  file = "data/loopApaWtHic.rda")
save(loopApaFsHic,  file = "data/loopApaFsHic.rda")
