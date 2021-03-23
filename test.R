library(hictoolsr)

## Use annotation hub to get chipseq data
library(AnnotationHub)
ah <- AnnotationHub()

## Search for data sets
query(ah, c("Homo", "Gm12878", "ctcf", "narrowPeak"))

## Get some narrowPeak data
ctcf <- ah[["AH27499"]]


## Calculate interactions
intra <- calcBedpe(bed = ctcf, res = 10000, x = 100,
                   interactions = "intrachromosomal")

inter <- calcBedpe(bed = ctcf, res = 10000, x = 100,
                   interactions = "interchromosomal")


## Filter intrachromosomal interactions
intra <- filterBedpe(bedpe = intra, res = 10000, buffer = 5)


## Calculate APA
apa <- calcApa(bedpe = intra,
               hicFile = "../../projects/CHON/data/hic/CHON_HiC_C28_WT_NA_0_S_inter_30.hic",
               norm = 'NONE', res = 10000, buffer = 5)



library(BentoBox)

plotApa(apa = apa, palette = colorRampPalette(brewer.pal(n=9,"YlGnBu")))
plotApa(apa = apa)

bb_pageCreate(width = 4, height = 4)

p <- bb_params(x = 0.5, y = 0.5, width = 3, height = 3,
               palette = colorRampPalette(brewer.pal(n=9,"YlGnBu")))

plotApa(apa = apa, params = p)
