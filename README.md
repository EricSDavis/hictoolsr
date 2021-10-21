# hictoolsr: Hi-C analysis in R <img src="man/figures/logo.png" id="hictoolsr_logo" align="right" width="125"/>

`hictoolsr` is an R package with a collection of useful tools for Hi-C data analysis. The tools fall into the following broad categories:

-   Extracting Hi-C Counts

    Functions that are optimized to quickly extract and aggregate data from `.hic` files.

-   Generating Paired Interactions

    Methods to easily generate combinations of paired interactions.

-   Merging/Filtering/Binning Interactions

    A collection of functions for converting to `GInteractions` objects and manipulating anchors.

-   Plotting Functions

    Additional visualization functions that are compatible with the `plotgardener` package.

## Installation

This package can be installed via github:

```{r}
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("EricSDavis/hictoolsr")
```

## Workflows

We have provided a series of workflows that demonstrate how `hictoolsr` can be integrated and used with existing packages to analyze genomic data.

### Differential Loop Calling

This workflow demonstrates how to identify and visualize differential loops between two biological conditions using the `hictoolsr`, `DESeq2`, and `plotgardener` packages.

### Inter/Intrachromosomal Aggregate Analysis

(Coming soon) This workflow will demonstrate how to aggregate Hi-C signal from interactions between proteins bound on DNA.

## Future updates

Improve function naming:

-   `binBedpe()` args `a1Pos` and `a2Pos`
-   `binBedpe()` to `mapBedpeToBin()`, `mapBedpeToBin()`, `assignBedpeToBin()`, `convertBedpeToBin()`?
-   More general names (i.e. `binAnchor()` -\> `binGRanges()`)?

`binBedpe()`

-   Improve documentation and examples

`calcApa()`

-   Improve documentation and examples

`calcBepde()`

-   Improve documentation and examples
-   Change output from `data.table` to `GInteractions` object.
-   Replace internal binning method with other binning functions

`extractCounts()`

-   Improve documentation and examples

-   Change name to `extractHicCounts()`

-   Remove excess comments in function

-   Move binning warning to its own function

-   Utilize the chromosome extraction method from `calcApa()`

    -   This might remove the chromosomes argument as well.

`filterBedpe()`

-   Improve documentation and examples

-   Currently does not filter out regions that are too close to the starts and ends of chromosomes. This is not a problem (it doesn't affect how straw imports the data), but it might be desirable to have those interactions filtered out.

    -   Specifically, remove `starts - res*buffer` and `end + res*buffer` in filtering.
    -   Will require `GRanges` to have `seqinfo`

-   Move binning warning to its own function

-   Currently removes mcols after filtering - make sure to add those back in.

`mergeBedpe()`

-   Improve documentation and examples.
-   Change input to handle R objects (data.tables/GInteraction objects).
-   Replace internal binning method with other binning functions.
-   Change output from `data.table` to `GInteractions` object (or make this an argument).

`plotApa()`
