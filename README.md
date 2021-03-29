# hictoolsr: Hi-C analysis in R

hictoolsr is an R package with a collection of useful tools for Hi-C data analysis. The tools in this package are optimized to quickly extract and apply data from `.hic` files without the use of external tools.


## Future updates

Improve function naming:

* `binBedpe()` args `a1Pos` and `a2Pos`
* `binBedpe()` to `mapBedpeToBin()`, `mapBedpeToBin()`, `assignBedpeToBin()`, `convertBedpeToBin()`?
* More general names (i.e. `binAnchor()` -> `binGRanges()`)?


`binBedpe()`

* Improve documentation and examples


`calcApa()`

* Improve documentation and examples


`calcBepde()`

* Improve documentation and examples
* Change output from `data.table` to `GInteractions` object.
* Replace internal binning method with other binning functions


`extractCounts()`

* Improve documentation and examples
* Change name to `extractHicCounts()`
* Remove excess comments in function
* Move binning warning to its own function
* Utilize the chromosome extraction method from `calcApa()`
  * This might remove the chromosomes argument as well.


`filterBedpe()`

* Improve documentation and examples
* Currently does not filter out regions that are too close to the starts and ends of chromosomes. This is not a problem (it doesn't affect how straw imports the data), but it might be desirable to have those interactions filtered out.
  * Specifically, remove `starts - res*buffer` and `end + res*buffer` in filtering.
  * Will require `GRanges` to have `seqinfo`
* Move binning warning to its own function


`mergeBedpe()`

* Improve documentation and examples.
* Change input to handle R objects (data.tables/GInteraction objects).
* Replace internal binning method with other binning functions.
* Change output from `data.table` to `GInteractions` object (or make this an argument).


`plotApa()`
