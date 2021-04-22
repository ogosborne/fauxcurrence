# fauxcurrence

This R package generates sets of randomised occurrences (long & lat coordinates) for one or more species with a spatial structure (within- and between-species distances) matching that of an observed set of coordinates. These can be used to define null expectations to test hypotheses in ecology and biogeography.

To install, use:

```R
devtools::install_github("ogosborne/fauxcurrence")
```

To learn more, see the man pages and run the example in the vignette:

```R
library(fauxcurrence)
?fauxcurrence
?make.distmat
vignette("Using-Fauxcurrence")
```
