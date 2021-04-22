# fauxcurrence

This R package generates sets of randomised occurrences (long & lat coordinates) for one or more species with a spatial structure (within- and between-species distances) matching that of an observed set of coordinates. These can be used to define null expectations to test hypotheses in ecology and biogeography.

To install, use:

```R
devtools::install_github("ogosborne/fauxcurrence",build_vignettes = TRUE)
```

To learn more, see the man pages and run the example in the vignette:

```R
# load fauxcurrence
library(fauxcurrence)

# view manual files
?fauxcurrence
?make.distmat

# view vignette
vignette("Using-Fauxcurrence")
```
