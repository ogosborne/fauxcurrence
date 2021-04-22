## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(fauxcurrence)
library(raster)

## -----------------------------------------------------------------------------
data(chameleons)
data(madagascar)

## -----------------------------------------------------------------------------
head(chameleons)

## ---- fig.dim = c(6, 7), out.extra = 'style="border:0"'-----------------------
# set up a colour palette to represent the species
my.cols <- rainbow(length(unique(chameleons$species)))
names(my.cols) <- sort(unique(chameleons$species))
# plot the raster of madagascar
plot(madagascar,col=c("white","gray85"),box=F,axes=F,legend=F)
# plot the points and a legend
points(chameleons,bg=my.cols[chameleons$species],pch=21)
legend("bottomright",legend=names(my.cols), pt.bg=my.cols,pch=21)

## ---- out.width = "50%", echo = FALSE, out.extra = 'style="border:0"'---------
knitr::include_graphics("flowchart.png")

## -----------------------------------------------------------------------------
set.seed(1234)

## -----------------------------------------------------------------------------
faux_intra <- fauxcurrence(coords = chameleons, rast = madagascar)

## -----------------------------------------------------------------------------
str(faux_intra)

## -----------------------------------------------------------------------------
head(faux_intra$points)

## ---- fig.dim = c(6, 8), out.extra = 'style="border:0"'-----------------------
# set the plotting device to use two panels so we can display the observed and simulated data together
par(mfrow=c(1,2))
# plot the raster of madagascar
plot(madagascar,col=c("white","gray85"),box=F,axes=F,legend=F,main="observed")
# plot the original points 
points(chameleons,bg=my.cols[chameleons$species],pch=21)
# plot the raster of madagascar
plot(madagascar,col=c("white","gray85"),box=F,axes=F,legend=F,main="intra")
# plot the simulated points and a legend
points(faux_intra$points,bg=my.cols[faux_intra$points$species],pch=21)
legend("bottomright",legend=names(my.cols), pt.bg=my.cols,pch=21)

## ---- fig.dim = c(6, 7), out.extra = 'style="border:0"'-----------------------
# set seed
set.seed(4321)
# run fauxcurrence 
faux_inter <- fauxcurrence(coords = chameleons, rast = madagascar, inter.spp = TRUE, verbose=FALSE)
# set the plotting device to use two panels so we can display the observed and simulated data together
par(mfrow=c(1,2))
# plot the raster of madagascar
plot(madagascar,col=c("white","gray85"),box=F,axes=F,legend=F,main="observed")
# plot the original points 
points(chameleons,bg=my.cols[chameleons$species],pch=21)
# plot the raster of madagascar
plot(madagascar,col=c("white","gray85"),box=F,axes=F,legend=F,main="inter")
# plot the simulated points and a legend
points(faux_inter$points,bg=my.cols[faux_inter$points$species],pch=21)
legend("bottomright",legend=names(my.cols), pt.bg=my.cols,pch=21)

## ---- fig.dim = c(6, 7), out.extra = 'style="border:0"'-----------------------
# set seed
set.seed(54321)
# run fauxcurrence 
faux_interSep <- fauxcurrence(coords = chameleons, rast = madagascar, inter.spp = TRUE, sep.inter.spp = TRUE, verbose=FALSE)
# set the plotting device to use two panels so we can display the observed and simulated data together
par(mfrow=c(1,2))
# plot the raster of madagascar
plot(madagascar,col=c("white","gray85"),box=F,axes=F,legend=F,main="observed")
# plot the original points 
points(chameleons,bg=my.cols[chameleons$species],pch=21)
# plot the raster of madagascar
plot(madagascar,col=c("white","gray85"),box=F,axes=F,legend=F,main="inter-sep")
# plot the simulated points and a legend
points(faux_interSep$points,bg=my.cols[faux_interSep$points$species],pch=21)
legend("bottomright",legend=names(my.cols), pt.bg=my.cols,pch=21)

