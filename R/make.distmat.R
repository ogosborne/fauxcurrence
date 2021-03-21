#' Make a distance matrix between all non-NA raster cells
#'
#' @description
#' This function produces a matrix of distances between all non-NA cells in a raster for use by fauxcurrence().
#'
#' @details
#' If dist_meth is set to "distm", the geosphere::distm function is used to calculate distances between all non-NA cells of rast using the function specified by dist_fun, for example "distGeo", "distHaversine" or "distVincentyEllipsoid" from the geosphere package. If dist_meth is set to "distRcpp", distRcpp::dist_mtom is used. This is significantly faster than geosphere::distm, but only two distance functions can be used: "Haversine" or "Vincenty". If dist_meth is set to "costdist", relative cost-distance is calculated using the gdistance package. Rast is also used for cost-distance calculation. Therefore, if relative overland distance is desired, values in the study region (i.e. land) should all be set to the same positive value, but if an arbitrary cost-distance calculation is desired, they should be set to the appropriate travel cost.
#'
#' @param rast A raster object.
#' @param dist_meth A string indicating the distance method, either "distm", "distRcpp" or "costdist".
#' @param dist_fun A string indicating the distance function if dist_meth is "distm" or "distRcpp". Must be either "Haversine" or "Vincenty" if dist_meth == "distRcpp", can be the name of any function which takes the same input and produces the same output as geosphere::distHaversine if dist_meth == "distm".
#' @return a matrix containing the distances between the centres of each non-NA cell of rast.
#' @export
#' @examples
#' library(raster)
#' my.raster <- raster(system.file("external/test.grd", package="raster"))
#' my.raster[!is.na(my.raster)] <- 1
#' my.distmat.overland <- make.distmat(rast=my.raster,dist_meth="costdist")
#' @importFrom geosphere distCosine distGeo distHaversine distm distVincentyEllipsoid distVincentySphere
#' @export
make.distmat <- function(rast,dist_meth="distRcpp",dist_fun="Haversine"){
  # find the xy coordinates of the centre of each cell
  rast.cent <- raster::xyFromCell(rast,raster::Which(!is.na(rast),cells=T))
  if(dist_meth == "distm"){
    # calculate a matrix of distances bwteen each cell using the function in
    distmat <- geosphere::distm(rast.cent,fun=get(dist_fun))
  } else if (dist_meth == "distRcpp"){
    if(!(dist_fun %in% c("Haversine", "Vincenty"))) stop('error: if dist_meth is "distRcpp", dist_fun must be either "Haversine" or "Vincenty"')
    distmat <- distRcpp::dist_mtom(rast.cent[,1],rast.cent[,2],rast.cent[,1],rast.cent[,2],dist_function = dist_fun)
  } else if(dist_meth == "costdist"){
    # make a geocorrected transition layer
    rast.tr <- gdistance::geoCorrection(gdistance::transition(rast,mean,16))
    # Use the transition layer to calculate a matrix of relative distances between all grid cell centres
    distmat <-gdistance::costDistance(rast.tr,rast.cent,rast.cent)
  } else {
    stop('error: dist_meth must be either "distm", "distRcpp" or "costdist"')
  }
  return(distmat)
}
