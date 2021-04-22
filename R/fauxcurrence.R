

KL.div<-function(vec1,vec2,break.num,breaks=NULL) {
  #This function calculates the Kullback-Leibler divergence between the distribution of observed interpoint distances and the null interpoint distances. Breaks are defined based on the range of vec1
  if (is.null(breaks)) {
    rng<-range(vec1)
    breaks<-seq(rng[1],rng[2],diff(rng)/break.num)
  }
  counts1<-hist(vec1,breaks=breaks,plot=F,)$count
  counts2<-hist(vec2,breaks=breaks,plot=F,)$count
  d1<-(counts1+1)/(length(vec1)+length(counts1))
  d2<-(counts2+1)/(length(vec2)+length(counts2))
  KL<-sum(d1*log10(d1/d2))
  return(list(KL,breaks))
} # end


get_dists <- function(coords,species.vec,dist_meth="distRcpp",dist_fun="Haversine",inter.spp=TRUE,sep.inter.spp=TRUE,combs=NULL){
  # initiate output
  dist.list <- list()
  # get distances
  if(dist_meth == "distm"){
    dists <- geosphere::distm(coords,fun=get(dist_fun))
  } else if (dist_meth == "distRcpp"){
    if(!(dist_fun %in% c("Haversine", "Vincenty"))) stop('if dist_meth is "distRcpp", dist_fun must be either "Haversine" or "Vincenty"')
    dists <- distRcpp::dist_mtom(coords[,"x"],coords[,"y"],coords[,"x"],coords[,"y"],dist_function = dist_fun)
  } else {
    stop('if use.distmat is FALSE, dist_meth must be either "distm" or "distRcpp"')
  }
  diag(dists)<-NA
  # extract intraspecific distances
  for (sp in sort(unique(species.vec))){
    in.coords <- which(species.vec == sp)
    dist.list[[paste("intra",sp,sep="_")]] <- dists[in.coords,in.coords][lower.tri(dists[in.coords,in.coords])]
  }
  # extract interspecific distances
  if(inter.spp){
    if(sep.inter.spp){
      for (comb in colnames(combs)){
        #coords.1 <- which(species.vec == combs[1,comb])
        #coords.2 <- which(species.vec == combs[2,comb])
        dist.list[[comb]] <- as.vector(dists[which(species.vec == combs[1,comb]),which(species.vec == combs[2,comb])])
      }
    } else {
      for (sp in sort(unique(species.vec))){
        in.coords <-which(species.vec == sp)
        out.coords <- which(species.vec != sp)
        dist.list[[paste("inter",sp,"all",sep="_")]] <- as.vector(dists[in.coords,out.coords])
      }
    }
  }
  return(dist.list)
} # end

get_dists_from_mat <- function(mat,indices,species.vec,inter.spp=TRUE,sep.inter.spp=TRUE,combs=NULL){
  # arguments
  # mat: matrix of distances for between all non-NA points in a raster
  # indices: indices of mat which correspond to all occurrence points
  # species.vec: vector of species names for all occurrence points
  #
  dist.list <- list()
  ## get intraspecific distances
  for (sp in sort(unique(species.vec))){
    my.in.idx <- indices[which(species.vec == sp)]
    my.in.submat <- mat[my.in.idx,my.in.idx]
    dist.list[[paste("intra",sp,sep="_")]] <- my.in.submat[lower.tri(my.in.submat,diag = FALSE)]
  }
  ## get interspecific distances
  if(inter.spp){
    if(sep.inter.spp){
      for (comb in colnames(combs)){
        my.1.idx <- indices[which(species.vec == combs[1,comb])]
        my.2.idx <- indices[which(species.vec == combs[2,comb])]
        my.out.submat <- mat[my.1.idx,my.2.idx]
        dist.list[[comb]] <- as.vector(my.out.submat)
      }
    } else {
      for (sp in unique(species.vec)){
        my.in.idx <- indices[which(species.vec == sp)]
        my.out.idx <- indices[which(species.vec != sp)]
        my.out.submat <- mat[my.in.idx,my.out.idx]
        dist.list[[paste("inter",sp,"all",sep="_")]] <- as.vector(my.out.submat)
      }
    }
  }
  return(dist.list)
} # end


#' Simulate species occurrence data
#'
#' @description
#' This function generates a set of randomised occurrences (long & lat coordinates) for one or more species with a spatial structure (within- and between-species distances) matching that of an observed set of coordinates.
#'
#' @details

#' Within-species distances are always used. If inter.spp is TRUE, it also uses interspecies distances, and if sep.inter.spp is TRUE it separates interspecific distances into all pairwise sets of species. Thus, occurrences for each species can be simulated independently (inter.spp=FALSE, sep.inter.spp=FALSE), general distance to heterospecifics can be taken into account (inter.spp=TRUE, sep.inter.spp=FALSE), or distance relationships between individual pairs of species can be preserved (inter.spp=TRUE, sep.inter.spp=TRUE).
#'
#' The algorithm starts by generating an initial set of randomised occurrences which fit within the ranges of observed interpoint distances. There are two stages to this process. In stage 1, it places one point for each species on the map provided with rast and, if inter.spp is TRUE, checks that the distances between species are within the observed range (or within the central init.range quantile of the real distribution if trim.init.range is TRUE). If interspecific distances are outside the real distribution (or init.range), it iteratively improves the initial points by randomly replacing a point and rechecking the distances until they are within the desired range. If fix.seed.pts is provided, stage 1 is skipped and these points are used instead. The points in fix.seed.points are not replaced in later stages of the algorithm, so they can be used to define the rough centroids of the species distributions. In stage 2, more points are added (using the method set by new.pt.meth.stg.2) until the correct number of points for each species is reached. As with stage 1, after each point is added, distances are checked to ensure they are within the observed range of distances. To stop stages 1 and 2 becoming stuck with a previous set of points which make step-wise iterative improvement impossible, iter.max.stg1 and iter.max.stg2 set an upper limit on the number of iterations without improvement. If these are reached, initial point generation is restarted.
#'
#' The algorithm then implements an iterative improvement procedure to improve the similarity of the spatial structure between the observed and simulated points. For each iteration, one (or several if switch.n > 1) change(s) is made to the points. This consists of randomly replacing a point (using the method set by new.pt.meth.stg.3) and, if the change improves the match of the spatial structure between the null and observed points, it is kept. If the match is not improved, it is discarded and the original is retained. This procedure is repeated a for a user-defined number of iterations until either no improvement has been made for div.n.flat iterations or the total number of iterations has reached iter.max.stg3. The match between spatial structures is evaluated using the probability distribution of the interpoint distances. The comparison is made using the discrete version of the Kullback-Leibler divergence, where smaller values indicate a better match between the simulated and observed points. Since there are multiple interpoint distance distributions (i.e. intra- and inter-specific distances for multiple species or pairs of species), a weighted mean of divergence statistics across all distance distributions is used, weighted such that interpoint and intrapoint distances contribute equally. The divergence statistic should be plotted against number of iterations to ensure the statistic has been minimised (i.e. the curve has gone flat), a basic ASCII plot showing N iterations vs divergence is printed to the console (or to a log file if logfile!=NULL) if verbose == TRUE.
#'
#'
#' @param coords A dataframe of observed occurrences. There should be columns for longitude and latitude, named "x" and "y" respectively, and an optional third column named "species" with species identities - required only if there is more than one species.
#' @param rast A raster object defining the extent of the study area. Also used for distance matrix calculation if use.distmat is TRUE but distmat is not provided (see \code{?make.distmat} for more details).
#' @param distmat A matrix of distances between all non-NA cells of the raster in order, as can be created with make.distmat(). Used if use.distmat == TRUE.
#' @param use.distmat Logical indicating whether a precomputed distance matrix (produced with make.distmat()) should be used for inter-point distances. If TRUE, and distmat is not provided, one is created with make.distmat(). This is usually faster than use.distmat==FALSE, but the RAM requirements for larger rasters can make it computationally inviable.
#' @param inter.spp Logical indicating whether interspecific distances should be used in addition to intraspecific distances
#' @param sep.inter.spp Logical indicating whether interspecific distances should be seperated into pairwise species combinations. If FALSE, there is one interspecific distance measure for each species which is the distance from that species' points to all heterospecific points, whereas if TRUE there is a separate distance measure for each pairwise combination of species. Ignored if inter.spp is FALSE.
#' @param fix.seed.pts A dataframe of points that will be used as the first point for each species and will not be moved during stage 3. This should be in the same format as coords with one point for each species.
#' @param div.n.flat Integer specifying the number of iterations for which divergence values should remain unchanged before stopping iterations in the the stage 3 iterative improvement procedure. Should be high enough for divergence to minimise.
#' @param allow.ident.conspec Logical indicating whether conspecific points in the same place should be allowed. This should match the input data i.e. if the input is "presence only" data, it should be set to FALSE but if it is raw occurrence data with multiple conspecific points in the same locations, it should be set to TRUE.
#' @param dist_meth A string indicating the distance method to use. Should be either "distRcpp", "distm" or "costdist" ("costdist" is only implemented when use.distmat is TRUE).
#' @param dist_fun A string indicating the distance function to be used for distance calculations by geosphere::distm (if dist_meth is "distm") or distRcpp::dist_mtom (if dist_methis "distRcpp"). For dist_meth=="distRcpp", this should either be "Haversine" or "Vincenty". For dist_meth=="distm", it can be the name of any loaded function which takes the same input and produces the same output as geosphere::distHaversine.
#' @param trim.init.range Logical indicating whether the initial points generated in stage 1 should be within the central init.range quantile of the observed distribution rather than just the total range. This can speed up stage 2, especially if many species are present.
#' @param init.range Numeric between 0 and 1 defining the central range to constrain initial points to if trim.init.range is TRUE.
#' @param new.pt.meth.stg.2 String indicating the method of sampling a new point in stage 2. Either "sample", which samples a random point from the raster; "dist.obs" which randomly chooses an existing point and places a new point of the same species \emph{D} distance away, where \emph{D} is sampled from the observed intraspecific distances for that species; or "dist.dens" (recommended) which is similar to dist.obs but \emph{D} is sampled from a density object constructed from the observed intraspecific distances for that species.
#' @param new.pt.meth.stg.3 As with new.pt.meth.stg.2 but for stage 3.
#' @param switch.n An integer defining the number of changes which are made in each iteration of the the stage 3 iterative improvement procedure. Recommendation is to leave at 1.
#' @param iter.max.stg1 The maximum number of iterations with no improvement during stage 1 before initial point generation is restarted.
#' @param iter.max.stg2 The maximum number of iterations with no improvement during stage 2 before initial point generation is restarted.
#' @param iter.max.stg3 The maximum number of iterations for the the stage 3 iterative improvement procedure before the algorithm is stopped if div.n.flat has not been reached.
#' @param div.int Integer indicating the interval at which to record divergence values during the stage 3 iterative improvement procedure, so that the trace can be examined to determine if it has minimised.
#' @param break.num Integer specifying the number of breaks (to delineate bins in a histogram) for the KL calculation. Defaults to 20, but this is arbitrary.
#' @param ret.seed.pts Logical indicating whether to output the set of 'seed points' produced in stage 2. Useful for illustrating the advantage of the stage 3 iterative improvement procedure.
#' @param ret.all.iter Logical indicating whether to retain the points from every iteration of stage 3. Can be used to examine the stage 3 iterative improvement procedure in detail.
#' @param logfile A string with a filename to output progress information to, if NULL it is printed to the console. Useful if running multiple runs in parallel to prevent progress info from multiple runs being mixed.
#' @param verbose Logical indicating whether to print progress information.
#' @return Returns a list with the following elements:
#' \itemize{
#' \item{points:}{ the simulated points}
#' \item{div.vecs:}{ a list containing vectors of the divergence metric sampled every div.int iterations. There is one for each distance measure (i.e. intraspecific distances for each species as well as interspecific distances if inter.spp == TRUE) and weighted mean divergence across all distance measures.}
#' \item{dist.obs:}{ a list containing vectors of interpoint distances for each distance measure for the observed points}
#' \item{dist.sim:}{ a list containing vectors of interpoint distances for each distance measure for the simulated points, can be compared to dist.obs}
#' \item{seed.pts:}{ the seed points produced in stage 2, only outputted if ret.seed.pts == TRUE}
#' \item{pts.progress:}{ a list containing every iteration of points from stage 3. Only outputted if ret.all.iter == TRUE}
#' \item{Nflat.vec:}{ a vector containing Nflat - the number of iterations since an improvement has been made - for every iteration of stage 3. Only outputted if ret.all.iter == TRUE}
#' }
#' @importFrom geosphere distCosine distGeo distHaversine distm distVincentyEllipsoid distVincentySphere
#' @importFrom graphics hist
#' @importFrom stats density quantile rnorm runif
#' @importFrom utils combn
#' @export
#' @examples
#' \dontrun{
#' # intraspecific distances only
#' my.sim <- fauxcurrence(coords=my.coords,rast=my.raster)
#' # intraspecific distances and interspecific distances
#' my.sim <- fauxcurrence(coords=my.coords,rast=my.raster,inter.spp=TRUE)
#' }
#' @export
fauxcurrence <-function(coords,rast,distmat=NULL,use.distmat=FALSE,inter.spp=FALSE,sep.inter.spp=FALSE,fix.seed.pts=NULL,div.n.flat=1000,allow.ident.conspec=FALSE,dist_meth="distRcpp",dist_fun="Haversine",trim.init.range=FALSE,init.range=0.9,new.pt.meth.stg.2="dist.dens",new.pt.meth.stg.3="dist.dens",switch.n=1,iter.max.stg1=10000,iter.max.stg2=2000,iter.max.stg3=100000,div.int=10,break.num=20,ret.seed.pts=FALSE,ret.all.iter=FALSE,logfile=NULL,verbose=TRUE) {
  # record start time
  start_time <- Sys.time()
  # start outputting to logfile if provided
  if(!is.null(logfile)){
    sink(logfile)
  }
  on.exit(if(!is.null(logfile)) sink())
  # check inputs and options:
  ## coords
  if(class(coords) != "data.frame") stop("coords must be a dataframe")
  if(!(all(c("x","y") %in% colnames(coords))) | ncol(coords)>3) stop("coords must have columns 'x' (longitude) and 'y' (latitude) and (optionally) 'species'")
  if(ncol(coords) == 3){
    if(min(table(coords$species)) < 3) stop("species with fewer than three occurrences are currently unsupported")
  }
  ## rast
  if(class(rast) != "RasterLayer") stop("rast must be a RasterLayer")
  ## distmat
  if(!is.null(distmat)){
    if(class(distmat)[1] != "matrix") stop("distmat must be a matrix")
    # check sizes of rast and distmat match
    if(use.distmat & nrow(distmat) != length(which(!is.na(raster::values(rast))))) {
      stop("number of non-NA cells in rast and number of columns/rows in distmat do not match")
    }
  }
  ## dist_meth
  if(!(dist_meth %in% c("distm","costdist","distRcpp"))) stop('dist_meth must be either "distm", "distRcpp" or "costdist"')
  ## init_range
  if(init.range <= 0 | init.range > 1) stop("init.range must be between 0 and 1")
  ## fix.seed.pts
  if(!is.null(fix.seed.pts)){
    if(class(fix.seed.pts) != "data.frame") stop("fix.seed.pts must be a dataframe")
    if(!(all(c("x","y","species") %in% colnames(fix.seed.pts))) | ncol(fix.seed.pts)!=3) stop("fix.seed.pts must have columns 'x' (longitude) and 'y' (latitude) and 'species'")
    if(any(sort(fix.seed.pts[,"species"]) != sort(unique(coords[,"species"])))) stop("fix seed points must have columns 'x' (longitude), 'y' (latitude) and 'species' and contain one point for each unique species in coords")
  }
  # If verbose, print message with options
  if(verbose){
    cat("Starting fauxcurrence\n")
    cat("Input options:\n")
    cat("   coords:",deparse(substitute(coords)),"\n")
    cat("   rast:",deparse(substitute(rast)),"\n")
    cat("   distmat:",deparse(substitute(distmat)),"\n")
    cat("   fix.seed.pts:",deparse(substitute(fix.seed.pts)),"\n\n")
    cat("Main algorithm options:\n")
    cat("   inter.spp:",inter.spp,"\n")
    cat("   sep.inter.spp:",sep.inter.spp,"\n")
    cat("   div.n.flat:",div.n.flat,"\n")
    cat("   iter.max.stg3:",iter.max.stg3,"\n")
    cat("   use.distmat:",use.distmat,"\n")
    cat("   dist_meth:",dist_meth,"\n")
    cat("   dist_fun:",dist_fun,"\n")
    cat("   allow.ident.conspec:",allow.ident.conspec,"\n\n")
    cat("Advanced and tuning options:\n")
    cat("   new.pt.meth.stg.2:",new.pt.meth.stg.2,"\n")
    cat("   new.pt.meth.stg.3:",new.pt.meth.stg.3,"\n")
    cat("   iter.max.stg1:",iter.max.stg1,"\n")
    cat("   iter.max.stg2:",iter.max.stg2,"\n")
    cat("   switch.n:",switch.n,"\n")
    cat("   trim.init.range:",trim.init.range,"\n")
    cat("   init.range:",init.range,"\n")
    cat("   break.num:",break.num,"\n\n")
    cat("Output options:\n")
    cat("   div.int:",div.int,"\n")
    cat("   ret.seed.pts:",ret.seed.pts,"\n")
    cat("   ret.all.iter:",ret.all.iter,"\n")
    cat("   logfile:",deparse(substitute(logfile)),"\n")
    cat("   verbose:",verbose,"\n\n")
  }
  # Initiate output
  out <- list()
  # if use.distmat is TRUE and distmat is NULL, create distance matrix.
  if(is.null(distmat) & use.distmat){
    distmat <- make.distmat(rast = rast, dist_meth = dist_meth, dist_fun = dist_fun)
  }
  # make dummy species vector if there is no species column in coords, otherwise take species column, also make vector of unique species
  if(!("species" %in% colnames(coords))){
    species.vec <- rep("sp1",nrow(coords))
  } else {
    species.vec <- as.character(coords[,"species"])
  }
  coords=as.matrix(coords[,c("x","y")])
  species <- sort(unique(species.vec))
  # check species names
  if(length(grep("_",species)) > 0 | length(grep(" ",species)) > 0) stop("species names cannot contain the characters '_' or ' '")
  if(length(species) == 1 & any(inter.spp,sep.inter.spp)){
    inter.spp <- FALSE
    sep.inter.spp <- FALSE
    if(verbose) cat("warning: automatically setting inter.species and sep.inter.species to FALSE because only one species is present.\n")
  }
  # set sep.inter.spp to FALSE if there are too few species
  if(length(species) == 2){
    sep.inter.spp <- TRUE
    if(verbose) cat("warning: automatically setting sep.inter.species to TRUE because only two species are present. This prevents two identical vectors of interspecific distances being created.\n")
  }
  # if sep.inter.spp = TRUE, get all combinations of species
  if(sep.inter.spp){
    combs <- combn(species,m=2,simplify=T)
    colnames(combs) <- paste("inter",combs[1,],combs[2,],sep="_")
  } else {
    combs <- NULL
  }
  # get observed distances using the appropriate method. If use.distmat is TRUE, also make named vectors to convert between indices of raster cells and distance matrix.
  if(use.distmat){
    all.nonNA <- 1:length(raster::Which(!is.na(rast),cells=T))
    all.nonNA.rev <- raster::Which(!is.na(rast),cells=T)
    names(all.nonNA) <-  as.character(all.nonNA.rev)
    names(all.nonNA.rev) <-  as.character(unname(all.nonNA))
    dist.orig <- get_dists_from_mat(mat = distmat,indices = unname(all.nonNA[as.character(raster::cellFromXY(rast,coords))]),species.vec=species.vec,inter.spp=inter.spp,sep.inter.spp=sep.inter.spp,combs=combs)
  } else {
    dist.orig <- get_dists(coords=coords,species.vec=species.vec,dist_meth=dist_meth,dist_fun=dist_fun,inter.spp=inter.spp,sep.inter.spp=sep.inter.spp,combs=combs)
  }
  out[["dist.obs"]] <-  dist.orig
  # if either new.pt.meth is "dist.dens", make density objects for each sample
  if(new.pt.meth.stg.2 == "dist.dens" | new.pt.meth.stg.3 == "dist.dens"){
    dist.dens.orig <- list()
    for(n in names(dist.orig)){
      dist.dens.orig[[n]] <- density(dist.orig[[n]])
    }
  }
  # get lists of min and max distance for each distance measure
  min.list <- list()
  max.list <- list()
  for (dis in names(dist.orig)){
      min.list[[dis]] <- min(dist.orig[[dis]])[1]
      max.list[[dis]] <- max(dist.orig[[dis]])[1]
  }
  # if trim.init.range is TRUE, get lists of lower and upper quantiles to constrain initial points to.
  if(trim.init.range){
    quant.min <- list()
    quant.max <- list()
    for (dis in names(dist.orig)){
      quant.min[[dis]] <- quantile(dist.orig[[dis]],0.5-init.range/2)
      quant.max[[dis]] <- quantile(dist.orig[[dis]],0.5+init.range/2)
    }
  }
  # get initial set of points
  if(verbose){
    cat("Starting seed point generation...")
    seed.start <- Sys.time()
  }
  # initiate counters for number of restarts in stage 1 and 2
  n.restart.1 <- 0
  n.restart.2 <- 0
  repeat{
    ######## STAGE 1
    # generate initial points randomly if fix.seed.pts is NULL, or use fix.seed.pts if supplied
    if(is.null(fix.seed.pts)){
      if(use.distmat){
        pts <- unname(sample(all.nonNA,length(species)))
      } else {
        pts <- raster::xyFromCell(rast,sample(which(!is.na(raster::values(rast))),length(species)))
      }
      pts.species <- species
      cnt.wrong.old <- NULL
      iter.unchanged <- 0
      # this loop replaces the starting points (one for each species), until the interpoint distances are within the range of the real data (and within the middle <init.range> quantile if trim.init.range is TRUE)
      repeat{
        # if the number of iterations since a change was last made (iter.unchanged) is above the maximum limit (iter.max.stg1), start again with a new set of random points. This prevents the algorithm getting stuck with an initial set of points which can't be made to sit within the real distribution through stepwise changes.
        if(iter.unchanged > iter.max.stg1){
          if(use.distmat){
            pts <- unname(sample(all.nonNA,length(species)))
          } else {
            pts <- raster::xyFromCell(rast,sample(which(!is.na(raster::values(rast))),length(species)))
          }
          if(verbose) n.restart.1 <- n.restart.1+1
          cnt.wrong.old <- NULL
          iter.unchanged <- 0
        }
        pts.new <- pts
        # unless this is the first iteration, replace a point
        if(!is.null(cnt.wrong.old)){
          if(use.distmat){
            old.pt <- sample(length(pts),1)
            new.pt <- unname(sample(all.nonNA,1))
            pts.new[old.pt] <- new.pt
          } else {
            old.pt <- sample(nrow(pts),1)
            new.pt <- raster::xyFromCell(rast,sample(which(!is.na(raster::values(rast))),1))
            pts.new[old.pt,] <- new.pt
          }
        }
        # get distances
        if(use.distmat){
          dist.pts <- get_dists_from_mat(mat=distmat,indices=pts.new,species.vec=pts.species,inter.spp=inter.spp,sep.inter.spp=sep.inter.spp,combs=combs)
        } else {
          dist.pts <- get_dists(coords=pts.new,species.vec=pts.species,dist_meth=dist_meth,dist_fun=dist_fun,inter.spp=inter.spp,sep.inter.spp=sep.inter.spp,combs=combs)
        }
        # count the distance measures which are outside the range of the real distances (or outside the trimmed range if trim.init.range = TRUE).
        cnt.wrong <- 0
        for(dis in names(dist.pts)){
          if(length(dist.pts[[dis]]) > 0){
            if(trim.init.range){
              if(min(dist.pts[[dis]]) < quant.min[[dis]] | max(dist.pts[[dis]]) > quant.max[[dis]]){
                cnt.wrong <- cnt.wrong+1
              }
            } else {
              if(min(dist.pts[[dis]]) < min.list[[dis]] | max(dist.pts[[dis]]) > max.list[[dis]]){
                cnt.wrong <- cnt.wrong+1
              }
            }
          }
        }
        # if all distances are within the real range (or trimmed real range if trim.init.range = TRUE), end loop
        if(cnt.wrong == 0){
          pts <- pts.new
          break
        }
        # if this is the first iteration set count.wrong.old
        if(is.null(cnt.wrong.old)){
          cnt.wrong.old <- cnt.wrong
        }
        # if more distances are within the real range than in the previous iteration, retain the change and reset iter.unchanged. If not, keep old points and add one to iter.unchanged.
        if(cnt.wrong < cnt.wrong.old){
          cnt.wrong.old <- cnt.wrong
          pts <- pts.new
          iter.unchanged <- 0
        } else {
          iter.unchanged <- iter.unchanged+1
        }
      }
    } else {
      if(use.distmat){
        pts <- unname(all.nonNA[as.character(raster::cellFromXY(rast,as.matrix(fix.seed.pts[,c("x","y")])))])
      } else {
        pts <- as.matrix(fix.seed.pts[,c("x","y")])
      }
      pts.species <- as.character(fix.seed.pts[,"species"])
      if(verbose) cat("Stage 1: initial seed points taken from fix.seed.pts")
    }
    ######## STAGE 1 END
    ########  STAGE 2
    # This loop adds new points until the correct number of points for each species has been generated
    # get number of points for each species
    len.sp <- list()
    for(sp in species){
      len.sp[[sp]] <- length(which(species.vec == sp))
    }
    # make pts and pts.species vectors the right size
    if(use.distmat){
      pts <- c(pts,rep(NA,nrow(coords)-length(pts)))
    } else {
      pts <- rbind(pts,matrix(NA,nrow=nrow(coords)-nrow(pts),ncol=2))
    }
    pts.species <- c(pts.species,rep(NA,nrow(coords)-length(pts.species)))
    iter.unchanged <- 0
    restart<- FALSE
    while(length(which(is.na(pts.species))) > 0){
      # if maximum iterations for stage 2 has been reached, break the loop so the algorithm can restart initial point generation
      if(iter.unchanged > iter.max.stg2){
        restart <- TRUE
        break
      }
      #sample random species
      sp <- sample(species,1)
      # skip if enough points have already been generated for the species
      if(length(which(pts.species == sp)) == len.sp[[sp]]) next
      # for "dist.obs" or "dist.dens" get an orig.pt and distance
      if(new.pt.meth.stg.2 == "dist.obs" | new.pt.meth.stg.2 == "dist.dens"){
        # sample a random starting point from pts, if there is only one then use that rather than attempting to pick one with sample.
        if(use.distmat){
          if(length(which(pts.species == sp)) == 1){
            orig.pt <- pts[which(pts.species == sp)]
          } else {
            orig.pt <- pts[sample(which(pts.species == sp),1)]
          }
        } else {
          if(length(which(pts.species == sp)) == 1){
            orig.pt <- pts[which(pts.species == sp),]
          } else {
            orig.pt <- pts[sample(which(pts.species == sp),1),]
          }
        }
        # if not using distmat, sample random bearing
        if(!use.distmat) bear<-runif(1,0,360)
        # sample intraspecific distance from real intraspecific distances, either from the vector of distances for "dist.obs", or from the density object for "dist.dens"
        if(new.pt.meth.stg.2 == "dist.obs"){
          new.dist <- sample(dist.orig[[paste("intra",sp,sep="_")]],1)
        } else if (new.pt.meth.stg.2 == "dist.dens"){
          new.dist <- sample(dist.orig[[paste("intra",sp,sep="_")]],1) + rnorm(1,0,dist.dens.orig[[paste("intra",sp,sep="_")]]$bw)
        }
        # generate new point
        if(use.distmat){
          # generate a new point by finding the point in distmat whos distance from orig.pt is most similar to new.dist, ties are broken at random.
          new.pt <- max.col(matrix(0-abs(new.dist - distmat[orig.pt,]),nrow=1))
        } else {
          new.pt <- geosphere::destPoint(p=orig.pt,b=bear,d=(new.dist))
        }
        # if new.pt.meth.stg.2 == "sample", just sample a random point from the raster
      } else if (new.pt.meth.stg.2 == "sample"){
        if(use.distmat){
          new.pt <-  unname(sample(all.nonNA,1))
        } else {
          new.pt <- raster::xyFromCell(rast,sample(which(!is.na(raster::values(rast))),1))
        }
      }
      # If using distance matrix and if allow.ident.conspec == FALSE, check new.pt is not the same as any conspecific points, otherwise skip the iteration. If not using distance matrix, check new point is within the study area and if allow.ident.conspec == FALSE, check new.pt is not the same as any conspecific points, otherwise skip the iteration.
      if(use.distmat){
        if(new.pt %in% pts[which(pts.species == sp)] & !(allow.ident.conspec)) next
      } else {
        if(is.na(raster::extract(rast,new.pt))) next
        if(nrow(unique(rbind(pts[which(pts.species == sp),],new.pt),1,incomparables = F)) != nrow(rbind(pts[which(pts.species == sp),],new.pt)) & !(allow.ident.conspec)) next
      }
      # find next NA point and assign new point to it
      if(use.distmat){
        next.pt <- which(is.na(pts))[1]
        pts[next.pt] <- new.pt
      } else {
        next.pt <- which(is.na(pts[,1]))[1]
        pts[next.pt,] <- new.pt
      }
      pts.species[next.pt] <- sp
      # calculate distances
      if(use.distmat){
        dist.pts <- get_dists_from_mat(mat = distmat,indices = pts[which(!is.na(pts))], species.vec=pts.species[which(!is.na(pts.species))],inter.spp=inter.spp,sep.inter.spp=sep.inter.spp,combs=combs)
      } else {
        dist.pts <- get_dists(coords=pts[which(!is.na(pts[,1])),],species.vec=pts.species[which(!is.na(pts.species))],dist_meth=dist_meth,dist_fun=dist_fun,inter.spp=inter.spp,sep.inter.spp=sep.inter.spp,combs=combs)
      }
      # if any of the distances are outside the range of the real distances, remove the new points
      dist.wrong <- FALSE
      for(dis in names(dist.pts)){
        if(length(dist.pts[[dis]]) > 0){
         if(min(dist.pts[[dis]]) < min.list[[dis]] | max(dist.pts[[dis]]) > max.list[[dis]]){
           dist.wrong <- TRUE
         }
        }
      }
      if(dist.wrong){
        if(use.distmat){
          pts[next.pt] <- NA
        } else {
          pts[next.pt,] <- NA
        }
        pts.species[next.pt] <- NA
        iter.unchanged <- iter.unchanged+1
      } else {
        iter.unchanged <- 0
      }
    }
    # if the loop has ended because the maximum number of iterations for stage 2 has been reached (i.e. restart == TRUE), restart the loop, otherwise it must have ended because a full set of points has been created, so break the loop
    if(restart){
      if(verbose) n.restart.2 <- n.restart.2+1
      next
    } else {
      break
    }
  }
  if(verbose){
    seed.time <- difftime(Sys.time(),seed.start,units="secs")
    cat("completed in ",round(as.numeric(seed.time),2)," ",units(seed.time),"\n","Seed point generation restarted ",n.restart.1," times at stage 1 and ",n.restart.2," times at stage 2.\nIf the number of restarts is high, consider adjusting iter.max.stg1 and iter.max.stg2\n\n",sep="")
  }
  ######## STAGE 2 END
  # get distances of simulated points
  if(use.distmat){
    dist.pts <- get_dists_from_mat(mat=distmat,indices=pts,species.vec=pts.species,inter.spp=inter.spp,sep.inter.spp=sep.inter.spp,combs=combs)
  } else {
    dist.pts <- get_dists(coords=pts,species.vec=pts.species,dist_meth=dist_meth,dist_fun=dist_fun,inter.spp=inter.spp,sep.inter.spp=sep.inter.spp,combs=combs)
  }
  # get initial divergence between distances for observed and simulated points and breaks
  div.curr <- rep(NA,length(dist.orig))
  names(div.curr) <- names(dist.orig)
  orig.breaks <- list()
  for(n in names(dist.orig)){
    div <- KL.div(dist.orig[[n]],dist.pts[[n]],break.num=break.num)
    div.curr[[n]] <- div[[1]]
    orig.breaks[[n]] <- div[[2]]
  }
  div.curr[["mean"]] <- mean(c(mean(div.curr[grep("inter_",names(div.curr))]),mean(div.curr[grep("intra_",names(div.curr))])),na.rm=T)
  # set up intervals for recording divergence and points
  div.int.seq<-seq(0,iter.max.stg3,div.int)
  div.int.seq[1]<-1
  # initiate outputs for recording progress of divergence and points
  div.vecs <- list()
  for(n in names(dist.orig)){
    div.vecs[[n]] <- rep(NA,length(div.int.seq))
    div.vecs[[n]][1] <- div.curr[[n]]
  }
  div.vecs[["mean"]] <- rep(NA,length(div.int.seq))
  div.vecs[["mean"]][1] <- div.curr[["mean"]]
  # If ret.all.iter and/or ret.seed.points is TRUE, save seed points. If use.distmat is TRUE, convert them back to coordinates first. If ret.all.iter is TRUE, also set up a vector to record Nflat.
  if(ret.all.iter | ret.seed.pts){
    if(ret.all.iter) pts.progress <- list()
    if(use.distmat){
      new.coords <- as.data.frame(raster::xyFromCell(rast,all.nonNA.rev[pts]))
      new.coords$species <- pts.species
      if(ret.all.iter) pts.progress[[1]] <- new.coords
      if(ret.seed.pts) out[["seed.pts"]] <- new.coords
    } else {
      if(ret.all.iter) pts.progress[[1]] <- cbind(as.data.frame(pts),species=pts.species)
      if(ret.seed.pts) out[["seed.pts"]] <- cbind(as.data.frame(pts),species=pts.species)
    }
    if(ret.all.iter){
      Nflat.vec <- rep(NA,length(div.int.seq))
      Nflat.vec[[1]] <- 0
    }
  }
  ######## STAGE 3
  # set up counters
  int.num <- 2
  L <- 1
  Nflat <- 0
  if(verbose){
    cat("Starting stage 3: iterative improvement...")
    stg.3.start <- Sys.time()
    progress.seq <- seq(0,iter.max.stg3,iter.max.stg3/10)
    progress.seq <- progress.seq[2:(length(progress.seq)-1)]
  }
  # This loop replaces points, checks if the change improves the fit of interpoint distances to observed interpoint distances, and if so keeps the change. It runs until no improvements have been made for Nflat iterations or iter.max.stg3 has been reached.
  while(L <= iter.max.stg3 & Nflat <= div.n.flat){
    if(verbose){
      if(L %in% progress.seq) cat(L,"...",sep="")
    }
    # make new points
    new.pts<-pts
    # get random points to replace, do not replace seeded points if used
    if(is.null(fix.seed.pts)){
      if(use.distmat){
        to.replace<-sample(length(pts),switch.n)
      } else {
        to.replace<-sample(nrow(pts),switch.n)
      }
    } else {
      if(use.distmat){
        to.replace<-sample((nrow(fix.seed.pts)+1):length(pts),switch.n)
      } else {
        to.replace<-sample((nrow(fix.seed.pts)+1):nrow(pts),switch.n)
      }
    }
    # this loop replaces each point in to.replace with a new point
    for (P in 1:switch.n) {
      # for "dist.obs" or "dist.dens" get an orig.pt and distance
      if(new.pt.meth.stg.3 == "dist.obs" | new.pt.meth.stg.3 == "dist.dens"){
        # get species
        my.sp <- pts.species[to.replace[P]]
        # sample a random starting point from the correct species
        if(use.distmat){
          orig.pt <- pts[sample(which(pts.species == my.sp),1)]
        } else {
          orig.pt <- pts[sample(which(pts.species == my.sp),1),]
        }
        # if not using distmat, sample random bearing
        if(!use.distmat) bear<-runif(1,0,360)
        # sample intraspecific distance from real intraspecific distances, either from the vector of distances for "dist.obs", or from the density object for "dist.dens"
        if(new.pt.meth.stg.3 == "dist.obs"){
          new.dist <- sample(dist.orig[[paste("intra",my.sp,sep="_")]],1)
        } else if (new.pt.meth.stg.3 == "dist.dens"){
          new.dist <- sample(dist.orig[[paste("intra",my.sp,sep="_")]],1) + rnorm(1,0,dist.dens.orig[[paste("intra",my.sp,sep="_")]]$bw)
        }
        # generate new point
        if(use.distmat){
          # generate a new point by finding the point in distmat whos distance from orig.pt is most similar to new.dist, ties are broken at random.
          new.pt <- max.col(matrix(0-abs(new.dist - distmat[orig.pt,]),nrow=1))
        } else {
          new.pt <- geosphere::destPoint(p=orig.pt,b=bear,d=new.dist)
        }
        # if new.pt.meth.stg.3 == "sample", just sample a random point from the raster
      } else if (new.pt.meth.stg.3 == "sample"){
        if(use.distmat){
          new.pt <-  unname(sample(all.nonNA,1))
        } else {
          new.pt <- raster::xyFromCell(rast,sample(which(!is.na(raster::values(rast))),1))
        }
      }
      # If using distance matrix and if allow.ident.conspec == FALSE, check new.pt is not the same as any conspecific points, otherwise skip the iteration. If not using distance matrix, check new point is within the study area and if allow.ident.conspec == FALSE, check new.pt is not the same as any conspecific points, otherwise skip the iteration.
      if(use.distmat){
        if(new.pt %in% pts[which(pts.species == sp)] & !(allow.ident.conspec)) next
      } else {
        if(is.na(raster::extract(rast,new.pt))) next
        if(nrow(unique(rbind(pts[which(pts.species == sp),],new.pt),1,incomparables = F)) != nrow(rbind(pts[which(pts.species == sp),],new.pt)) & !(allow.ident.conspec)) next
      }
      if(use.distmat){
        new.pts[to.replace[P]] <- new.pt
      } else {
        new.pts[to.replace[P],] <- new.pt
      }
      # get new distances
      if(use.distmat){
        dist.pts <- get_dists_from_mat(mat=distmat,indices=new.pts,species.vec=pts.species,inter.spp=inter.spp,sep.inter.spp=sep.inter.spp,combs=combs)
      } else {
        dist.pts <- get_dists(coords=new.pts,species.vec=pts.species,dist_meth=dist_meth,dist_fun=dist_fun,inter.spp=inter.spp,sep.inter.spp=sep.inter.spp,combs=combs)
      }
      # if any of the distances are outside the range of the real distances, remove the new point
      for(dis in names(dist.pts)){
        if(min(dist.pts[[dis]]) < min.list[[dis]] | max(dist.pts[[dis]]) > max.list[[dis]]){
          if(use.distmat){
            new.pts[to.replace[P]] <- pts[to.replace[P]]
          } else {
            new.pts[to.replace[P],] <- pts[to.replace[P],]
          }
        }
      }
    } #P loop
    # calculate new divergence. If divergence fails, skip iteration
    div.new <- rep(NA,length(dist.orig))
    names(div.new) <- names(dist.orig)
    div.fail <- FALSE
    for(n in names(dist.orig)){
      div <- try(KL.div(dist.orig[[n]],dist.pts[[n]],break.num=break.num),silent = T)
      if (class(div)=="try-error"){
        div.fail <- TRUE
      } else {
        div.new[[n]] <- div[[1]]
      }
    }
    if(div.fail){
      if(L%in%div.int.seq){
        for(n in names(dist.orig)){
          div.vecs[[n]][int.num] <- div.curr[[n]]
        }
        div.vecs[["mean"]][int.num] <- div.curr[["mean"]]
        if(ret.all.iter){
          if(use.distmat){
            pts.progress[[int.num]] <- cbind(as.data.frame(raster::xyFromCell(rast,all.nonNA.rev[pts])),species=pts.species)
          } else {
            pts.progress[[int.num]] <- cbind(as.data.frame(pts),species=pts.species)
          }
          Nflat.vec[int.num] <- Nflat
        }
        int.num <- int.num+1
      }
    } else {
      div.new[["mean"]] <- mean(c(mean(div.new[grep("inter_",names(div.new))]),mean(div.new[grep("intra_",names(div.new))])),na.rm=T)
      if(div.new[["mean"]] < div.curr[["mean"]]){
        pts <- new.pts
        div.curr <- div.new
      }
      if(L%in%div.int.seq){
        for(n in names(dist.orig)){
          div.vecs[[n]][int.num] <- div.curr[[n]]
        }
        div.vecs[["mean"]][int.num] <- div.curr[["mean"]]
        if(ret.all.iter){
          if(use.distmat){
            pts.progress[[int.num]] <- cbind(as.data.frame(raster::xyFromCell(rast,all.nonNA.rev[pts])),species=pts.species)
          } else {
            pts.progress[[int.num]] <- cbind(as.data.frame(pts),species=pts.species)
          }
          Nflat.vec[int.num] <- Nflat
        }
        int.num <- int.num+1
      }
    }
    # update Nflat, div.old and L
    if(L > 1){
      if(div.curr[["mean"]] == div.old){
        Nflat <- Nflat+1
      } else {
        Nflat <- 0
      }
    }
    div.old <- div.curr[["mean"]]
    L <- L+1
  } # L loop
  if(verbose){
    stg.3.time <- difftime(Sys.time(),stg.3.start,units="secs")
    cat("completed in",round(as.numeric(stg.3.time),2),units(stg.3.time),"with",L-1,"iterations","\n\n")
    txtplot::txtplot(1:length(div.vecs[["mean"]])*div.int,div.vecs[["mean"]],height = 12,width=50,xlab = "N iterations",ylab="D")
    cat("\n")
  }
  if(verbose & Nflat <= div.n.flat){
    cat("warning: div.n.flat was not reached, please check D is minimised and consider increasing iter.max.stg3\n")
  }
  # trim div.vec and assign to output list
  for(n in names(div.vecs)){
    div.vecs[[n]] <- div.vecs[[n]][!is.na(div.vecs[[n]])]
  }
  out[["div.vecs"]] <- div.vecs
  # if ret.all.iter is TRUE, trim Nflat.vec and assign to output list
  if(ret.all.iter){
    Nflat.vec <- Nflat.vec[!is.null(Nflat.vec)]
    out[["Nflat.vec"]] <- Nflat.vec
    out[["pts.progress"]] <- pts.progress
  }
  # assign simulated distances to output list
  out[["dist.sim"]] <- dist.pts
  # assign final points to output list
  if(use.distmat){
    out[["points"]] <- cbind(as.data.frame(raster::xyFromCell(rast,all.nonNA.rev[pts])),species=pts.species)
  } else {
    out[["points"]] <- new.coords <- cbind(as.data.frame(pts),species=pts.species)
  }
  if(verbose){
    total.time <- difftime(Sys.time(),start_time,units="secs")
    cat("Total analysis completed in",round(as.numeric(total.time),2),units(total.time),"\n")
  }
  return(out)
} # end

