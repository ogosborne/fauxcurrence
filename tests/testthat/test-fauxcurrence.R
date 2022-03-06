library(raster)
library(sp)
library(geosphere)

test_that("fauxcurrence produces the same output with the same random seed, regardless of output settings", {
  ## run function
  # intra model
  coords <- chameleons[1:3,]
  set.seed(1234)
  intra <- fauxcurrence(coords = coords, rast = madagascar,iter.max.stg3 = 10, verbose = F)
  set.seed(1234)
  capture.output(
    intra_v <- fauxcurrence(coords = coords, rast = madagascar,iter.max.stg3 = 10, verbose = T)
  )
  set.seed(1234)
  intra_rsp <- fauxcurrence(coords = coords, rast = madagascar,iter.max.stg3 = 10, ret.seed.pts = T, verbose = F)
  set.seed(1234)
  intra_rai <- fauxcurrence(coords = coords, rast = madagascar,iter.max.stg3 = 10, ret.all.iter = T, verbose = F)
  file <- tempfile()
  set.seed(1234)
  intra_lf <- fauxcurrence(coords = coords, rast = madagascar,iter.max.stg3 = 10, logfile = file, verbose = F)
  # inter model
  coords <- chameleons[c(1:3,9:11,61:63),]
  set.seed(1234)
  inter <- fauxcurrence(coords = coords, rast = madagascar,iter.max.stg3 = 10, inter.spp = TRUE, verbose = F)
  set.seed(1234)
  capture.output(
    inter_v <- fauxcurrence(coords = coords, rast = madagascar,iter.max.stg3 = 10, inter.spp = TRUE, verbose = T)
  )
  set.seed(1234)
  inter_rsp <- fauxcurrence(coords = coords, rast = madagascar,iter.max.stg3 = 10, ret.seed.pts = T, inter.spp = TRUE, verbose = F)
  set.seed(1234)
  inter_rai <- fauxcurrence(coords = coords, rast = madagascar,iter.max.stg3 = 10, ret.all.iter = T, inter.spp = TRUE, verbose = F)
  file <- tempfile()
  set.seed(1234)
  inter_lf <- fauxcurrence(coords = coords, rast = madagascar,iter.max.stg3 = 10, logfile = file, inter.spp = TRUE, verbose = F)
  set.seed(1234)
  # interSep model
  coords <- chameleons[c(1:3,9:11),]
  set.seed(1234)
  interS <- fauxcurrence(coords = coords, rast = madagascar,iter.max.stg3 = 10, inter.spp = TRUE, sep.inter.spp = TRUE, verbose = F)
  set.seed(1234)
  capture.output(
    interS_v <- fauxcurrence(coords = coords, rast = madagascar,iter.max.stg3 = 10, inter.spp = TRUE, sep.inter.spp = TRUE, verbose = T)
  )
  set.seed(1234)
  interS_rsp <- fauxcurrence(coords = coords, rast = madagascar,iter.max.stg3 = 10, ret.seed.pts = T, inter.spp = TRUE, sep.inter.spp = TRUE, verbose = F)
  set.seed(1234)
  interS_rai <- fauxcurrence(coords = coords, rast = madagascar,iter.max.stg3 = 10, ret.all.iter = T, inter.spp = TRUE, sep.inter.spp = TRUE, verbose = F)
  file <- tempfile()
  set.seed(1234)
  interS_lf <- fauxcurrence(coords = coords, rast = madagascar,iter.max.stg3 = 10, logfile = file, inter.spp = TRUE, sep.inter.spp = TRUE, verbose = F)
  # test that outputs that only differ by output options are identical
  expect_identical(intra$points, intra_v$points)
  expect_identical(intra$points, intra_rsp$points)
  expect_identical(intra$points, intra_rai$points)
  expect_identical(intra$points, intra_lf$points)
  expect_identical(inter$points, inter_v$points)
  expect_identical(inter$points, inter_rsp$points)
  expect_identical(inter$points, inter_rai$points)
  expect_identical(inter$points, inter_lf$points)
  expect_identical(interS$points, interS_v$points)
  expect_identical(interS$points, interS_rsp$points)
  expect_identical(interS$points, interS_rai$points)
  expect_identical(interS$points, interS_lf$points)
})

test_that("make.distmat makes matrix of correct size", {
  rast <- raster::crop(madagascar,sp::SpatialPoints(chameleons[1:5,1:2]))
  dm <- make.distmat(rast)
  expect_equal(ncol(dm), dim(rast)[1] * dim(rast)[2])
  expect_equal(nrow(dm), dim(rast)[1] * dim(rast)[2])
  expect_true(is.matrix(dm))
})

test_that("get_dists_coords", {
  res <- get_dists_coords(coords = as.matrix(chameleons[1:5,1:2]), species.vec = chameleons[1:5,3], dist_meth = "distRcpp", dist_fun = "Haversine", inter.spp = F, sep.inter.spp = F, combs = NULL)
  expect_length(res, 1)

})

test_that("get_dists_distmat and get_dists_coords give similar results", {
  # the results are not expected to be identical since get_dists_coords returns distances between the points and get_dists_distmat returns distances between the closest raster cell centers to the points
  pts <- sp::SpatialPoints(chameleons[1:3,1:2], proj4string = crs(madagascar))
  rast <- raster::crop(madagascar,pts)
  dm <- make.distmat(rast, dist_meth = "distRcpp", dist_fun = "Haversine")
  all.nonNA <- get_nonNA(use.distmat = T, rast = rast)[["all.nonNA"]]
  idx <- coords2indices(all.nonNA = all.nonNA, rast = rast, coords = as.matrix(chameleons[1:3,1:2]))
  dist.dm <- get_dists_distmat(mat = dm, indices = idx, species.vec = chameleons[1:3,3], inter.spp = F, sep.inter.spp = F, combs = NULL)
  dist.co <- get_dists_coords(coords = as.matrix(chameleons[1:3,1:2]), species.vec = chameleons[1:3,3], dist_meth = "distRcpp", dist_fun = "Haversine", inter.spp = F, sep.inter.spp = F, combs = NULL)
  expect_true(all(rank(dist.co$intra_antimena) == rank(dist.dm$intra_antimena)))
})


test_that("get.dist.dens returns NULL if neither new.pt.meth option is 'dist.dens'", {
  dist.list <- get_dists(use.distmat = F, distmat = NULL, pts = as.matrix(chameleons[1:8,1:2]), indices = NULL, species.vec = chameleons[1:8,3], inter.spp = F, sep.inter.spp = F, dist_meth = "distRcpp", dist_fun = "Haversine", combs = NULL)
  res <- get.dist.dens(new.pt.meth.stg.2 = "sample", new.pt.meth.stg.3 = "sample", dist.list = dist.list)
  expect_null(res)
})

test_that("get.dist.dens returns a list of density objects if either new.pt.meth option is 'dist.dens'", {
  dist.list <- get_dists(use.distmat = F, distmat = NULL, pts = as.matrix(chameleons[1:8,1:2]), indices = NULL, species.vec = chameleons[1:8,3], inter.spp = F, sep.inter.spp = F, dist_meth = "distRcpp", dist_fun = "Haversine", combs = NULL)
  res1 <- get.dist.dens(new.pt.meth.stg.2 = "dist.dens", new.pt.meth.stg.3 = "sample", dist.list = dist.list)
  res2 <- get.dist.dens(new.pt.meth.stg.2 = "sample", new.pt.meth.stg.3 = "dist.dens", dist.list = dist.list)
  expect_true(is.list(res1))
  expect_true(is.list(res2))
  expect_s3_class(res1[[1]], "density")
  expect_equal(length(dist.list), length(res1))
})

test_that("get_nonNA produces a named list each containing NULL if use.distmat = F", {
  udmF <- get_nonNA(use.distmat = F, rast = madagascar)
  expect_null(udmF[[1]])
  expect_null(udmF[[2]])
  expect_true(is.list(udmF))
  expect_true(all(names(udmF) == c("all.nonNA", "all.nonNA.rev")))
})

test_that("get_nonNA produces a named list of 2 vectors of the right length if use.distmat = T", {
  udmT <- get_nonNA(use.distmat = T, rast = madagascar)
  expect_true(is.list(udmT))
  expect_true(all(names(udmT) == c("all.nonNA", "all.nonNA.rev")))
  expect_equal(length(udmT$all.nonNA), length(udmT$all.nonNA.rev))
  expect_equal(length(udmT$all.nonNA), length(which(!is.na(raster::values(madagascar)))))
})

test_that("get_complete_pts selects correct points", {
  udmT <- get_complete_pts(use.distmat = T, pts = c(1:10,rep(NA,10)))
  udmF <- get_complete_pts(use.distmat = F, pts = matrix(c(1:10,rep(NA,10)), ncol=2, byrow = T))
  expect_equal(length(udmT), 10)
  expect_equal(nrow(udmF), 5)
  expect_false(any(is.na(udmT)))
  expect_false(any(is.na(udmF)))
})

test_that("sample.pts.rand produces a correctly formatted output with the correct number of points", {
  # run function
  udT.10 <- sample.pts.rand(use.distmat = T, all.nonNA = which(!is.na(raster::values(madagascar))), n = 10, rast = madagascar)
  udF.10 <- sample.pts.rand(use.distmat = F, n = 10, rast = madagascar)
  # test types
  expect_type(udT.10, "integer")
  expect_type(udF.10, "double")
  expect_true(is.matrix(udF.10))
  # test colnames
  expect_identical(colnames(udF.10),c("x","y"))
  # test number of points
  expect_equal(nrow(udF.10),10)
  expect_equal(length(udT.10),10)
})

test_that("choose.pts.rand doesn't choose seed points", {
  tab <- table(replicate(choose.pts.rand(pts = as.matrix(chameleons[1:12,1:2]), fix.seed.pts = as.matrix(chameleons[1:10,1:2]),use.distmat = F, n = 1),n=100))
  expect_false(10 %in% names(tab))
})

test_that("get.combs produces null if sep.inter.spp == FALSE", {
  # run function
  sisF <-  get.combs(sep.inter.spp = F, species = letters[1:5])
  # test types
  expect_null(sisF)
})

test_that("get.combs produces a matrix of the correct size if if sep.inter.spp == TRUE", {
  # run function
  sisT5 <- get.combs(sep.inter.spp = T, species = letters[1:5])
  # test types
  expect_true(is.matrix(sisT5))
  expect_type(sisT5, "character")
  # test sizes
  expect_equal(dim(sisT5),c(2,10))
})

test_that("get.dist.minmax produces a list of the correct dimensions", {
  dist.list <- get_dists(use.distmat = F, distmat = NULL, pts = as.matrix(chameleons[1:8,1:2]), indices = NULL, species.vec = chameleons[1:8,3], inter.spp = F, sep.inter.spp = F, dist_meth = "distRcpp", dist_fun = "Haversine", combs = NULL)
  mm <- get.dist.minmax(dist.list)
  expect_length(mm, 2)
  expect_length(mm$min, length(dist.list))
  expect_length(mm$max, length(dist.list))
  expect_gte(mm$max[[1]], mm$min[[1]])
  expect_equal(names(mm$min), names(dist.list))
})

test_that("get.dist.quant produces a minmax list if trim.init.range == FALSE", {
  dist.list <- get_dists(use.distmat = F, distmat = NULL, pts = as.matrix(chameleons[1:8,1:2]), indices = NULL, species.vec = chameleons[1:8,3], inter.spp = F, sep.inter.spp = F, dist_meth = "distRcpp", dist_fun = "Haversine", combs = NULL)
  tirF <- get.dist.quant(trim.init.range = F, dist.list = dist.list)
  mm <- get.dist.minmax(dist.list)
  expect_identical(tirF, mm)
})

test_that("get.dist.quant produces correct results if trim.init.range == TRUE", {
  dist.list <- get_dists(use.distmat = F, distmat = NULL, pts = as.matrix(chameleons[1:8,1:2]), indices = NULL, species.vec = chameleons[1:8,3], inter.spp = F, sep.inter.spp = F, dist_meth = "distRcpp", dist_fun = "Haversine", combs = NULL)
  tirT <- get.dist.quant(trim.init.range = T, dist.list = dist.list, init.range = 0.5)
  expect_equal(tirT$min[[1]],quantile(dist.list[[1]], 0.25))
  expect_equal(tirT$max[[1]],quantile(dist.list[[1]], 0.75))
})

test_that("replace.pt replaces point with NA if new.pt is set to NA", {
  udmT <- replace.pt(use.distmat = T, new.pts = 10:1, new.pt = NA, i = 10)
  udmF <- replace.pt(use.distmat = F, new.pts = matrix(c(1:10), ncol=2, byrow = T), new.pt = NA, i = 4)
  expect_equal(is.na(udmT[10]), TRUE)
  expect_equal(any(is.na(udmT[-10])), FALSE)
  expect_equal(all(is.na(udmF[4,])),TRUE)
  expect_equal(any(is.na(udmF[-4,])),FALSE)
})

test_that("replace.pt replaces point correctly", {
  udmT <- replace.pt(use.distmat = T, new.pts = 10:1, new.pt = 13, i = 10)
  udmF <- replace.pt(use.distmat = F, new.pts = matrix(c(1:10), ncol=2, byrow = T), new.pt = c(13,14), i = 4)
  expect_equal(udmT[10], 13)
  expect_equal(udmT[-10], c(10:1)[-10])
  expect_true(all(udmF[4,] == c(13,14)))
  expect_equal(udmF[-4,], matrix(c(1:10), ncol=2, byrow = T)[-4,])
})

test_that("get.next.pt selects correct point", {
  udmT <- get.next.pt(use.distmat = T, pts = c(1:10,rep(NA,10)))
  udmF <- get.next.pt(use.distmat = F, pts = matrix(c(1:10,rep(NA,10)), ncol=2, byrow = T))
  expect_equal(udmT, 11)
  expect_equal(udmF, 6)
})

test_that("if there is only one point available for a species, sample.orig.pt chooses that point rather than trying to use sample", {
  udmT <- sample.orig.pt(use.distmat = T, pts = 10:1, pts.species = c(rep("a",9),"b"), sp = "b")
  udmF <- sample.orig.pt(use.distmat = F, pts = matrix(c(1:10), ncol=2, byrow = T), pts.species = c(rep("a",4),"b"), sp = "b")
  expect_equal(udmT, 1)
  expect_true(all(udmF ==  c(9,10)))
})

test_that("sample.orig.pt chooses a point of the correct species", {
  udmT <- replicate(sample.orig.pt(use.distmat = T, pts = 10:1, pts.species = rep(letters[1:5], 2), sp = "a"), n = 100)
  udmF <- replicate(sample.orig.pt(use.distmat = F, pts = matrix(c(1:20), ncol=2, byrow = T), pts.species = rep(letters[1:5], 2), sp = "a"), n = 100)
  expect_true(all(sort(unique(udmT)) %in% c(5, 10)))
  expect_true(all(sort(unique(udmF[1,])) %in% c(1, 11)))
  expect_true(all(sort(unique(udmF[2,])) %in% c(2, 12)))
})


test_that("sample.intra.dist samples an observed distance of the correct species when new.pt.meth == 'dist.obs'", {
  dist.list <- get_dists(use.distmat = F, distmat = NULL, pts = as.matrix(chameleons[c(1:3,9:11,61:63),1:2]), indices = NULL, species.vec = chameleons[c(1:3,9:11,61:63),3], inter.spp = F, sep.inter.spp = F, dist_meth = "distRcpp", dist_fun = "Haversine", combs = NULL)
  res <- replicate(sample.intra.dist(new.pt.meth = "dist.obs", dist.list = dist.list, dist.dens = NULL, sp = "oustaleti"), n = 100)
  expect_true(all(res %in% dist.list$intra_oustaleti))
})

test_that("sample.bear returns a single number between 0 and 360 if make.distmat is FALSE", {
  res <- replicate(sample.bear(use.distmat = F), n = 100)
  expect_length(res, 100)
  expect_true(all(res >= 0 & res <= 360))
})

test_that("sample.bear returns NULL if make.distmat is TRUE", {
  udmT <- sample.bear(use.distmat = T)
  expect_null(udmT)
})

test_that("gen.dist.pt produces point the expected distance and bearing away", {
  d <- 10000
  b <- 90
  udmF <- gen.dist.pt(use.distmat = F, distmat = NULL, orig.pt = as.vector(chameleons[1,1:2]), dist = 10000, bear = 90)
  expect_equal(geosphere::distGeo(chameleons[1,1:2], udmF), d)
  expect_equal(geosphere::bearing(chameleons[1,1:2], udmF), b)
})

test_that("check.pt returns FALSE if the point already exists only if allow.ident.conspec == FALSE", {
  udmFaicF <- check.pt(use.distmat = F, allow.ident.conspec = F, new.pt = as.vector(chameleons[1,1:2]), sp = chameleons[1,3], pts = as.matrix(chameleons[,1:2]), pts.species = chameleons$species, rast = madagascar)
  udmFaicT <- check.pt(use.distmat = F, allow.ident.conspec = T, new.pt = as.vector(chameleons[1,1:2]), sp = chameleons[1,3], pts = as.matrix(chameleons[,1:2]), pts.species = chameleons$species, rast = madagascar)
  udmTaicF <- check.pt(use.distmat = T, allow.ident.conspec = F, new.pt = 1, sp = chameleons[1,3], pts = 1:nrow(chameleons), pts.species = chameleons$species, rast = madagascar)
  udmTaicT <- check.pt(use.distmat = T, allow.ident.conspec = T, new.pt = 1, sp = chameleons[1,3], pts = 1:nrow(chameleons), pts.species = chameleons$species, rast = madagascar)
  expect_true(udmFaicT)
  expect_false(udmFaicF)
  expect_true(udmTaicT)
  expect_false(udmTaicF)
})

test_that("check.pt returns TRUE if the point exists, but only for a different species", {
  udmFaicFr <- check.pt(use.distmat = F, allow.ident.conspec = F, new.pt = as.vector(chameleons[1,1:2]), sp = chameleons[1,3], pts = as.matrix(chameleons[,1:2]), pts.species = rev(chameleons$species), rast = madagascar)
  udmTaicFr <- check.pt(use.distmat = T, allow.ident.conspec = F, new.pt = 1, sp = chameleons[1,3], pts = 1:nrow(chameleons), pts.species = rev(chameleons$species), rast = madagascar)
  expect_true(udmFaicFr)
  expect_true(udmTaicFr)
})

test_that("check.pt returns FALSE if the point has value NA", {
  res <- check.pt(use.distmat = F, allow.ident.conspec = F, new.pt = c(43.25 -11.98333), sp = chameleons[1,3], pts = as.matrix(chameleons[,1:2]), pts.species = chameleons$species, rast = madagascar)
  expect_false(res)
})

test_that("check.dists count the correct number of dist sets outside the range", {
  dist.list <- get_dists(use.distmat = F, distmat = NULL, pts = as.matrix(chameleons[,1:2]), indices = NULL, species.vec = chameleons$species, inter.spp = T, sep.inter.spp = T, dist_meth = "distRcpp", dist_fun = "Haversine", combs = get.combs(sep.inter.spp = T, sort(unique(chameleons$species))))
  minmax <- get.dist.minmax(dist.list = dist.list)
  quant <- get.dist.quant(trim.init.range = T, dist.list = dist.list, init.range = 0.5)
  expect_equal(check.dists(dist.list = dist.list, minmax.list = minmax),0)
  expect_equal(check.dists(dist.list = dist.list, minmax.list = quant),length(dist.list))
})

test_that("coords2indices returns the correct indices", {
  all.nonNA <- get_nonNA(use.distmat = T, rast = madagascar)[["all.nonNA"]]
  coords <- raster::xyFromCell(madagascar, raster::Which(!is.na(madagascar),cells=T)[1:3])
  res <- coords2indices(all.nonNA = all.nonNA, rast = madagascar, coords = coords)
  expect_true(all(res == c(1,2,3)))
})

test_that("indices2coords returns the correct coords", {
  all.nonNArev <- get_nonNA(use.distmat = T, rast = madagascar)[["all.nonNA.rev"]]
  exp <- raster::xyFromCell(madagascar, raster::Which(!is.na(madagascar),cells=T)[1:3])
  res <- indices2coords(all.nonNArev = all.nonNArev, rast = madagascar, indices = 1:3, species.vec = rep("a",3))
  expect_equal(res[,"x"], exp[,"x"])
  expect_equal(res[,"y"], exp[,"y"])
})

test_that("coords2pts produces correct output given use.distmat", {
  all.nonNA <- get_nonNA(use.distmat = T, rast = madagascar)[["all.nonNA"]]
  udmF <- coords2pts(use.distmat = F, coords = chameleons)
  udmT <- coords2pts(use.distmat = T, all.nonNA = all.nonNA, rast = madagascar, coords = chameleons)
  expect_true(is.matrix(udmF))
  expect_type(udmT, "integer")
})

test_that("pts2coords produces correct output", {
  all.nonNA <- get_nonNA(use.distmat = T, rast = madagascar)
  pts.udmF <- coords2pts(use.distmat = F, coords = chameleons)
  pts.udmT <- coords2pts(use.distmat = T, all.nonNA = all.nonNA[["all.nonNA"]], rast = madagascar, coords = chameleons)
  res.udmF <- pts2coords(use.distmat = F, pts = pts.udmF, pts.species = chameleons$species)
  res.udmT <- pts2coords(use.distmat = F, pts = pts.udmF, pts.species = chameleons$species, all.nonNArev = all.nonNA[["all.nonNA.rev"]], rast = madagascar)
  expect_equal(res.udmF$x, chameleons$x)
  expect_equal(res.udmF$y, chameleons$y)
  expect_equal(res.udmF$species, chameleons$species)
  expect_equal(res.udmT$x, chameleons$x)
  expect_equal(res.udmT$y, chameleons$y)
  expect_equal(res.udmT$species, chameleons$species)
})

test_that("KL.div returns 0 for identical distributions", {
  expect_equal(KL.div(1:100,1:100,20)[[1]],0)
})

test_that("KL.div returns higher values for more dissimilar distributions", {
  d1 <- c(1,2:51,112)
  d2 <- c(1,22:71,112)
  d3 <- c(1,42:91,112)
  d4 <- c(1,62:111,112)
  expect_true(all(rank(c(KL.div(d1,d2,20)[[1]], KL.div(d1,d3,20)[[1]], KL.div(d1,d4,20)[[1]])) == 1:3))
  expect_true(all(rank(c(KL.div(d4,d3,20)[[1]], KL.div(d4,d2,20)[[1]], KL.div(d4,d1,20)[[1]])) == 1:3))
})

test_that("KL.div.all produces expected output format", {
  combs <- get.combs(sep.inter.spp = T, sort(unique(chameleons$species)))
  dist.list.1 <- get_dists(use.distmat = F, distmat = NULL, pts = as.matrix(chameleons[,1:2]), indices = NULL, species.vec = chameleons$species, inter.spp = T, sep.inter.spp = T, dist_meth = "distRcpp", dist_fun = "Haversine", combs = combs)
  dist.list.2 <- get_dists(use.distmat = F, distmat = NULL, pts = as.matrix(chameleons[,1:2]), indices = NULL, species.vec = chameleons$species, inter.spp = T, sep.inter.spp = T, dist_meth = "distRcpp", dist_fun = "Haversine", combs = combs)
  exp <- rep(0,10)
  names(exp) <- c(paste("intra",sort(unique(chameleons$species)),sep = "_"),colnames(combs))
  res <- KL.div.all(dist.list.1 = dist.list.1, dist.list.2 = dist.list.2, break.num = 20)
  expect_identical(res, exp)
})
