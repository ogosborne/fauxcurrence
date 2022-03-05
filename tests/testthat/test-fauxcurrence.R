library(raster)

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

test_that("sample.bear returns a single number between 0 and 360 if make.distmat is FALSE", {
  udmF <- sample.bear(use.distmat = F)
  expect_length(udmF, 1)
  expect_lte(udmF, 360)
  expect_gte(udmF, 0)
})

test_that("sample.bear returns NULL if make.distmat is TRUE", {
  udmT <- sample.bear(use.distmat = T)
  expect_null(udmT)
})
