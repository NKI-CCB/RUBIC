library(R.matlab, quietly=T, warn.conflict=F)

context("Parameters estimation - Case 0")

test_that("cna.shuffle has the same output as the model", {
  mat.in <- readMat('../referenceio/0/shuffleCNA.genPermProfiles.1.in.mat')
  mat.out <- readMat('../referenceio/0/shuffleCNA.genPermProfiles.1.out.mat')
  
  r.out <- shuffle.cna(mat.in$cna, 
                       as.vector(mat.out$shiftI))
  
  expect_equal(r.out, mat.out$cna)
})

test_that("gen.perm.profiles has the same output as the model", {
  mat.in <- readMat('../referenceio/0/genPermProfiles.in.mat')
  mat.out <- readMat('../referenceio/0/genPermProfiles.out.mat')
  
  rand.mat <- read.table('../referenceio/0/randi.genperm.csv')
  test.env <- new.env()
  test.env$random.matrix <- data.matrix(rand.mat)
  
  r.out <- gen.perm.profiles(mat.in$cna, 
                             as.numeric(mat.in$numIter),
                             test.env=test.env)

  expect_equal(r.out$agr.profiles, mat.out$agrProfiles)
  expect_equal(r.out$agr.cum.profiles, mat.out$agrCumProfiles)
})

test_that("calc.single.scale.roughness has the same output as the model", {
  mat.in <- readMat('../referenceio/0/calcSingleScaleRoughness.in.mat')
  mat.out <- readMat('../referenceio/0/calcSingleScaleRoughness.out.mat')
  
  r.out <- calc.single.scale.roughness(mat.in$lrCum, 
                                       as.numeric(mat.in$kw),
                                       as.numeric(mat.in$kwBreak),
                                       as.numeric(mat.in$P))
  
  expect_equal(matrix(list(r.out),nrow=1,ncol=1), 
               as.r.per.sample(list(list(mat.out$perSamp)), 1, 1))
})

test_that("single.samp.param.estimation has the same output as the model", {
  mat.in <- readMat('../referenceio/0/getSingleSampParamEsts.in.mat')
  mat.out <- readMat('../referenceio/0/getSingleSampParamEsts.out.mat')

  # cum.perm.profiles is taken from the output struct, even though it's an input
  # it can be done because it is only generated once from the inputs and never changed    
  r.out <- single.samp.param.estimation(mat.out$params[,,1]$agrCumIter, 
                                        as.numeric(mat.out$maxKw), 
                                        as.numeric(mat.in$data[,,1]$P))
  
  nr <- nrow(r.out$per.sample)
  nc <- ncol(r.out$per.sample)
  
  expect_equal(r.out$kws, as.vector(mat.out$kws))
  expect_equal(r.out$per.sample, 
               as.r.per.sample(mat.out$perSamp, nr, nc))
  
})

test_that("single.samp.param.estimation.abs has the same output as the model", {
  mat.in <- readMat('../referenceio/0/getSingleSampParamEstsAbs.in.mat')
  mat.out <- readMat('../referenceio/0/getSingleSampParamEstsAbs.out.mat')
  
  mat.maploc <- read.rubic.maploc(mat.in)
  
  r.out <- single.samp.param.estimation.abs(max.chrom.length(mat.maploc),
                                            ncol(mat.in$data[,,1]$cna),
                                            mat.in$params[,,1]$agrPCumIter,
                                            mat.in$params[,,1]$agrNCumIter,
                                            as.numeric(mat.in$tSign))
  
  nr <- nrow(r.out$per.sample)
  nc <- ncol(r.out$per.sample)
  
  expect_equal(r.out$kws, as.vector(mat.out$kws))
  expect_equal(r.out$per.sample, 
               as.r.per.sample(mat.out$perSamp, nr, nc))
      
})

test_that("all.samp.param.estimation has the same output as the model", {
  mat.in <- readMat('../referenceio/0/getAllSampParamEsts.in.mat')
  mat.out <- readMat('../referenceio/0/getAllSampParamEsts.out.mat')
  
  mat.in.kws <- as.vector(mat.in$params[,,1]$kws)
  mat.in.numProbes <- as.numeric(mat.in$data[,,1]$P)  
  mat.in.per.sample <- mat.in$params[,,1]$perSamp
  
  nr <- length(mat.in.kws)
  nc <- length(mat.in.per.sample)/nr
  r.out <- all.samp.param.estimation(mat.in.kws, 
                                     as.r.per.sample(mat.in.per.sample, nr, nc), 
                                     mat.in.numProbes)
  
  expect_equal(r.out$mus, mat.out$params[,,1]$mus)
  expect_equal(r.out$vars, mat.out$params[,,1]$vars)
  expect_equal(r.out$int.sqrt.vs, mat.out$params[,,1]$intSqrtVs)
  
})

test_that("cancel.recurrent.breaks has the same output as the model", {
  mat.in <- readMat('../referenceio/0/cancelRecurrentBreaks.in.mat')
  mat.out <- readMat('../referenceio/0/cancelRecurrentBreaks.out.mat') 
  
  r.out <- cancel.recurrent.breaks(mat.in$cna, 
                                   lapply(mat.in$modules, as.r.module.long))
  
  expect_equal(r.out, mat.out$cna)
})

test_that("iter.extract.single.samp.params has the same output as the model", {
  mat.in <- readMat('../referenceio/0/iterExtractSingleSampParams.in.mat')
  mat.out <- readMat('../referenceio/0/iterExtractSingleSampParams.out.mat')
  
  mat.maploc <- read.rubic.maploc(mat.in)
  mat.map.loc.agr <- sum.map.loc(mat.maploc)
  mat.max.chrom.len <- max.chrom.length(mat.maploc)
  
  rand.mat <- read.table('../referenceio/0/randi.all.csv')
  test.env <- new.env()
  test.env$random.matrix <- data.matrix(rand.mat)
  
  r.out <- iter.extract.single.samp.params(mat.in$data[,,1]$cna,
                                           mat.map.loc.agr,
                                           mat.max.chrom.len,
                                           as.numeric(mat.in$pValue),
                                           as.vector(mat.in$sampsUse),
                                           test.env=test.env)
  
  expect_equal(r.out$params, list(kws=as.vector(mat.out$params[,,1]$kws),
                                  mus=mat.out$params[,,1]$mus,
                                  vars=mat.out$params[,,1]$vars,
                                  int.sqrt.vs=mat.out$params[,,1]$intSqrtVs))
  # The iteration of the function cancel.recurrent.breaks produce a numerical error in the order of 1.14e-15
  expect_true(all.equal(r.out$canc.cna.matrix, mat.out$dataCanc[,,1]$cna, check.attributes=FALSE, tolerance=3e-15))
  
  mat.modules <- lapply(mat.out$modules, as.r.module.long)
  for (i in 1:length(mat.modules)) {
    mat.modules[[i]]$p <- NULL
  }
  
  for (i in 1:length(r.out$modules)) {
    r.out$modules[[i]] <- r.out$modules[[i]][c("I","kw","probes.ignore","l","r","amp","p")]
    r.out$modules[[i]]$p <- NULL
  }
  
  expect_equal(r.out$modules, mat.modules)
})

test_that("cna.matrix.to.segs has the same output as the model", {
  mat.in <- readMat('../referenceio/0/cna2segs.in.mat')
  mat.out <- readMat('../referenceio/0/cna2segs.out.mat')  
  
  r.out <- cna.matrix.to.segs(mat.in$cna,
                              lapply(mat.in$segments, as.r.module.long))

  expect_equal(r.out, mat.out$segMatrix)  
})

test_that("compute.gauss.kernel has the same output as the model", {
  mat.in <- readMat('../referenceio/0/getGaussKernel.in.mat')
  mat.out <- readMat('../referenceio/0/getGaussKernel.out.mat')  
  
  r.out <- compute.gauss.kernel(as.numeric(mat.in$numStds),
                                as.numeric(mat.in$kw))

  expect_equal(r.out, as.vector(mat.out$gaussKernel))
})

test_that("max.mode has the same output as the model", {
  mat.in <- readMat('../referenceio/0/getMaxMode.in.mat')
  mat.out <- readMat('../referenceio/0/getMaxMode.out.mat')  
  
  r.out <- max.mode(as.vector(mat.in$profile), 
                    as.numeric(mat.in$filterStd))

  expect_equal(r.out$x, as.vector(mat.out$x))
  expect_equal(r.out$y.orig, as.vector(mat.out$yOrig))
  expect_equal(r.out$z, as.vector(mat.out$z))
  expect_equal(r.out$peak, as.vector(mat.out$peak))
})

test_that("gen.perm.abs.profiles has the same output as the model", {
  mat.in <- readMat('../referenceio/0/genPermAbsProfiles.in.mat')
  mat.out <- readMat('../referenceio/0/genPermAbsProfiles.out.mat')
  
  rand.mat <- read.table('../referenceio/0/randi.all.csv')
  test.env <- new.env()
  test.env$random.matrix <- data.matrix(rand.mat)
  
  r.out <- gen.perm.abs.profiles(mat.in$cna,
                                 as.vector(mat.in$dataAgr),
                                 as.numeric(mat.in$numIter),
                                 as.numeric(mat.in$ampLevel),
                                 as.numeric(mat.in$ampSign),
                                 test.env=test.env)
  
  expect_equal(r.out$agr.profiles, mat.out$agrProfiles)
  expect_equal(r.out$agr.cum.profiles, mat.out$agrCumProfiles)
})

test_that("call.sig.abs.segments has the same output as the model", {
  mat.in <- readMat('../referenceio/0/callSigAbsSegments.in.mat')
  mat.out <- readMat('../referenceio/0/callSigAbsSegments.out.mat')
  
  r.out <- call.sig.abs.segments(as.vector(mat.in$cnaAbsAgr),
                                 mat.in$agrCumPerms,
                                 lapply(mat.in$segments, as.r.module.long),
                                 as.numeric(mat.in$pValue),
                                 as.numeric(mat.in$ampSign))
  
  expect_equivalent(r.out, as.logical(mat.out$sigIs))
})

test_that("null.seg.matrix has the same output as the model", {
  mat.in <- readMat('../referenceio/0/getNullSegMatrix.in.mat')
  mat.out <- readMat('../referenceio/0/getNullSegMatrix.out.mat')
  
  r.out <- null.seg.matrix(mat.in$segMatrix,
                           lapply(mat.in$segments, as.r.module.long),
                           as.vector(mat.in$sigI))
  
  expect_equivalent(r.out, mat.out$segMatrixNull)
})

test_that("segs.to.cna.matrix has the same output as the model", {
  mat.in <- readMat('../referenceio/0/segs2cna.in.mat')
  mat.out <- readMat('../referenceio/0/segs2cna.out.mat') 
  
  r.out <- segs.to.cna.matrix(mat.in$segMatrix,
                              lapply(mat.in$segments, as.r.module.long),
                              as.numeric(mat.in$P))
  
  expect_equal(r.out, mat.out$cna)
})

test_that("update.abs.cna has the same output as the model", {
  mat.in <- readMat('../referenceio/0/updateAbsCNA.in.mat')
  mat.out <- readMat('../referenceio/0/updateAbsCNA.out.mat')
  
  rand.mat <- read.table('../referenceio/0/randi.all.csv')
  test.env <- new.env()
  test.env$random.matrix <- data.matrix(rand.mat)
  
  r.out <- update.abs.cna(mat.in$data[,,1]$cna,
                          as.numeric(mat.in$data[,,1]$ampLevel),
                          as.numeric(mat.in$data[,,1]$delLevel),
                          mat.in$segMatrix,
                          mat.in$cna,
                          as.vector(mat.in$cnaAgrP),
                          as.vector(mat.in$cnaAgrN),
                          lapply(mat.in$segments, as.r.module.long),
                          as.numeric(mat.in$pValue),
                          as.numeric(mat.in$data[,,1]$P),
                          test.env=test.env)

  expect_equivalent(r.out$cna.matrix, mat.out$cna)
  expect_equivalent(r.out$num.sig.segs, as.numeric(mat.out$numSigSegs))
})

test_that("iter.update.abs.cna has the same output as the model", {
  mat.in <- readMat('../referenceio/0/iterUpdateAbsCNA.in.mat')
  mat.out <- readMat('../referenceio/0/iterUpdateAbsCNA.out.mat') 
  
  mat.maploc <- read.rubic.maploc(mat.in)
  
  rand.mat <- read.table('../referenceio/0/randi.all.csv')
  test.env <- new.env()
  test.env$random.matrix <- data.matrix(rand.mat)
  
  r.out <- iter.update.abs.cna(mat.in$data[,,1]$cna,
                               as.numeric(mat.in$data[,,1]$ampLevel),
                               as.numeric(mat.in$data[,,1]$delLevel),
                               lapply(mat.in$segments, as.r.module.long),
                               as.numeric(mat.in$pValue),
                               test.env=test.env)
  
  expect_equal(r.out$cna.matrix, mat.out$cna)
  expect_equal(r.out$agr.p.perms, mat.out$agrPPerms)
  expect_equal(r.out$agr.n.perms, mat.out$agrNPerms)
  expect_equal(r.out$agr.p.cum.perms, mat.out$agrPCumPerms)
  expect_equal(r.out$agr.n.cum.perms, mat.out$agrNCumPerms)
})

test_that("estimate.parameters has the same output as the model", {
  mat.in <- readMat('../referenceio/0/estimateParametersBlock.in.mat')
  mat.out <- readMat('../referenceio/0/estimateParametersBlock.out.mat') 
  
  mat.maploc <- read.rubic.maploc(mat.in)
  
  rand.mat <- read.table('../referenceio/0/randi.all.csv')
  test.env <- new.env()
  test.env$random.matrix <- data.matrix(rbind(rand.mat,rand.mat))
  
  mat.cna.matrix <- extract.matrix(mat.maploc)
  mat.map.loc.agr <- sum.map.loc(mat.maploc)
  mat.max.chrom.len <- max.chrom.length(mat.maploc)
  
  r.out <- estimate.parameters(mat.cna.matrix, mat.map.loc.agr,
                               mat.max.chrom.len,
                               as.numeric(mat.in$data[,,1]$ampLevel),
                               as.numeric(mat.in$data[,,1]$delLevel),
                               as.numeric(mat.in$FDR),
                               test.env=test.env)
  
  
  expect_equal(r.out$params.p$kws, as.numeric(mat.out$paramsP[,,1]$kws))
  expect_equal(r.out$params.p$int.sqrt.vs, mat.out$paramsP[,,1]$intSqrtVs)
  expect_equal(r.out$params.p$mus, mat.out$paramsP[,,1]$mus)
  expect_equal(r.out$params.p$vars, mat.out$paramsP[,,1]$vars)
  expect_equal(r.out$params.p$agr.cum.iter, mat.out$paramsP[,,1]$agrCumIter)
  expect_equal(r.out$params.p$per.sample, as.r.per.sample(mat.out$paramsP[,,1]$perSamp,
                                                          nrow(r.out$params.p$per.sample),
                                                          ncol(r.out$params.p$per.sample)))
  
  expect_equal(r.out$params.n$kws, as.numeric(mat.out$paramsN[,,1]$kws))
  expect_equal(r.out$params.n$int.sqrt.vs, mat.out$paramsN[,,1]$intSqrtVs)
  expect_equal(r.out$params.n$mus, mat.out$paramsN[,,1]$mus)
  expect_equal(r.out$params.n$vars, mat.out$paramsN[,,1]$vars)
  expect_equal(r.out$params.n$agr.cum.iter, mat.out$paramsN[,,1]$agrCumIter)
  expect_equal(r.out$params.n$per.sample, as.r.per.sample(mat.out$paramsN[,,1]$perSamp,
                                                          nrow(r.out$params.n$per.sample),
                                                          ncol(r.out$params.n$per.sample)))
})

context("Parameters estimation - Case 1")

test_that("cna.shuffle has the same output as the model", {
  mat.in <- readMat('../referenceio/1/genPermProfiles.in.mat')
  mat.out <- readMat('../referenceio/1/genPermProfiles.out.mat')

  rand.mat <- read.table('../referenceio/1/randi.all.csv')
  rand.vect <- data.matrix(rand.mat)[10,]
  r.out <- shuffle.cna(mat.in$cna, 
                       random.vector=rand.vect)
  
  expect_equal(r.out, mat.out$cnaShuff)
})

test_that("gen.perm.profiles has the same output as the model", {
  mat.in <- readMat('../referenceio/1/genPermProfiles.in.mat')
  mat.out <- readMat('../referenceio/1/genPermProfiles.out.mat')
  
  rand.mat <- read.table('../referenceio/1/randi.all.csv')
  test.env <- new.env()
  test.env$random.matrix <- data.matrix(rand.mat)
  
  r.out <- gen.perm.profiles(mat.in$cna, 
                             as.numeric(mat.in$numIter),
                             test.env=test.env)
  
  expect_equal(r.out$agr.profiles, mat.out$agrProfiles)
  expect_equal(r.out$agr.cum.profiles, mat.out$agrCumProfiles)
})

test_that("calc.single.scale.roughness has the same output as the model", {
  mat.in <- readMat('../referenceio/1/calcSingleScaleRoughness.in.mat')
  mat.out <- readMat('../referenceio/1/calcSingleScaleRoughness.out.mat')
  
  r.out <- calc.single.scale.roughness(mat.in$lrCum, 
                                       as.numeric(mat.in$kw),
                                       as.numeric(mat.in$kwBreak),
                                       as.numeric(mat.in$P))
  
  expect_equal(matrix(list(r.out),nrow=1,ncol=1), 
               as.r.per.sample(list(list(mat.out$perSamp)), 1, 1))
})

test_that("single.samp.param.estimation has the same output as the model", {
  mat.in <- readMat('../referenceio/1/getSingleSampParamEsts.in.mat')
  mat.out <- readMat('../referenceio/1/getSingleSampParamEsts.out.mat')
  
  # cum.perm.profiles is taken from the output struct, even though it's an input
  # it can be done because it is only generated once from the inputs and never changed    
  r.out <- single.samp.param.estimation(mat.out$params[,,1]$agrCumIter, 
                                        as.numeric(mat.out$maxKw), 
                                        as.numeric(mat.in$data[,,1]$P))
  
  nr <- nrow(r.out$per.sample)
  nc <- ncol(r.out$per.sample)
  
  expect_equal(r.out$kws, as.vector(mat.out$kws))
  expect_equal(r.out$per.sample, 
               as.r.per.sample(mat.out$perSamp, nr, nc))
  
})

test_that("single.samp.param.estimation.abs has the same output as the model", {
  mat.in <- readMat('../referenceio/1/getSingleSampParamEstsAbs.in.mat')
  mat.out <- readMat('../referenceio/1/getSingleSampParamEstsAbs.out.mat')
  
  mat.maploc <- read.rubic.maploc(mat.in)
  
  r.out <- single.samp.param.estimation.abs(max.chrom.length(mat.maploc),
                                            ncol(mat.in$data[,,1]$cna),
                                            mat.in$params[,,1]$agrPCumIter,
                                            mat.in$params[,,1]$agrNCumIter,
                                            as.numeric(mat.in$tSign))
  
  nr <- nrow(r.out$per.sample)
  nc <- ncol(r.out$per.sample)
  
  expect_equal(r.out$kws, as.vector(mat.out$kws))
  expect_equal(r.out$per.sample, 
               as.r.per.sample(mat.out$perSamp, nr, nc))
  
})

test_that("all.samp.param.estimation has the same output as the model", {
  mat.in <- readMat('../referenceio/1/getAllSampParamEsts.in.mat')
  mat.out <- readMat('../referenceio/1/getAllSampParamEsts.out.mat')
  
  mat.in.kws <- as.vector(mat.in$params[,,1]$kws)
  mat.in.numProbes <- as.numeric(mat.in$data[,,1]$P)  
  mat.in.per.sample <- mat.in$params[,,1]$perSamp
  
  nr <- length(mat.in.kws)
  nc <- length(mat.in.per.sample)/nr
  r.out <- all.samp.param.estimation(mat.in.kws, 
                                     as.r.per.sample(mat.in.per.sample, nr, nc), 
                                     mat.in.numProbes)
  
  expect_equal(r.out$mus, mat.out$params[,,1]$mus)
  expect_equal(r.out$vars, mat.out$params[,,1]$vars)
  expect_equal(r.out$int.sqrt.vs, mat.out$params[,,1]$intSqrtVs)
  
})

test_that("cancel.recurrent.breaks has the same output as the model", {
  mat.in <- readMat('../referenceio/1/cancelRecurrentBreaks.in.mat')
  mat.out <- readMat('../referenceio/1/cancelRecurrentBreaks.out.mat') 
  
  r.out <- cancel.recurrent.breaks(mat.in$cna, 
                                   lapply(mat.in$modules, as.r.module.long))
  
  expect_equal(r.out, mat.out$cna)
})

test_that("iter.extract.single.samp.params has the same output as the model", {
  mat.in <- readMat('../referenceio/1/iterExtractSingleSampParams.in.mat')
  mat.out <- readMat('../referenceio/1/iterExtractSingleSampParams.out.mat')
  
  mat.maploc <- read.rubic.maploc(mat.in)
  mat.map.loc.agr <- sum.map.loc(mat.maploc)
  mat.max.chrom.len <- max.chrom.length(mat.maploc)
  
  rand.mat <- read.table('../referenceio/1/randi.all.csv')
  test.env <- new.env()
  test.env$random.matrix <- data.matrix(rand.mat)
  
  r.out <- iter.extract.single.samp.params(mat.in$data[,,1]$cna,
                                           mat.map.loc.agr,
                                           mat.max.chrom.len,
                                           as.numeric(mat.in$pValue),
                                           as.vector(mat.in$sampsUse),
                                           test.env=test.env)
  
  expect_equal(r.out$params, list(kws=as.vector(mat.out$params[,,1]$kws),
                                  mus=mat.out$params[,,1]$mus,
                                  vars=mat.out$params[,,1]$vars,
                                  int.sqrt.vs=mat.out$params[,,1]$intSqrtVs))
  # The iteration of the function cancel.recurrent.breaks produce a numerical error in the order of 1.14e-15
  expect_true(all.equal(r.out$canc.cna.matrix, mat.out$dataCanc[,,1]$cna, check.attributes=FALSE, tolerance=3e-15))
  
  mat.modules <- lapply(mat.out$modules, as.r.module.long)
  for (i in 1:length(mat.modules)) {
    mat.modules[[i]]$p <- NULL
  }
  
  for (i in 1:length(r.out$modules)) {
    r.out$modules[[i]] <- r.out$modules[[i]][c("I","kw","probes.ignore","l","r","amp","p")]
    r.out$modules[[i]]$p <- NULL
  }
  
  expect_equal(r.out$modules, mat.modules)
})

test_that("cna.matrix.to.segs has the same output as the model", {
  mat.in <- readMat('../referenceio/1/cna2segs.in.mat')
  mat.out <- readMat('../referenceio/1/cna2segs.out.mat')  
  
  r.out <- cna.matrix.to.segs(mat.in$cna,
                              lapply(mat.in$segments, as.r.module.long))
  
  expect_equal(r.out, mat.out$segMatrix)  
})



test_that("compute.gauss.kernel has the same output as the model", {
  mat.in <- readMat('../referenceio/1/getGaussKernel.in.mat')
  mat.out <- readMat('../referenceio/1/getGaussKernel.out.mat')  
  
  r.out <- compute.gauss.kernel(as.numeric(mat.in$numStds),
                                as.numeric(mat.in$kw))
  
  expect_equal(r.out, as.vector(mat.out$gaussKernel))
})

test_that("max.mode has the same output as the model", {
  mat.in <- readMat('../referenceio/1/getMaxMode.in.mat')
  mat.out <- readMat('../referenceio/1/getMaxMode.out.mat')  
  
  r.out <- max.mode(as.vector(mat.in$profile), 
                    as.numeric(mat.in$filterStd))
  
  expect_equal(r.out$x, as.vector(mat.out$x))
  expect_equal(r.out$y.orig, as.vector(mat.out$yOrig))
  expect_equal(r.out$z, as.vector(mat.out$z))
  expect_equal(r.out$peak, as.vector(mat.out$peak))
})

test_that("gen.perm.abs.profiles has the same output as the model", {
  mat.in <- readMat('../referenceio/1/genPermAbsProfiles.in.mat')
  mat.out <- readMat('../referenceio/1/genPermAbsProfiles.out.mat')
  
  rand.mat <- read.table('../referenceio/1/randi.all.csv')
  test.env <- new.env()
  test.env$random.matrix <- data.matrix(rand.mat)
  
  r.out <- gen.perm.abs.profiles(mat.in$cna,
                                 as.vector(mat.in$dataAgr),
                                 as.numeric(mat.in$numIter),
                                 as.numeric(mat.in$ampLevel),
                                 as.numeric(mat.in$ampSign),
                                 test.env=test.env)
  
  expect_equal(r.out$agr.profiles, mat.out$agrProfiles)
  expect_equal(r.out$agr.cum.profiles, mat.out$agrCumProfiles)
})

test_that("call.sig.abs.segments has the same output as the model", {
  mat.in <- readMat('../referenceio/1/callSigAbsSegments.in.mat')
  mat.out <- readMat('../referenceio/1/callSigAbsSegments.out.mat')
  
  r.out <- call.sig.abs.segments(as.vector(mat.in$cnaAbsAgr),
                                 mat.in$agrCumPerms,
                                 lapply(mat.in$segments, as.r.module.long),
                                 as.numeric(mat.in$pValue),
                                 as.numeric(mat.in$ampSign))
  
  expect_equivalent(r.out, as.logical(mat.out$sigIs))
})

test_that("null.seg.matrix has the same output as the model", {
  mat.in <- readMat('../referenceio/1/getNullSegMatrix.in.mat')
  mat.out <- readMat('../referenceio/1/getNullSegMatrix.out.mat')
  
  r.out <- null.seg.matrix(mat.in$segMatrix,
                           lapply(mat.in$segments, as.r.module.long),
                           as.vector(mat.in$sigI))
  
  expect_equivalent(r.out, mat.out$segMatrixNull)
})

test_that("segs.to.cna.matrix has the same output as the model", {
  mat.in <- readMat('../referenceio/1/segs2cna.in.mat')
  mat.out <- readMat('../referenceio/1/segs2cna.out.mat') 
  
  r.out <- segs.to.cna.matrix(mat.in$segMatrix,
                              lapply(mat.in$segments, as.r.module.long),
                              as.numeric(mat.in$P))
  
  expect_equal(r.out, mat.out$cna)
})

test_that("update.abs.cna has the same output as the model", {
  mat.in <- readMat('../referenceio/1/updateAbsCNA.in.mat')
  mat.out <- readMat('../referenceio/1/updateAbsCNA.out.mat')
  
  rand.mat <- read.table('../referenceio/1/randi.all.csv')
  test.env <- new.env()
  test.env$random.matrix <- data.matrix(rand.mat)
  
  r.out <- update.abs.cna(mat.in$data[,,1]$cna,
                          as.numeric(mat.in$data[,,1]$ampLevel),
                          as.numeric(mat.in$data[,,1]$delLevel),
                          mat.in$segMatrix,
                          mat.in$cna,
                          as.vector(mat.in$cnaAgrP),
                          as.vector(mat.in$cnaAgrN),
                          lapply(mat.in$segments, as.r.module.long),
                          as.numeric(mat.in$pValue),
                          as.numeric(mat.in$data[,,1]$P),
                          test.env=test.env)
  
  expect_equivalent(r.out$cna.matrix, mat.out$cna)
  expect_equivalent(r.out$num.sig.segs, as.numeric(mat.out$numSigSegs))
})

test_that("iter.update.abs.cna has the same output as the model", {
  mat.in <- readMat('../referenceio/1/iterUpdateAbsCNA.in.mat')
  mat.out <- readMat('../referenceio/1/iterUpdateAbsCNA.out.mat') 
  
  mat.maploc <- read.rubic.maploc(mat.in)
  
  rand.mat <- read.table('../referenceio/1/randi.all.csv')
  test.env <- new.env()
  test.env$random.matrix <- data.matrix(rand.mat)
  
  r.out <- iter.update.abs.cna(mat.in$data[,,1]$cna,
                               as.numeric(mat.in$data[,,1]$ampLevel),
                               as.numeric(mat.in$data[,,1]$delLevel),
                               lapply(mat.in$segments, as.r.module.long),
                               as.numeric(mat.in$pValue),
                               test.env=test.env)
  
  expect_equal(r.out$cna.matrix, mat.out$cna)
  expect_equal(r.out$agr.p.perms, mat.out$agrPPerms)
  expect_equal(r.out$agr.n.perms, mat.out$agrNPerms)
  expect_equal(r.out$agr.p.cum.perms, mat.out$agrPCumPerms)
  expect_equal(r.out$agr.n.cum.perms, mat.out$agrNCumPerms)
})

test_that("estimate.parameters has the same output as the model", {
  mat.in <- readMat('../referenceio/1/estimateParametersBlock.in.mat')
  mat.out <- readMat('../referenceio/1/estimateParametersBlock.out.mat') 
  
  mat.maploc <- read.rubic.maploc(mat.in)
  
  rand.mat <- read.table('../referenceio/1/randi.all.csv')
  test.env <- new.env()
  test.env$random.matrix <- data.matrix(rbind(rand.mat,rand.mat))
  
  mat.cna.matrix <- extract.matrix(mat.maploc)
  mat.map.loc.agr <- sum.map.loc(mat.maploc)
  mat.max.chrom.len <- max.chrom.length(mat.maploc)
  
  r.out <- estimate.parameters(mat.cna.matrix, mat.map.loc.agr,
                               mat.max.chrom.len,
                               as.numeric(mat.in$data[,,1]$ampLevel),
                               as.numeric(mat.in$data[,,1]$delLevel),
                               as.numeric(mat.in$FDR),
                               test.env=test.env)
  
  
  expect_equal(r.out$params.p$kws, as.numeric(mat.out$paramsP[,,1]$kws))
  expect_equal(r.out$params.p$int.sqrt.vs, mat.out$paramsP[,,1]$intSqrtVs)
  expect_equal(r.out$params.p$mus, mat.out$paramsP[,,1]$mus)
  expect_equal(r.out$params.p$vars, mat.out$paramsP[,,1]$vars)
  expect_equal(r.out$params.p$agr.cum.iter, mat.out$paramsP[,,1]$agrCumIter)
  expect_equal(r.out$params.p$per.sample, as.r.per.sample(mat.out$paramsP[,,1]$perSamp,
                                                          nrow(r.out$params.p$per.sample),
                                                          ncol(r.out$params.p$per.sample)))
  
  expect_equal(r.out$params.n$kws, as.numeric(mat.out$paramsN[,,1]$kws))
  expect_equal(r.out$params.n$int.sqrt.vs, mat.out$paramsN[,,1]$intSqrtVs)
  expect_equal(r.out$params.n$mus, mat.out$paramsN[,,1]$mus)
  expect_equal(r.out$params.n$vars, mat.out$paramsN[,,1]$vars)
  expect_equal(r.out$params.n$agr.cum.iter, mat.out$paramsN[,,1]$agrCumIter)
  expect_equal(r.out$params.n$per.sample, as.r.per.sample(mat.out$paramsN[,,1]$perSamp,
                                                          nrow(r.out$params.n$per.sample),
                                                          ncol(r.out$params.n$per.sample)))
})