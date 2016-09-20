library(R.matlab, quietly=T, warn.conflict=F)

# Autotesters for Case 0
# 40 samples tested
context("Call events - Case 0")

test_that("min.kw.peak.index has the same output as the model", {
  mat.in <- readMat('../referenceio/0/getMinKwPeak.in.mat')
  mat.out <- readMat('../referenceio/0/getMinKwPeak.out.mat')
  
  r.out <- min.kw.peak.index(lapply(mat.in$segments, as.r.module.long), 
                             as.numeric(mat.in$tSign))
  
  expect_equal(r.out, as.numeric(mat.out$minKwPeakI))
})


test_that("percentile.v2 has the same output as the model", {
  mat.in <- readMat('../referenceio/0/getPercentileV2.in.mat')
  mat.out <- readMat('../referenceio/0/getPercentileV2.out.mat')
  
  r.out <- percentile.v2(as.numeric(mat.in$t),
                         as.numeric(mat.in$P),
                         mat.in$agrCumIter,
                         as.numeric(mat.in$kw),
                         as.numeric(mat.in$tSign))
  
  expect_equal(r.out$perc, as.numeric(mat.out$perc))
  expect_equal(r.out$z, as.numeric(mat.out$z))
})


test_that("is.peak.valid has the same output as the model", {
  mat.in <- readMat('../referenceio/0/isPeakValid.in.mat')
  mat.out <- readMat('../referenceio/0/isPeakValid.out.mat') 
  
  r.out <- is.peak.valid(mat.in$params[,,1]$P, 
                         mat.in$params[,,1]$agrCumIter,
                         lapply(mat.in$segments, as.r.module.long),
                         as.numeric(mat.in$minKwPeakI),
                         as.vector(mat.in$cnaAgr),
                         as.numeric(mat.in$pValue),
                         as.numeric(mat.in$tSign))
  
  expect_equal(r.out$peak.valid, as.logical(mat.out$peakValid))
  expect_equal(r.out$sig.info$amp, as.numeric(mat.out$sigInfo[,,1]$amp))
  expect_equal(r.out$sig.info$p, as.numeric(mat.out$sigInfo[,,1]$p))
})


test_that("broadest.adj has the same output as the model", {
  mat.in <- readMat('../referenceio/0/getBroadestAdj.in.mat')
  mat.out <- readMat('../referenceio/0/getBroadestAdj.out.mat') 
  
  r.out <- broadest.adj(mat.in$focalEvents, 
                            as.numeric(mat.in$segLoc),
                            as.numeric(mat.in$type))
  
  expect_equal(r.out, as.numeric(mat.out$I))
  
})


test_that("adjust.boundaries has the same output as the model", {
  mat.in <- readMat('../referenceio/0/adjustBoundaries.in.mat')
  mat.out <- readMat('../referenceio/0/adjustBoundaries.out.mat') 
  
  r.out <- adjust.boundaries(mat.in$focalEvents, 
                             as.r.module.long(list(mat.in$segment)))
  
  expect_equal(r.out, as.r.module.long(list(mat.out$segment)))
  
})


test_that("join.adj.segs has the same output as the model", {
  mat.in <- readMat('../referenceio/0/joinAdjSegs.in.mat')
  mat.out <- readMat('../referenceio/0/joinAdjSegs.out.mat') 
  
  r.out <- join.adj.segs(mat.in$cna,
                         lapply(mat.in$segments, as.r.module.long),
                         as.numeric(mat.in$relPeakI),
                         list(kws=as.vector(mat.in$params[,,1]$kws),
                              mus=mat.in$params[,,1]$mus,
                              vars=mat.in$params[,,1]$vars, 
                              int.sqrt.vs=mat.in$params[,,1]$intSqrtVs,
                              agr.cum.iter=mat.in$params[,,1]$agrCumIter),
                         as.numeric(mat.in$E))
  
  expect_equal(r.out, lapply(mat.out$segments, as.r.module.long))
})


test_that("call.events has the same output as the model", {
  mat.in <- readMat('../referenceio/0/callEvents.in.mat')
  mat.out <- readMat('../referenceio/0/callEvents.out.mat') 
  
  r.out <- call.events(mat.in$data[,,1]$cna,
                       as.numeric(mat.in$data[,,1]$ampLevel),
                       as.numeric(mat.in$data[,,1]$delLevel),
                       lapply(mat.in$segments, as.r.module.long),
                       list(kws=as.vector(mat.in$params[,,1]$kws),
                            mus=mat.in$params[,,1]$mus,
                            vars=mat.in$params[,,1]$vars,
                            int.sqrt.vs=mat.in$params[,,1]$intSqrtVs,
                            agr.cum.iter=mat.in$params[,,1]$agrCumIter),
                       as.numeric(mat.in$FDR),
                       as.numeric(mat.in$E),
                       as.numeric(mat.in$tSign),
                       samps.use=as.logical(mat.in$sampsUse))
  
  expect_equal(r.out, lapply(mat.out$focalEvents, as.r.module.long))
})


# Autotesters for Case 1
# 141 samples tested
context("Call events - Case 1")

test_that("min.kw.peak.index has the same output as the model", {
  mat.in <- readMat('../referenceio/1/getMinKwPeak.in.mat')
  mat.out <- readMat('../referenceio/1/getMinKwPeak.out.mat')
  
  r.out <- min.kw.peak.index(lapply(mat.in$segments, as.r.module.long), 
                             as.numeric(mat.in$tSign))
  
  expect_equal(r.out, as.numeric(mat.out$minKwPeakI))
})


test_that("percentile.v2 has the same output as the model", {
  mat.in <- readMat('../referenceio/1/getPercentileV2.in.mat')
  mat.out <- readMat('../referenceio/1/getPercentileV2.out.mat')
  
  r.out <- percentile.v2(as.numeric(mat.in$t),
                         as.numeric(mat.in$P),
                         mat.in$agrCumIter,
                         as.numeric(mat.in$kw),
                         as.numeric(mat.in$tSign))
  
  expect_equal(r.out$perc, as.numeric(mat.out$perc))
  expect_equal(r.out$z, as.numeric(mat.out$z))
})


test_that("is.peak.valid has the same output as the model", {
  mat.in <- readMat('../referenceio/1/isPeakValid.in.mat')
  mat.out <- readMat('../referenceio/1/isPeakValid.out.mat') 
  
  r.out <- is.peak.valid(mat.in$params[,,1]$P, 
                         mat.in$params[,,1]$agrCumIter,
                         lapply(mat.in$segments, as.r.module.long),
                         as.numeric(mat.in$minKwPeakI),
                         as.vector(mat.in$cnaAgr),
                         as.numeric(mat.in$pValue),
                         as.numeric(mat.in$tSign))
  
  expect_equal(r.out$peak.valid, as.logical(mat.out$peakValid))
  expect_equal(r.out$sig.info$amp, as.numeric(mat.out$sigInfo[,,1]$amp))
  expect_equal(r.out$sig.info$p, as.numeric(mat.out$sigInfo[,,1]$p))
})


test_that("broadest.adj has the same output as the model", {
  mat.in <- readMat('../referenceio/1/getBroadestAdj.in.mat')
  mat.out <- readMat('../referenceio/1/getBroadestAdj.out.mat') 
  
  r.out <- broadest.adj(mat.in$focalEvents, 
                        as.numeric(mat.in$segLoc),
                        as.numeric(mat.in$type))
  
  expect_equal(r.out, as.numeric(mat.out$I))
  
})


test_that("adjust.boundaries has the same output as the model", {
  mat.in <- readMat('../referenceio/1/adjustBoundaries.in.mat')
  mat.out <- readMat('../referenceio/1/adjustBoundaries.out.mat') 
  
  r.out <- adjust.boundaries(mat.in$focalEvents, 
                             as.r.module.long(list(mat.in$segment)))
  
  expect_equal(r.out, as.r.module.long(list(mat.out$segment)))
  
})


test_that("join.adj.segs has the same output as the model", {
  mat.in <- readMat('../referenceio/1/joinAdjSegs.in.mat')
  mat.out <- readMat('../referenceio/1/joinAdjSegs.out.mat') 
  
  r.out <- join.adj.segs(mat.in$cna,
                         lapply(mat.in$segments, as.r.module.long),
                         as.numeric(mat.in$relPeakI),
                         list(kws=as.vector(mat.in$params[,,1]$kws),
                              mus=mat.in$params[,,1]$mus,
                              vars=mat.in$params[,,1]$vars, 
                              int.sqrt.vs=mat.in$params[,,1]$intSqrtVs,
                              agr.cum.iter=mat.in$params[,,1]$agrCumIter),
                         as.numeric(mat.in$E))

  expect_equal(r.out, lapply(mat.out$segments, as.r.module.long))
})


test_that("call.events has the same output as the model", {
  mat.in <- readMat('../referenceio/1/callEvents.in.mat')
  mat.out <- readMat('../referenceio/1/callEvents.out.mat') 
  
  r.out <- call.events(mat.in$data[,,1]$cna,
                       as.numeric(mat.in$data[,,1]$ampLevel),
                       as.numeric(mat.in$data[,,1]$delLevel),
                       lapply(mat.in$segments, as.r.module.long),
                       list(kws=as.vector(mat.in$params[,,1]$kws),
                            mus=mat.in$params[,,1]$mus,
                            vars=mat.in$params[,,1]$vars,
                            int.sqrt.vs=mat.in$params[,,1]$intSqrtVs,
                            agr.cum.iter=mat.in$params[,,1]$agrCumIter),
                       as.numeric(mat.in$FDR),
                       as.numeric(mat.in$E),
                       as.numeric(mat.in$tSign),
                       samps.use=as.logical(mat.in$sampsUse))
  
  expect_equal(r.out, lapply(mat.out$focalEvents, as.r.module.long))
})