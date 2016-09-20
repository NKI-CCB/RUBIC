library(R.matlab, quietly=T, warn.conflict=F)

context("Aggregate segmentation - Case 0")

test_that("segments.at.fdr.abs has the same output as the model", {
  mat.in <- readMat('../referenceio/0/getSegmentsAtFDRAbs.in.mat')
  mat.out <- readMat('../referenceio/0/getSegmentsAtFDRAbs.out.mat')
  
  mat.maploc <- read.rubic.maploc(mat.in)
  map.loc.agr <- sum.map.loc(mat.maploc)
  
  r.out <- segments.at.fdr.abs(mat.in$data[,,1]$cna,
                               map.loc.agr,
                               as.numeric(mat.in$data[,,1]$ampLevel),
                               as.numeric(mat.in$data[,,1]$delLevel),
                               as.vector(mat.in$sampsUse),
                               list(kws=as.vector(mat.in$params[,,1]$kws),
                                    mus=mat.in$params[,,1]$mus,
                                    vars=mat.in$params[,,1]$vars, 
                                    int.sqrt.vs=mat.in$params[,,1]$intSqrtVs,
                                    agr.cum.iter=mat.in$params[,,1]$agrCumIter),
                               as.numeric(mat.in$FDR),
                               as.numeric(mat.in$tSign))
  
  # reorder the individual modules by name according to the matlab convention
  for (i in 1:length(r.out$modules)) {
    r.out$modules[[i]] <- r.out$modules[[i]][c("I","kw","probes.ignore","l","r","amp","p")]
  }
  
  expect_equal(r.out$e, as.numeric(mat.out$E))
  expect_equal(r.out$modules, lapply(mat.out$modules, as.r.module.long))
  expect_equal(r.out$break.info, as.r.break.info.chrom(mat.out$breakInfo))
})

test_that("aggregate.segments has the same output as the model", {
  mat.in <- readMat('../referenceio/0/aggregateSegmentsBlock.in.mat')
  mat.out <- readMat('../referenceio/0/aggregateSegmentsBlock.out.mat')
  
  mat.maploc <- read.rubic.maploc(mat.in)
  map.loc.agr <- sum.map.loc(mat.maploc)
  
  r.out <- aggregate.segments(mat.in$data[,,1]$cna,
                              map.loc.agr,
                              as.numeric(mat.in$data[,,1]$ampLevel),
                              as.numeric(mat.in$data[,,1]$delLevel),
                              list(kws=as.vector(mat.in$paramsP[,,1]$kws),
                                   mus=mat.in$paramsP[,,1]$mus,
                                   vars=mat.in$paramsP[,,1]$vars, 
                                   int.sqrt.vs=mat.in$paramsP[,,1]$intSqrtVs,
                                   agr.cum.iter=mat.in$paramsP[,,1]$agrCumIter),
                              list(kws=as.vector(mat.in$paramsN[,,1]$kws),
                                   mus=mat.in$paramsN[,,1]$mus,
                                   vars=mat.in$paramsN[,,1]$vars, 
                                   int.sqrt.vs=mat.in$paramsN[,,1]$intSqrtVs,
                                   agr.cum.iter=mat.in$paramsN[,,1]$agrCumIter),
                              as.numeric(mat.in$FDR))
  
  # reorder the individual modules by name according to the matlab convention
  for (i in 1:length(r.out$segments.p)) {
    r.out$segments.p[[i]] <- r.out$segments.p[[i]][c("I","kw","probes.ignore","l","r","amp","p")]
  }
  for (i in 1:length(r.out$segments.n)) {
    r.out$segments.n[[i]] <- r.out$segments.n[[i]][c("I","kw","probes.ignore","l","r","amp","p")]
  }
  
  expect_equal(r.out$e.p, as.numeric(mat.out$EP))
  expect_equal(r.out$e.n, as.numeric(mat.out$EN))
  expect_equal(r.out$segments.p, lapply(mat.out$segmentsP, as.r.module.long))
  expect_equal(r.out$segments.n, lapply(mat.out$segmentsN, as.r.module.long))
})

##########
context("Aggregate segmentation - Case 1")

test_that("segments.at.fdr.abs has the same output as the model", {
  mat.in <- readMat('../referenceio/1/getSegmentsAtFDRAbs.in.mat')
  mat.out <- readMat('../referenceio/1/getSegmentsAtFDRAbs.out.mat')
  
  mat.maploc <- read.rubic.maploc(mat.in)
  map.loc.agr <- sum.map.loc(mat.maploc)
  
  r.out <- segments.at.fdr.abs(mat.in$data[,,1]$cna,
                               map.loc.agr,
                               as.numeric(mat.in$data[,,1]$ampLevel),
                               as.numeric(mat.in$data[,,1]$delLevel),
                               as.vector(mat.in$sampsUse),
                               list(kws=as.vector(mat.in$params[,,1]$kws),
                                    mus=mat.in$params[,,1]$mus,
                                    vars=mat.in$params[,,1]$vars, 
                                    int.sqrt.vs=mat.in$params[,,1]$intSqrtVs,
                                    agr.cum.iter=mat.in$params[,,1]$agrCumIter),
                               as.numeric(mat.in$FDR),
                               as.numeric(mat.in$tSign))
  
  # reorder the individual modules by name according to the matlab convention
  for (i in 1:length(r.out$modules)) {
    r.out$modules[[i]] <- r.out$modules[[i]][c("I","kw","probes.ignore","l","r","amp","p")]
  }
  
  expect_equal(r.out$e, as.numeric(mat.out$E))
  expect_equal(r.out$modules, lapply(mat.out$modules, as.r.module.long))
  expect_equal(r.out$break.info, as.r.break.info.chrom(mat.out$breakInfo))
})

test_that("aggregate.segments has the same output as the model", {
  mat.in <- readMat('../referenceio/1/aggregateSegmentsBlock.in.mat')
  mat.out <- readMat('../referenceio/1/aggregateSegmentsBlock.out.mat')
  
  mat.maploc <- read.rubic.maploc(mat.in)
  map.loc.agr <- sum.map.loc(mat.maploc)
  
  r.out <- aggregate.segments(mat.in$data[,,1]$cna,
                              map.loc.agr,
                              as.numeric(mat.in$data[,,1]$ampLevel),
                              as.numeric(mat.in$data[,,1]$delLevel),
                              list(kws=as.vector(mat.in$paramsP[,,1]$kws),
                                   mus=mat.in$paramsP[,,1]$mus,
                                   vars=mat.in$paramsP[,,1]$vars, 
                                   int.sqrt.vs=mat.in$paramsP[,,1]$intSqrtVs,
                                   agr.cum.iter=mat.in$paramsP[,,1]$agrCumIter),
                              list(kws=as.vector(mat.in$paramsN[,,1]$kws),
                                   mus=mat.in$paramsN[,,1]$mus,
                                   vars=mat.in$paramsN[,,1]$vars, 
                                   int.sqrt.vs=mat.in$paramsN[,,1]$intSqrtVs,
                                   agr.cum.iter=mat.in$paramsN[,,1]$agrCumIter),
                              as.numeric(mat.in$FDR))
  
  # reorder the individual modules by name according to the matlab convention
  for (i in 1:length(r.out$segments.p)) {
    r.out$segments.p[[i]] <- r.out$segments.p[[i]][c("I","kw","probes.ignore","l","r","amp","p")]
  }
  for (i in 1:length(r.out$segments.n)) {
    r.out$segments.n[[i]] <- r.out$segments.n[[i]][c("I","kw","probes.ignore","l","r","amp","p")]
  }
  
  expect_equal(r.out$e.p, as.numeric(mat.out$EP))
  expect_equal(r.out$e.n, as.numeric(mat.out$EN))
  expect_equal(r.out$segments.p, lapply(mat.out$segmentsP, as.r.module.long))
  expect_equal(r.out$segments.n, lapply(mat.out$segmentsN, as.r.module.long))
})
