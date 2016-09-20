library(R.matlab, quietly=T, warn.conflict=F)

context("Focal Events - Case 0")

test_that("read.genes.info.tsv has the same output as the model", {
  mat.out <- readMat('../referenceio/0/getGeneInfoFromTxt.out.mat')
  
  r.out <- read.genes.info.tsv('../referenceio/hg19.tsv')
  
  expect_equal(r.out, read.rubic.geneinfo(mat.out))
})


test_that("gene.locs has the same output as the model", {
  mat.in <- readMat('../referenceio/0/getLocs.in.mat')
  mat.out <- readMat('../referenceio/0/getLocs.out.mat')
  mat.maploc <- read.rubic.maploc(mat.in)
  mat.maploc.agr <- sum.map.loc(mat.maploc)
  
  markers <- read.markers('../referenceio/markers.tsv', header=F)
  
  r.out <- gene.locs(mat.maploc.agr,
                     lapply(mat.in$calledEvents, as.r.module.long),
                     markers=markers)
  
  expect_equal(r.out, lapply(mat.out$calledEvents, as.r.focal.event))
})


test_that("filter.overlaps has the same output as the model", {
  mat.in <- readMat('../referenceio/0/filterOverlaps.in.mat')
  mat.out <- readMat('../referenceio/0/filterOverlaps.out.mat')
  
  r.out <- filter.overlaps(lapply(mat.in$focalEvents, as.r.focal.event))
  
  expect_equal(r.out, lapply(mat.out$focalEvents, as.r.focal.event))
})


test_that("called.to.genes has the same output as the model", {
  mat.in <- readMat('../referenceio/0/called2Genes.in.mat')
  mat.out <- readMat('../referenceio/0/called2Genes.out.mat')
  mat.genes.info <- read.rubic.geneinfo(mat.in)
  
  mat.maploc <- read.rubic.maploc(mat.in)
  mat.maploc.agr <- sum.map.loc(mat.maploc)
  
  markers <- read.markers('../referenceio/markers.tsv', header=F)
  
  r.out <- called.to.genes(mat.maploc.agr,
                           lapply(mat.in$calledEvents, as.r.module.long),
                           as.numeric(mat.in$maxLength),
                           mat.genes.info,
                           markers=markers)
  
  r.out <- r.out[order(sapply(r.out, function(x) x$I))]
  
  mat.focal.events <- lapply(mat.out$focalEvents, as.r.focal.event.long)
  mat.focal.events <- mat.focal.events[order(sapply(mat.focal.events, function(x) x$I))]
  
  expect_equal(r.out, mat.focal.events)
})


test_that("calc.break.qvalues has the same output as the model", {
  mat.in <- readMat('../referenceio/0/calcBreakQvalues.in.mat')
  mat.out <- readMat('../referenceio/0/calcBreakQvalues.out.mat')
  
  r.out <- calc.break.qvalues(lapply(mat.in$focalPEvents, as.r.focal.event.long),
                              lapply(mat.in$focalNEvents, as.r.focal.event.long))
  
  expect_equal(r.out$focal.p.events, lapply(mat.out$focalPEvents, as.r.focal.event.final))
  expect_equal(r.out$focal.n.events, lapply(mat.out$focalNEvents, as.r.focal.event.final))
  expect_equal(r.out$q.all, as.numeric(mat.out$qAll))
})


test_that("sort.regions.on.genome has the same output as the model", {
  mat.in <- readMat('../referenceio/0/sortRegionsOnGenome.in.mat')
  mat.out <- readMat('../referenceio/0/sortRegionsOnGenome.out.mat')
  
  r.out <- sort.regions.on.genome(lapply(mat.in$focalEvents, as.r.focal.event.final))
  
  expect_equal(r.out, lapply(mat.out$focalEvents, as.r.focal.event.final))
})

context("Focal Events - Case 1")

test_that("read.genes.info.tsv has the same output as the model", {
  mat.out <- readMat('../referenceio/1/getGeneInfoFromTxt.out.mat')
  
  r.out <- read.genes.info.tsv('../referenceio/hg19.tsv')
  
  expect_equal(r.out, read.rubic.geneinfo(mat.out))
})


test_that("gene.locs has the same output as the model", {
  mat.in <- readMat('../referenceio/1/getLocs.in.mat')
  mat.out <- readMat('../referenceio/1/getLocs.out.mat')
  mat.maploc <- read.rubic.maploc(mat.in)
  mat.maploc.agr <- sum.map.loc(mat.maploc)
  
  markers <- read.markers('../referenceio/markers.tsv', header=F)
  
  r.out <- gene.locs(mat.maploc.agr,
                     lapply(mat.in$calledEvents, as.r.module.long),
                     markers=markers)
  
  expect_equal(r.out, lapply(mat.out$calledEvents, as.r.focal.event))
})


test_that("filter.overlaps has the same output as the model", {
  mat.in <- readMat('../referenceio/1/filterOverlaps.in.mat')
  mat.out <- readMat('../referenceio/1/filterOverlaps.out.mat')
  
  r.out <- filter.overlaps(lapply(mat.in$focalEvents, as.r.focal.event))
  
  expect_equal(r.out, lapply(mat.out$focalEvents, as.r.focal.event))
})


test_that("called.to.genes has the same output as the model", {
  mat.in <- readMat('../referenceio/1/called2Genes.in.mat')
  mat.out <- readMat('../referenceio/1/called2Genes.out.mat')
  mat.genes.info <- read.rubic.geneinfo(mat.in)
  
  mat.maploc <- read.rubic.maploc(mat.in)
  mat.maploc.agr <- sum.map.loc(mat.maploc)
  
  markers <- read.markers('../referenceio/markers.tsv', header=F)
  
  r.out <- called.to.genes(mat.maploc.agr,
                           lapply(mat.in$calledEvents, as.r.module.long),
                           as.numeric(mat.in$maxLength),
                           mat.genes.info,
                           markers=markers)
  
  r.out <- r.out[order(sapply(r.out, function(x) x$I))]
  
  mat.focal.events <- lapply(mat.out$focalEvents, as.r.focal.event.long)
  mat.focal.events <- mat.focal.events[order(sapply(mat.focal.events, function(x) x$I))]
  
  expect_equal(r.out, mat.focal.events)
})


test_that("calc.break.qvalues has the same output as the model", {
  mat.in <- readMat('../referenceio/1/calcBreakQvalues.in.mat')
  mat.out <- readMat('../referenceio/1/calcBreakQvalues.out.mat')
  
  r.out <- calc.break.qvalues(lapply(mat.in$focalPEvents, as.r.focal.event.long),
                              lapply(mat.in$focalNEvents, as.r.focal.event.long))
  
  expect_equal(r.out$focal.p.events, lapply(mat.out$focalPEvents, as.r.focal.event.final))
  expect_equal(r.out$focal.n.events, lapply(mat.out$focalNEvents, as.r.focal.event.final))
  expect_equal(r.out$q.all, as.numeric(mat.out$qAll))
})


test_that("sort.regions.on.genome has the same output as the model", {
  mat.in <- readMat('../referenceio/1/sortRegionsOnGenome.in.mat')
  mat.out <- readMat('../referenceio/1/sortRegionsOnGenome.out.mat')
  
  r.out <- sort.regions.on.genome(lapply(mat.in$focalEvents, as.r.focal.event.final))
  
  expect_equal(r.out, lapply(mat.out$focalEvents, as.r.focal.event.final))
})

