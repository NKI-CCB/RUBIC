library(R.matlab, quietly=T, warn.conflict=F)

context("Preprocessing - Case 0")

test_that("RUBIC preprocess.map.loc has the same output as the model", {
  samples <- setNames(fread('../referenceio/0/samples.tsv', header=F), 'Sample')

  r.out <- preprocess.map.loc(read.seg.file('../referenceio/segfile.tsv', log.ratio=6), 
                              read.markers('../referenceio/markers.tsv', header=F), 
                              samples=samples$Sample, 
                              min.seg.markers = 1, 
                              min.mean = -2.5, 
                              max.mean = 2.5)
  
  mat.out <- readMat('../referenceio/0/loadDataCNA.out.mat') 
  mat.maploc <- read.rubic.maploc(mat.out, chr.levels=TRUE) 
  
  expect_equivalent(r.out, mat.maploc)
})


test_that("RUBIC ensure.min.probes has the same output as the model", {
  mat.in <- readMat('../referenceio/0/ensureMinProbes.in.mat')
  maploc.in <- read.rubic.maploc(mat.in)
  
  r.out <- ensure.min.probes(maploc.in, 
                             read.markers('../referenceio/markers.tsv', header=F))
  
  mat.out <- readMat('../referenceio/0/ensureMinProbes.out.mat')
  mat.maploc <- read.rubic.maploc(mat.out) 
  
  expect_equivalent(r.out, mat.maploc)
})


context("Preprocessing - Case 1")

test_that("RUBIC preprocess.map.loc has the same output as the model", {
  samples <- setNames(fread('../referenceio/1/samples.tsv', header=F), 'Sample')
  
  r.out <- preprocess.map.loc(read.seg.file('../referenceio/segfile.tsv', log.ratio=6), 
                              read.markers('../referenceio/markers.tsv', header=F), 
                              samples=samples$Sample, 
                              min.seg.markers = 1, 
                              min.mean = -2.5, 
                              max.mean = 2.5)
  
  mat.out <- readMat('../referenceio/1/loadDataCNA.norep.out.mat') 
  mat.maploc <- read.rubic.maploc(mat.out, chr.levels=TRUE) 
  
  expect_equivalent(r.out, mat.maploc)
})

test_that("RUBIC preprocess.map.loc has the same output as the model when input sample file has repeated sample names", {
  samples <- setNames(fread('../referenceio/1/samples-rep.tsv', header=F), 'Sample')
  
  r.out <- preprocess.map.loc(read.seg.file('../referenceio/segfile.tsv', log.ratio=6), 
                              read.markers('../referenceio/markers.tsv', header=F), 
                              samples=samples$Sample, 
                              min.seg.markers = 1, 
                              min.mean = -2.5, 
                              max.mean = 2.5)
  
  mat.out <- readMat('../referenceio/1/loadDataCNA.rep.out.mat') 
  mat.maploc <- read.rubic.maploc(mat.out, chr.levels=TRUE) 
  
  expect_equivalent(r.out, mat.maploc)
})


test_that("RUBIC ensure.min.probes has the same output as the model", {
  mat.in <- readMat('../referenceio/1/ensureMinProbes.in.mat')
  maploc.in <- read.rubic.maploc(mat.in)
  
  r.out <- ensure.min.probes(maploc.in, 
                             read.markers('../referenceio/markers.tsv', header=F))
  
  mat.out <- readMat('../referenceio/1/ensureMinProbes.out.mat')
  mat.maploc <- read.rubic.maploc(mat.out) 
  
  expect_equivalent(r.out, mat.maploc)
})
