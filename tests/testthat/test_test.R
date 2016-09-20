# test_that("RUBIC preprocess.map.loc has the same output as the model", {
#   samples <- setNames(fread('C:/Users/Nicola/Documents/RUBIC Matlab/source/RUBIC_v0_0/exampleDataLarge/samples.tsv', header=F), 'Sample')
#   
#   r.out <- preprocess.map.loc(read.seg.file('C:/Users/Nicola/Documents/RUBIC Matlab/source/RUBIC_v0_0/exampleDataLarge/segfile.tsv', log.ratio=6), 
#                               read.markers('C:/Users/Nicola/Documents/RUBIC Matlab/source/RUBIC_v0_0/exampleDataLarge/markers.tsv', header=F), 
#                               samples=samples$Sample, 
#                               min.seg.markers = 1, 
#                               min.mean = -2.5, 
#                               max.mean = 2.5)
#   
#   mat.out <- readMat('C:/Users/Nicola/Documents/RUBIC Matlab/source/RUBIC_v0_0/exampleDataLarge/loadDataCNA.out.mat') 
#   mat.maploc <- read.rubic.maploc(mat.out, chr.levels=TRUE) 
#   
#   expect_equivalent(r.out, mat.maploc)
# })

test_that("RUBIC ensure.min.probes has the same output as the model", {
  samples <- setNames(fread('C:/Users/Nicola/Documents/RUBIC Matlab/source/RUBIC_v0_0/exampleDataLarge/samples.tsv', header=F), 'Sample')
  
  r.out <- preprocess.map.loc(read.seg.file('C:/Users/Nicola/Documents/RUBIC Matlab/source/RUBIC_v0_0/exampleDataLarge/segfile.tsv', log.ratio=6), 
                              read.markers('C:/Users/Nicola/Documents/RUBIC Matlab/source/RUBIC_v0_0/exampleDataLarge/markers_unique_sorted.tsv', header=F), 
                              samples=samples$Sample, 
                              min.seg.markers = 1, 
                              min.mean = -2.5, 
                              max.mean = 2.5)
  
  r.out <- ensure.min.probes(r.out, 
                             read.markers('C:/Users/Nicola/Documents/RUBIC Matlab/source/RUBIC_v0_0/exampleDataLarge/markers_unique_sorted.tsv', header=F))
  
  mat.out <- readMat('C:/Users/Nicola/Documents/RUBIC Matlab/source/RUBIC_v0_0/exampleDataLarge/ensureMinProbes.out.mat')
  mat.maploc <- read.rubic.maploc(mat.out) 
  
  expect_equivalent(r.out, mat.maploc)
})
