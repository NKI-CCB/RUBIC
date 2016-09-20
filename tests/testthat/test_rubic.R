context("RUBIC - generic")

# Autotesters

test.seg.cna <- read.seg.file('../referenceio/segfile.tsv', log.ratio=6)
test.markers <- read.markers('../referenceio/markers.tsv', header=F)
test.genes <- read.genes.info.tsv('../referenceio/hg19.tsv')
test.samples <- c('S0001', 'S0005')
test.rubic <- rubic(0.25, test.seg.cna, test.markers, test.samples)

test_that("RUBIC build reads samples file", {
  rbc <- rubic(0.25, test.seg.cna, test.markers, '../referenceio/0/samples.tsv')
  expect_equal(length(rbc$samples), 40)
})


test_that("RUBIC build cannot read sample file", {
  expect_error(rubic(0.25, test.seg.cna, test.markers, 'badidea.tsv'))
})


test_that("RUBIC build use samples vector", {
  rbc <- rubic(0.25, test.seg.cna, test.markers, test.samples)
  expect_equal(rbc$samples, test.samples)
  expect_equal(unique(rbc$map.loc[,Sample]), test.samples)
})


test_that("RUBIC build fill in genes locations from file", {
  rbc <- rubic(0.25, test.seg.cna, test.markers, test.samples,
               genes='../referenceio/hg19.tsv')
  expect_equal(NROW(rbc$genes), 40430)
})

test_that("RUBIC build fill in genes locations from data.table", {
  rbc <- rubic(0.25, test.seg.cna, test.markers, test.samples,
               genes=test.genes)
  expect_equal(NROW(rbc$genes), 40430)
})

test_that("RUBIC build dont fill in genes locations", {
  rbc <- rubic(0.25, test.seg.cna, test.markers, test.samples)
  expect_equal(NROW(rbc$genes), 0)
})


test_that("RUBIC build cant read genes locations", {
  expect_error(rubic(0.25, test.seg.cna, test.markers, test.samples, genes='badidea.tsv'))
})


test_that("RUBIC build creates object with fdr", {
  expect_error(rubic(0xbad1dea, test.seg.cna, test.markers, test.samples))
  expect_error(rubic(-0xbad1dea, test.seg.cna, test.markers, test.samples))
  rbc <- rubic(0.0001, test.seg.cna, test.markers, test.samples)
  expect_equal(rbc$fdr, 0.0001)
})


test_that("RUBIC build creates object with focal.threshold", {
  expect_error(rubic(0.25, test.seg.cna, test.markers, test.samples,
                     focal.threshold=0))
  expect_error(rubic(0.25, test.seg.cna, test.markers, test.samples,
                     focal.threshold=-0xbad1dea))
  rbc <- rubic(0.25, test.seg.cna, test.markers, test.samples,
               focal.threshold=0xbad1dea)
  expect_equal(rbc$focal.threshold, 0xbad1dea)
})


test_that("RUBIC build creates object with del.level", {
  expect_error(rubic(0.25, test.seg.cna, test.markers, test.samples,
                     del.level=0))
  expect_error(rubic(0.25, test.seg.cna, test.markers, test.samples,
                     del.level=0xbad1dea))
  rbc <- rubic(0.25, test.seg.cna, test.markers, test.samples,
               del.level=-0xbad1dea)
  expect_equal(rbc$del.level, -0xbad1dea)
})


test_that("RUBIC build creates object with amp.level", {
  expect_error(rubic(0.25, test.seg.cna, test.markers, test.samples,
                     amp.level=0))
  expect_error(rubic(0.25, test.seg.cna, test.markers, test.samples,
                     amp.level=-0xbad1dea))
  rbc <- rubic(0.25, test.seg.cna, test.markers, test.samples,
               amp.level=0xbad1dea)
  expect_equal(rbc$amp.level, 0xbad1dea)
})


test_that("RUBIC build creates object with min.seg.markers", {
  expect_error(rubic(0.25, test.seg.cna, test.markers, test.samples,
                     min.seg.markers=0))
  expect_error(rubic(0.25, test.seg.cna, test.markers, test.samples,
                     min.seg.markers=-0xbad1dea))
  rbc <- rubic(0.25, test.seg.cna, test.markers, test.samples,
               min.seg.markers=0xbad1dea)
  expect_equal(rbc$min.seg.markers, 0xbad1dea)
})


test_that("RUBIC build creates object with min.probes", {
  expect_error(rubic(0.25, test.seg.cna, test.markers, test.samples,
                     min.probes=0))
  expect_error(rubic(0.25, test.seg.cna, test.markers, test.samples,
                     min.probes=-0xbad1dea))
  rbc <- rubic(0.25, test.seg.cna, test.markers, test.samples,
               min.probes=0xbad1dea)
  expect_equal(rbc$min.probes, 0xbad1dea)
})


test_that("RUBIC build creates object with max/min.mean", {
  rbc <- rubic(0.25, test.seg.cna, test.markers, test.samples,
               max.mean=0xbad1dea, min.mean=0x2bad1dea)
  expect_equal(rbc$max.mean, 0xbad1dea)
  expect_equal(rbc$min.mean, 0x2bad1dea)
})


test_that("RUBIC build creates object (2 samples)", {
  rbc <- rubic(0.25, '../referenceio/segfile.tsv',
               '../referenceio/markers.tsv', test.samples,
               col.log.ratio=6, markers.header=F)
  expect_true(is.na(rbc$min.mean))
  expect_true(is.na(rbc$max.mean))
  expect_equal(rbc$min.probes, 2.6e5)
  expect_equal(rbc$min.seg.markers, 1)
  expect_equal(rbc$amp.level, 0.1)
  expect_equal(rbc$del.level, -0.1)
  expect_equal(rbc$focal.threshold, 10e6)
  expect_equal(rbc$fdr, 0.25)
  expect_equal(unique(rbc$map.loc[,Sample]), test.samples)
  expect_equal(NROW(rbc$map.loc), 20000)
  expect_equal(NROW(rbc$markers), 10000)
  expect_equal(NROW(rbc$genes), 0)
})


test_that("RUBIC build creates object with seg.cna data.table", {
  seg.cna <- read.seg.file('../referenceio/segfile.tsv', log.ratio=6)
  rbc <- rubic(0.25, seg.cna, '../referenceio/markers.tsv', test.samples,
               markers.header=F)
  expect_true(is.na(rbc$min.mean))
  expect_true(is.na(rbc$max.mean))
  expect_equal(rbc$min.probes, 2.6e5)
  expect_equal(rbc$min.seg.markers, 1)
  expect_equal(rbc$amp.level, 0.1)
  expect_equal(rbc$del.level, -0.1)
  expect_equal(rbc$focal.threshold, 10e6)
  expect_equal(rbc$fdr, 0.25)
  expect_equal(unique(rbc$map.loc[,Sample]), test.samples)
  expect_equal(NROW(rbc$map.loc), 20000)
  expect_equal(NROW(rbc$markers), 10000)
})


test_that("RUBIC build creates object with markers data.table", {
  markers <- read.markers('../referenceio/markers.tsv', header=F)
  rbc <- rubic(0.25, '../referenceio/segfile.tsv',
               markers, test.samples, col.log.ratio=6)
  expect_true(is.na(rbc$min.mean))
  expect_true(is.na(rbc$max.mean))
  expect_equal(rbc$min.probes, 2.6e5)
  expect_equal(rbc$min.seg.markers, 1)
  expect_equal(rbc$amp.level, 0.1)
  expect_equal(rbc$del.level, -0.1)
  expect_equal(rbc$focal.threshold, 10e6)
  expect_equal(rbc$fdr, 0.25)
  expect_equal(unique(rbc$map.loc[,Sample]), test.samples)
  expect_equal(NROW(rbc$map.loc), 20000)
  expect_equal(NROW(rbc$markers), 10000)
})


test_that("RUBIC build fails if it cant find files", {
  expect_error(rubic(0.25, 'badidea.tsv', test.markers, test.samples))
  expect_error(rubic(0.25, test.seg.cna, 'badidea.tsv', test.samples))
})


test_that("RUBIC build fails if it cant find data", {
  expect_error(rubic(0.25, test.seg.cna, test.markers, 'badidea'))
})


test_that("RUBIC estimate parameters", {
  rubic <- test.rubic$copy(shallow=TRUE)
  expect_true(isempty(rubic$params.p))
  expect_true(isempty(rubic$params.n))
  expect_true(isempty(rubic$cna.matrix))
  expect_true(isempty(rubic$map.loc.agr))
  rubic$estimate.parameters()
  expect_false(isempty(rubic$cna.matrix))
  expect_false(isempty(rubic$map.loc.agr))
  expect_false(isempty(rubic$params.p))
  expect_false(isempty(rubic$params.n))
})


test_that("RUBIC segments", {
  rubic <- test.rubic$copy(shallow=TRUE)
  rubic$estimate.parameters()
  expect_true(isempty(rubic$segments.p))
  expect_true(isempty(rubic$segments.n))
  rubic$segment()
  expect_false(is.na(rubic$e.p))
  expect_false(is.na(rubic$e.n))
  expect_false(isempty(rubic$segments.p))
  expect_false(isempty(rubic$segments.n))
})


test_that("RUBIC segments at point blank", {
  rubic <- test.rubic$copy(shallow=TRUE)
  expect_true(isempty(rubic$params.p))
  expect_true(isempty(rubic$params.n))
  expect_true(isempty(rubic$segments.p))
  expect_true(isempty(rubic$segments.n))
  rubic$segment()
  expect_false(isempty(rubic$params.p))
  expect_false(isempty(rubic$params.n))
  expect_false(is.na(rubic$e.p))
  expect_false(is.na(rubic$e.n))
  expect_false(isempty(rubic$segments.p))
  expect_false(isempty(rubic$segments.n))
})


test_that("RUBIC call.events", {
  rubic <- test.rubic$copy(shallow=TRUE)
  rubic$estimate.parameters()
  expect_false(isempty(rubic$params.p))
  expect_false(isempty(rubic$params.n))
  rubic$segment()
  expect_false(isempty(rubic$segments.p))
  expect_false(isempty(rubic$segments.n))
  expect_true(isempty(rubic$called.p.events))
  expect_true(isempty(rubic$called.n.events))
  rubic$call.events()
  expect_false(isempty(rubic$called.p.events))
  expect_false(isempty(rubic$called.n.events))
})


test_that("RUBIC call.events at point blank", {
  rubic <- test.rubic$copy(shallow=TRUE)
  expect_true(isempty(rubic$called.p.events))
  expect_true(isempty(rubic$called.n.events))
  rubic$call.events()
  expect_false(isempty(rubic$params.p))
  expect_false(isempty(rubic$params.n))
  expect_false(isempty(rubic$called.p.events))
  expect_false(isempty(rubic$called.n.events))
})


test_that("RUBIC call.focal.events", {
  rubic <- test.rubic$copy(shallow=TRUE)
  rubic$estimate.parameters()
  expect_false(isempty(rubic$params.p))
  expect_false(isempty(rubic$params.n))
  rubic$segment()
  expect_false(isempty(rubic$segments.p))
  expect_false(isempty(rubic$segments.n))
  expect_true(isempty(rubic$called.p.events))
  # expect_true(isempty(rubic$called.n.events))
  rubic$call.events()
  expect_false(isempty(rubic$called.p.events))
  # expect_false(isempty(rubic$called.n.events))
  expect_true(isempty(rubic$focal.p.events))
  # expect_true(isempty(rubic$focal.n.events))
  rubic$call.focal.events()
  expect_false(isempty(rubic$focal.p.events))
  # expect_false(isempty(rubic$focal.n.events))
})


test_that("RUBIC call.focal.events uses downloaded genes locations", {
  rubic <- test.rubic$copy(shallow=TRUE)
  expect_equal(NROW(rubic$genes), 0)
  #genes=test.genes
  expect_true(isempty(rubic$focal.p.events))
  expect_true(isempty(rubic$focal.n.events))
  rubic$call.focal.events()
  expect_true(NROW(rubic$genes) > 0)
  expect_false(isempty(rubic$focal.p.events))
  # expect_false(isempty(rubic$focal.n.events))
})


test_that("RUBIC call.focal.events uses custom genes locations", {
  rubic <- test.rubic$copy(shallow=TRUE)
  expect_equal(NROW(rubic$genes), 0)
  expect_true(isempty(rubic$focal.p.events))
  expect_true(isempty(rubic$focal.n.events))
  rubic$call.focal.events(test.genes[1:20000])
  expect_true(NROW(rubic$genes) == 20000)
  expect_false(isempty(rubic$focal.p.events))
  # expect_false(isempty(rubic$focal.n.events))
})


test_that("RUBIC call.focal.events uses preloaded genes locations", {
  rubic <- rubic(0.25, test.seg.cna, test.markers, test.samples,
                 genes='../referenceio/hg19.tsv')
  expect_equal(NROW(rubic$genes), 40430)
  expect_true(isempty(rubic$focal.p.events))
  expect_true(isempty(rubic$focal.n.events))
  rubic$call.focal.events()
  expect_equal(NROW(rubic$genes), 40430)
  expect_false(isempty(rubic$focal.p.events))
  # expect_false(isempty(rubic$focal.n.events))
})