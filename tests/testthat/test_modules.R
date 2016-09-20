library(R.matlab, quietly=T, warn.conflict=F)

context("Modules - Case 0")

test_that("sort.modules  has the same output as the model", {
  mat.in <- readMat('../referenceio/0/sortModules.in.mat')
  mat.out <- readMat('../referenceio/0/sortModules.out.mat')
  
  r.out <- sort.modules(lapply(mat.in$modules, function(x) as.r.module(x[[1]])),
                        as.r.chr(mat.in$chroms, levels=TRUE))
  
  expect_equal(r.out$modules, lapply(mat.out$modules, function(x) as.r.module(x[[1]])))
  expect_equal(r.out$chromosomes, as.r.chr(mat.out$modChroms, levels=TRUE))  
})


test_that("base.modules  has the same output as the model", {
  mat.in <- readMat('../referenceio/0/getBaseModules.in.mat')
  mat.in.maploc <- read.rubic.maploc(mat.in)
  mat.out <- readMat('../referenceio/0/getBaseModules.out.mat')
  
  r.out <- base.modules(sum.map.loc(mat.in.maploc))
  
  expect_equal(r.out$modules, lapply(mat.out$modules, function(x) as.r.module(x[[1]])))
  expect_equal(r.out$chromosomes, as.r.chr(mat.out$chroms, levels=TRUE))  
})


test_that("trunk.mods  has the same output as the model", {
  mat.in <- readMat('../referenceio/0/getTrunkMods.in.mat')
  mat.out <- readMat('../referenceio/0/getTrunkMods.out.mat')
    
  r.out <- trunk.mods(as.r.module(mat.in$module1), 
                      as.r.module(mat.in$module2), 
                      as.vector(mat.in$params[,,1]$kws))
  
  expect_equivalent(r.out$mod1.probes, as.vector(mat.out$mod1Probes))
  expect_equivalent(r.out$mod2.probes, as.vector(mat.out$mod2Probes))
})


test_that("params.at.kws has the same output as the model", {
  mat.in <- readMat('../referenceio/0/getParamsAtKws.in.mat')
  mat.out <- readMat('../referenceio/0/getParamsAtKws.out.mat')
  
  r.out <- params.at.kws(as.numeric(mat.in$kw1),
                         as.numeric(mat.in$kw2),
                         list(kws=as.vector(mat.in$params[,,1]$kws),
                              mus=mat.in$params[,,1]$mus,
                              vars=mat.in$params[,,1]$vars, 
                              int.sqrt.vs=mat.in$params[,,1]$intSqrtVs))
  
  expect_equal(r.out$mu.est, as.numeric(mat.out$muEst))
  expect_equal(r.out$var.est, as.numeric(mat.out$varEst))
  expect_equal(r.out$int.sqrt.v.est, as.numeric(mat.out$intSqrtVEst))
})


test_that("anal.neg.log.e.euler has the same output as the model", {
  mat.in <- readMat('../referenceio/0/getAnalNegLogEEuler.in.mat')
  mat.out <- readMat('../referenceio/0/getAnalNegLogEEuler.out.mat')
  
  r.out <- anal.neg.log.e.euler(as.numeric(mat.in$intSqrtV),
                                  as.numeric(mat.in$eulerNum),
                                  as.numeric(mat.in$t),
                                  as.numeric(mat.in$numTails))
  
  expect_equal(r.out$e.euler, as.numeric(mat.out$eEuler))
  expect_equal(r.out$e.max, as.numeric(mat.out$eMax))
})


test_that("break.sig.v2 has the same output as the model", {
  mat.in <- readMat('../referenceio/0/getBreakSigV2.in.mat')
  mat.out <- readMat('../referenceio/0/getBreakSigV2.out.mat')
  
  r.out <- break.sig.v2(mat.in$cna,
                         as.r.module(mat.in$module1),
                         as.r.module(mat.in$module2),
                         list(kws=as.vector(mat.in$params[,,1]$kws),
                              mus=mat.in$params[,,1]$mus,
                              vars=mat.in$params[,,1]$vars, 
                              int.sqrt.vs=mat.in$params[,,1]$intSqrtVs))
  
  expect_equal(r.out$reg.start, as.numeric(mat.out$breakInfo[,,1]$regStart))
  expect_equal(r.out$reg.end, as.numeric(mat.out$breakInfo[,,1]$regEnd))
  expect_equal(r.out$kw1, as.numeric(mat.out$breakInfo[,,1]$kw1))
  expect_equal(r.out$kw2, as.numeric(mat.out$breakInfo[,,1]$kw2))
  expect_equal(r.out$t, as.numeric(mat.out$breakInfo[,,1]$t))
  expect_equal(r.out$z, as.numeric(mat.out$breakInfo[,,1]$z))
  expect_equal(r.out$p, as.numeric(mat.out$breakInfo[,,1]$p))
})


test_that("adj.sig has the same output as the model", {
  mat.in <- readMat('../referenceio/0/getSigAdj.in.mat')
  mat.out <- readMat('../referenceio/0/getSigAdj.out.mat')
  
  r.out <- adj.sig(mat.in$cna,
                  lapply(mat.in$modules, function(x) as.r.module(x[[1]])),
                  as.numeric(mat.in$checkI), 
                  list(kws=as.vector(mat.in$params[,,1]$kws),
                       mus=mat.in$params[,,1]$mus,
                       vars=mat.in$params[,,1]$vars,
                       int.sqrt.vs=mat.in$params[,,1]$intSqrtVs))
  
  expect_equal(r.out$reg.start, as.numeric(mat.out$breakInfo[,,1]$regStart))
  expect_equal(r.out$reg.end, as.numeric(mat.out$breakInfo[,,1]$regEnd))
  expect_equal(r.out$kw1, as.numeric(mat.out$breakInfo[,,1]$kw1))
  expect_equal(r.out$kw2, as.numeric(mat.out$breakInfo[,,1]$kw2))
  expect_equal(r.out$t, as.numeric(mat.out$breakInfo[,,1]$t))
  expect_equal(r.out$z, as.numeric(mat.out$breakInfo[,,1]$z))
  expect_equal(r.out$p, as.numeric(mat.out$breakInfo[,,1]$p))
})


test_that("adj.sigs has the same output as the model", {
  mat.in <- readMat('../referenceio/0/getAdjSigs.in.mat')
  mat.out <- readMat('../referenceio/0/getAdjSigs.out.mat')
  
  r.out <- adj.sigs(mat.in$cna,
                   lapply(mat.in$modules, function(x) as.r.module(x[[1]])),
                   list(kws=as.vector(mat.in$params[,,1]$kws),
                        mus=mat.in$params[,,1]$mus,
                        vars=mat.in$params[,,1]$vars,
                        int.sqrt.vs=mat.in$params[,,1]$intSqrtVs))
  
  expect_equal(r.out$break.ps, as.vector(mat.out$breakPs))
  expect_equal(r.out$break.info, as.r.break.info(mat.out$breakInfo))
})
  

test_that("cna.cum.to.smooth has the same output as the model", {
  mat.in <- readMat('../referenceio/0/cnaCum2Smooth.in.mat')
  mat.out <- readMat('../referenceio/0/cnaCum2Smooth.out.mat')

  r.out <- cna.cum.to.smooth(mat.in$cnaCum,
                              as.numeric(mat.in$P),
                              as.numeric(mat.in$kw1),
                              as.numeric(mat.in$kw2))
  
  expect_equal(r.out, mat.out$cnaS)
})


test_that("iter.at.t has the same output as the model", {
  mat.in <- readMat('../referenceio/0/iterAtT.in.mat')
  mat.out <- readMat('../referenceio/0/iterAtT.out.mat')
  
  r.out <- iter.at.t(mat.in$agrProfilesCum,
                     as.numeric(mat.in$P),
                     as.numeric(mat.in$t),
                     as.numeric(mat.in$kw1),
                     as.numeric(mat.in$kw2))
  
  expect_equal(r.out, as.numeric(mat.out$Eemp))
})


test_that("updated.adj.sigs has the same output as the model", {
  mat.in <- readMat('../referenceio/0/updatedAdjSigs.in.mat')
  mat.out <- readMat('../referenceio/0/updatedAdjSigs.out.mat')
  
  r.out <- updated.adj.sigs(mat.in$cna,
                   lapply(mat.in$modules, function(x) as.r.module(x[[1]])),
                   as.numeric(mat.in$mergeI),
                   as.r.break.info(mat.in$breakInfo),
                   as.vector(mat.in$breakPs),
                   list(kws=as.vector(mat.in$params[,,1]$kws),
                        mus=mat.in$params[,,1]$mus,
                        vars=mat.in$params[,,1]$vars,
                        int.sqrt.vs=mat.in$params[,,1]$intSqrtVs))
  
  expect_equal(r.out$break.ps, as.vector(mat.out$breakPs))
  expect_equal(r.out$break.info, as.r.break.info(mat.out$breakInfo))
})


test_that("merge.adj.modules has the same output as the model", {
  mat.in <- readMat('../referenceio/0/mergeAdjModules.in.mat')
  mat.out <- readMat('../referenceio/0/mergeAdjModules.out.mat')

  r.out <- merge.adj.modules(mat.in$cna,
                             lapply(mat.in$modules, function(x) as.r.module(x[[1]])),
                             as.numeric(mat.in$mergeI),
                             as.r.break.info(mat.in$breakInfo),
                             as.vector(mat.in$breakPs),
                             list(kws=as.vector(mat.in$params[,,1]$kws),
                                  mus=mat.in$params[,,1]$mus,
                                  vars=mat.in$params[,,1]$vars,
                                  int.sqrt.vs=mat.in$params[,,1]$intSqrtVs))
  
  expect_equal(r.out$modules, lapply(mat.out$modules, function(x) as.r.module(x[[1]])))
  expect_equal(r.out$break.ps, as.vector(mat.out$breakPs))
  expect_equal(r.out$break.info, as.r.break.info(mat.out$breakInfo))         
})

test_that("mod.sigs has the same output as the model", {
  mat.in <- readMat('../referenceio/0/getModSigs.in.mat')
  mat.out <- readMat('../referenceio/0/getModSigs.out.mat')
  
  r.out <- mod.sigs(lapply(mat.in$modules, function(x) as.r.module(x[[1]])),
                    as.r.break.info(mat.in$breakInfo))
  
  mat.modules <- lapply(mat.out$modules, as.r.module.long)
  for (i in 1:length(mat.modules)) {
    mat.modules[[i]]$p <- NULL
  }
  
  # reorder the individual modules by name according to the matlab convention
  for (i in 1:length(r.out)) {
    r.out[[i]] <- r.out[[i]][c("I","kw","probes.ignore","l","r","amp","p")]
    r.out[[i]]$p <- NULL
  }
  
  expect_equal(r.out, mat.modules)
})

test_that("merge.modules.v2 has the same output as the model", {
  mat.in <- readMat('../referenceio/0/mergeModulesV2.in.mat')
  mat.out <- readMat('../referenceio/0/mergeModulesV2.out.mat')
  
  map.loc <- read.rubic.maploc(mat.in)
  map.loc.agr <- sum.map.loc(map.loc)
  
  # MATLAB call in iterExtractSingleSampParams      
  # [modules, breakInfo, sigMods] = mergeModulesV2(data, sampsUse, [], sigModsPrev, params, -log10(pValue));
  #
  # We collected i/o at this point raher than inside mergeModulesV2 itself
  # since the MATLAB implementation is recursive and would have made
  # difficult to collect output matching to input for the initial call
  
  r.out <- merge.modules.v2(mat.in$data[,,1]$cna,
                            map.loc.agr,
                            as.numeric(mat.in$data[,,1]$P),
                            as.vector(mat.in$sampsUse),
                            matrix(nrow=0,ncol=0),
                            as.list(mat.in$sigModPrev), 
                            mat.in$params[,,1]$agrCumIter,
                            list(kws=as.vector(mat.in$params[,,1]$kws),
                                 mus=mat.in$params[,,1]$mus,
                                 vars=mat.in$params[,,1]$vars, 
                                 int.sqrt.vs=mat.in$params[,,1]$intSqrtVs),
                            -log10(as.numeric(mat.in$pValue))) 
  
  
  mat.modules <- lapply(mat.out$modules, as.r.module.long)
  for (i in 1:length(mat.modules)) {
    mat.modules[[i]]$p <- NULL
  }
  
  # reorder the individual modules by name according to the matlab convention
  for (i in 1:length(r.out$modules)) {
    r.out$modules[[i]] <- r.out$modules[[i]][c("I","kw","probes.ignore","l","r","amp","p")]
    r.out$modules[[i]]$p <- NULL
  }
  
  expect_equal(r.out$modules, mat.modules)
  expect_equal(r.out$break.info, as.r.break.info.chrom(mat.out$breakInfo))    
  expect_equal(r.out$sig.mods$starts, as.vector(mat.out$sigMods[,,1]$starts))
  expect_equal(r.out$sig.mods$ends, as.vector(mat.out$sigMods[,,1]$ends))
})


test_that("merge.modules.v2 has the same output as the model with usePerm == FALSE", {
  mat.in <- readMat('../referenceio/0/mergeModulesV2.noperm.in.mat')
  mat.out <- readMat('../referenceio/0/mergeModulesV2.noperm.out.mat')
  
  mat.map.loc <- read.rubic.maploc(mat.in)
  mat.map.loc.agr <- sum.map.loc(mat.map.loc)
  
  r.out <- merge.modules.v2(mat.in$data[,,1]$cna,
                            mat.map.loc.agr,
                            as.numeric(mat.in$data[,,1]$P),
                            as.vector(mat.in$sampsUse),
                            matrix(nrow=0,ncol=0),
                            as.list(mat.in$sigModPrev), 
                            mat.in$params[,,1]$agrCumIter,
                            list(kws=as.vector(mat.in$params[,,1]$kws),
                                 mus=mat.in$params[,,1]$mus,
                                 vars=mat.in$params[,,1]$vars, 
                                 int.sqrt.vs=mat.in$params[,,1]$intSqrtVs),
                            -log10(0.25),
                            use.perm=F)
  
  mat.modules <- lapply(mat.out$modules, as.r.module.long)
  for (i in 1:length(mat.modules)) {
    mat.modules[[i]]$p <- NULL
  }
  
  # reorder the individual modules by name according to the matlab convention
  for (i in 1:length(r.out$modules)) {
    r.out$modules[[i]] <- r.out$modules[[i]][c("I","kw","probes.ignore","l","r","amp","p")]
    r.out$modules[[i]]$p <- NULL
  }
  
  expect_equal(r.out$modules, mat.modules)
})


context("Modules - Case 1")

test_that("sort.modules  has the same output as the model", {
  mat.in <- readMat('../referenceio/1/sortModules.in.mat')
  mat.out <- readMat('../referenceio/1/sortModules.out.mat')
  
  r.out <- sort.modules(lapply(mat.in$modules, function(x) as.r.module(x[[1]])),
                        as.r.chr(mat.in$chroms, levels=TRUE))
  
  expect_equal(r.out$modules, lapply(mat.out$modules, function(x) as.r.module(x[[1]])))
  expect_equal(r.out$chromosomes, as.r.chr(mat.out$modChroms, levels=TRUE))  
})


test_that("base.modules  has the same output as the model", {
  mat.in <- readMat('../referenceio/1/getBaseModules.in.mat')
  mat.in.maploc <- read.rubic.maploc(mat.in)
  mat.out <- readMat('../referenceio/1/getBaseModules.out.mat')
  
  r.out <- base.modules(sum.map.loc(mat.in.maploc))
  
  expect_equal(r.out$modules, lapply(mat.out$modules, function(x) as.r.module(x[[1]])))
  expect_equal(r.out$chromosomes, as.r.chr(mat.out$chroms, levels=TRUE))  
})


test_that("trunk.mods  has the same output as the model", {
  mat.in <- readMat('../referenceio/1/getTrunkMods.in.mat')
  mat.out <- readMat('../referenceio/1/getTrunkMods.out.mat')
  
  r.out <- trunk.mods(as.r.module(mat.in$module1), 
                      as.r.module(mat.in$module2), 
                      as.vector(mat.in$params[,,1]$kws))
  
  expect_equivalent(r.out$mod1.probes, as.vector(mat.out$mod1Probes))
  expect_equivalent(r.out$mod2.probes, as.vector(mat.out$mod2Probes))
})


test_that("params.at.kws has the same output as the model", {
  mat.in <- readMat('../referenceio/1/getParamsAtKws.in.mat')
  mat.out <- readMat('../referenceio/1/getParamsAtKws.out.mat')
  
  r.out <- params.at.kws(as.numeric(mat.in$kw1),
                         as.numeric(mat.in$kw2),
                         list(kws=as.vector(mat.in$params[,,1]$kws),
                              mus=mat.in$params[,,1]$mus,
                              vars=mat.in$params[,,1]$vars, 
                              int.sqrt.vs=mat.in$params[,,1]$intSqrtVs))
  
  expect_equal(r.out$mu.est, as.numeric(mat.out$muEst))
  expect_equal(r.out$var.est, as.numeric(mat.out$varEst))
  expect_equal(r.out$int.sqrt.v.est, as.numeric(mat.out$intSqrtVEst))
  
})


test_that("anal.neg.log.e.euler has the same output as the model", {
  mat.in <- readMat('../referenceio/1/getAnalNegLogEEuler.in.mat')
  mat.out <- readMat('../referenceio/1/getAnalNegLogEEuler.out.mat')
  
  r.out <- anal.neg.log.e.euler(as.numeric(mat.in$intSqrtV),
                                as.numeric(mat.in$eulerNum),
                                as.numeric(mat.in$t),
                                as.numeric(mat.in$numTails))
  
  expect_equal(r.out$e.euler, as.numeric(mat.out$eEuler))
  expect_equal(r.out$e.max, as.numeric(mat.out$eMax))
  
})


test_that("break.sig.v2 has the same output as the model", {
  mat.in <- readMat('../referenceio/1/getBreakSigV2.in.mat')
  mat.out <- readMat('../referenceio/1/getBreakSigV2.out.mat')
  
  r.out <- break.sig.v2(mat.in$cna,
                        as.r.module(mat.in$module1),
                        as.r.module(mat.in$module2),
                        list(kws=as.vector(mat.in$params[,,1]$kws),
                             mus=mat.in$params[,,1]$mus,
                             vars=mat.in$params[,,1]$vars, 
                             int.sqrt.vs=mat.in$params[,,1]$intSqrtVs))
  
  expect_equal(r.out$reg.start, as.numeric(mat.out$breakInfo[,,1]$regStart))
  expect_equal(r.out$reg.end, as.numeric(mat.out$breakInfo[,,1]$regEnd))
  expect_equal(r.out$kw1, as.numeric(mat.out$breakInfo[,,1]$kw1))
  expect_equal(r.out$kw2, as.numeric(mat.out$breakInfo[,,1]$kw2))
  expect_equal(r.out$t, as.numeric(mat.out$breakInfo[,,1]$t))
  expect_equal(r.out$z, as.numeric(mat.out$breakInfo[,,1]$z))
  expect_equal(r.out$p, as.numeric(mat.out$breakInfo[,,1]$p))
})


test_that("adj.sig has the same output as the model", {
  mat.in <- readMat('../referenceio/1/getSigAdj.in.mat')
  mat.out <- readMat('../referenceio/1/getSigAdj.out.mat')
  
  r.out <- adj.sig(mat.in$cna,
                   lapply(mat.in$modules, function(x) as.r.module(x[[1]])),
                   as.numeric(mat.in$checkI), 
                   list(kws=as.vector(mat.in$params[,,1]$kws),
                        mus=mat.in$params[,,1]$mus,
                        vars=mat.in$params[,,1]$vars,
                        int.sqrt.vs=mat.in$params[,,1]$intSqrtVs))
  
  expect_equal(r.out$reg.start, as.numeric(mat.out$breakInfo[,,1]$regStart))
  expect_equal(r.out$reg.end, as.numeric(mat.out$breakInfo[,,1]$regEnd))
  expect_equal(r.out$kw1, as.numeric(mat.out$breakInfo[,,1]$kw1))
  expect_equal(r.out$kw2, as.numeric(mat.out$breakInfo[,,1]$kw2))
  expect_equal(r.out$t, as.numeric(mat.out$breakInfo[,,1]$t))
  expect_equal(r.out$z, as.numeric(mat.out$breakInfo[,,1]$z))
  expect_equal(r.out$p, as.numeric(mat.out$breakInfo[,,1]$p))
})



test_that("adj.sigs has the same output as the model", {
  mat.in <- readMat('../referenceio/1/getAdjSigs.in.mat')
  mat.out <- readMat('../referenceio/1/getAdjSigs.out.mat')
  
  r.out <- adj.sigs(mat.in$cna,
                    lapply(mat.in$modules, function(x) as.r.module(x[[1]])),
                    list(kws=as.vector(mat.in$params[,,1]$kws),
                         mus=mat.in$params[,,1]$mus,
                         vars=mat.in$params[,,1]$vars,
                         int.sqrt.vs=mat.in$params[,,1]$intSqrtVs))
  
  expect_equal(r.out$break.ps, as.vector(mat.out$breakPs))
  expect_equal(r.out$break.info, as.r.break.info(mat.out$breakInfo))
})


test_that("cna.cum.to.smooth has the same output as the model", {
  mat.in <- readMat('../referenceio/1/cnaCum2Smooth.in.mat')
  mat.out <- readMat('../referenceio/1/cnaCum2Smooth.out.mat')
  
  r.out <- cna.cum.to.smooth(mat.in$cnaCum,
                             as.numeric(mat.in$P),
                             as.numeric(mat.in$kw1),
                             as.numeric(mat.in$kw2))
  
  expect_equal(r.out, mat.out$cnaS)
})


test_that("iter.at.t has the same output as the model", {
  mat.in <- readMat('../referenceio/1/iterAtT.in.mat')
  mat.out <- readMat('../referenceio/1/iterAtT.out.mat')
  
  r.out <- iter.at.t(mat.in$agrProfilesCum,
                     as.numeric(mat.in$P),
                     as.numeric(mat.in$t),
                     as.numeric(mat.in$kw1),
                     as.numeric(mat.in$kw2))
  
  expect_equal(r.out, as.numeric(mat.out$Eemp))
})


test_that("updated.adj.sigs has the same output as the model", {
  mat.in <- readMat('../referenceio/1/updatedAdjSigs.in.mat')
  mat.out <- readMat('../referenceio/1/updatedAdjSigs.out.mat')
  
  r.out <- updated.adj.sigs(mat.in$cna,
                            lapply(mat.in$modules, function(x) as.r.module(x[[1]])),
                            as.numeric(mat.in$mergeI),
                            as.r.break.info(mat.in$breakInfo),
                            as.vector(mat.in$breakPs),
                            list(kws=as.vector(mat.in$params[,,1]$kws),
                                 mus=mat.in$params[,,1]$mus,
                                 vars=mat.in$params[,,1]$vars,
                                 int.sqrt.vs=mat.in$params[,,1]$intSqrtVs))
  
  expect_equal(r.out$break.ps, as.vector(mat.out$breakPs))
  expect_equal(r.out$break.info, as.r.break.info(mat.out$breakInfo))
})


test_that("merge.adj.modules has the same output as the model", {
  mat.in <- readMat('../referenceio/1/mergeAdjModules.in.mat')
  mat.out <- readMat('../referenceio/1/mergeAdjModules.out.mat')
  
  r.out <- merge.adj.modules(mat.in$cna,
                             lapply(mat.in$modules, function(x) as.r.module(x[[1]])),
                             as.numeric(mat.in$mergeI),
                             as.r.break.info(mat.in$breakInfo),
                             as.vector(mat.in$breakPs),
                             list(kws=as.vector(mat.in$params[,,1]$kws),
                                  mus=mat.in$params[,,1]$mus,
                                  vars=mat.in$params[,,1]$vars,
                                  int.sqrt.vs=mat.in$params[,,1]$intSqrtVs))
  
  expect_equal(r.out$modules, lapply(mat.out$modules, function(x) as.r.module(x[[1]])))
  expect_equal(r.out$break.ps, as.vector(mat.out$breakPs))
  expect_equal(r.out$break.info, as.r.break.info(mat.out$breakInfo))         
})


test_that("mod.sigs has the same output as the model", {
  mat.in <- readMat('../referenceio/1/getModSigs.in.mat')
  mat.out <- readMat('../referenceio/1/getModSigs.out.mat')
  
  r.out <- mod.sigs(lapply(mat.in$modules, function(x) as.r.module(x[[1]])),
                    as.r.break.info(mat.in$breakInfo))
  
  mat.modules <- lapply(mat.out$modules, as.r.module.long)
  for (i in 1:length(mat.modules)) {
    mat.modules[[i]]$p <- NULL
  }
  
  # reorder the individual modules by name according to the matlab convention
  for (i in 1:length(r.out)) {
    r.out[[i]] <- r.out[[i]][c("I","kw","probes.ignore","l","r","amp","p")]
    r.out[[i]]$p <- NULL
  }
  
  expect_equal(r.out, mat.modules)
})

test_that("merge.modules.v2 has the same output as the model", {
  mat.in <- readMat('../referenceio/1/mergeModulesV2.in.mat')
  mat.out <- readMat('../referenceio/1/mergeModulesV2.out.mat')
  
  map.loc <- read.rubic.maploc(mat.in)
  map.loc.agr <- sum.map.loc(map.loc)
  
  # MATLAB call in iterExtractSingleSampParams      
  # [modules, breakInfo, sigMods] = mergeModulesV2(data, sampsUse, [], sigModsPrev, params, -log10(pValue));
  #
  # We collected i/o at this point raher than inside mergeModulesV2 itself
  # since the MATLAB implementation is recursive and would have made
  # difficult to collect output matching to input for the initial call
  
  r.out <- merge.modules.v2(mat.in$data[,,1]$cna,
                            map.loc.agr,
                            as.numeric(mat.in$data[,,1]$P),
                            as.vector(mat.in$sampsUse),
                            matrix(nrow=0,ncol=0),
                            as.list(mat.in$sigModPrev), 
                            mat.in$params[,,1]$agrCumIter,
                            list(kws=as.vector(mat.in$params[,,1]$kws),
                                 mus=mat.in$params[,,1]$mus,
                                 vars=mat.in$params[,,1]$vars, 
                                 int.sqrt.vs=mat.in$params[,,1]$intSqrtVs),
                            -log10(as.numeric(mat.in$pValue))) 
  
  
  mat.modules <- lapply(mat.out$modules, as.r.module.long)
  for (i in 1:length(mat.modules)) {
    mat.modules[[i]]$p <- NULL
  }
  
  # reorder the individual modules by name according to the matlab convention
  for (i in 1:length(r.out$modules)) {
    r.out$modules[[i]] <- r.out$modules[[i]][c("I","kw","probes.ignore","l","r","amp","p")]
    r.out$modules[[i]]$p <- NULL
  }
  
  expect_equal(r.out$modules, mat.modules)
  expect_equal(r.out$break.info, as.r.break.info.chrom(mat.out$breakInfo))    
  expect_equal(r.out$sig.mods$starts, as.vector(mat.out$sigMods[,,1]$starts))
  expect_equal(r.out$sig.mods$ends, as.vector(mat.out$sigMods[,,1]$ends))
})


test_that("merge.modules.v2 has the same output as the model with usePerm == FALSE", {
  mat.in <- readMat('../referenceio/1/mergeModulesV2.noperm.in.mat')
  mat.out <- readMat('../referenceio/1/mergeModulesV2.noperm.out.mat')
  
  mat.map.loc <- read.rubic.maploc(mat.in)
  mat.map.loc.agr <- sum.map.loc(mat.map.loc)
  
  r.out <- merge.modules.v2(mat.in$data[,,1]$cna,
                            mat.map.loc.agr,
                            as.numeric(mat.in$data[,,1]$P),
                            as.vector(mat.in$sampsUse),
                            matrix(nrow=0,ncol=0),
                            as.list(mat.in$sigModPrev), 
                            mat.in$params[,,1]$agrCumIter,
                            list(kws=as.vector(mat.in$params[,,1]$kws),
                                 mus=mat.in$params[,,1]$mus,
                                 vars=mat.in$params[,,1]$vars, 
                                 int.sqrt.vs=mat.in$params[,,1]$intSqrtVs),
                            -log10(0.25),
                            use.perm=F)
  
  mat.modules <- lapply(mat.out$modules, as.r.module.long)
  for (i in 1:length(mat.modules)) {
    mat.modules[[i]]$p <- NULL
  }
  
  # reorder the individual modules by name according to the matlab convention
  for (i in 1:length(r.out$modules)) {
    r.out$modules[[i]] <- r.out$modules[[i]][c("I","kw","probes.ignore","l","r","amp","p")]
    r.out$modules[[i]]$p <- NULL
  }
  
  expect_equal(r.out$modules, mat.modules)
})

