# Copyright 2015 Netherlands Cancer Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


sort.modules <- function(modules, chromosomes) {
  num.mod <- length(modules)
  starts <- rep(0, num.mod)
  mod.chroms <- vector('character', num.mod)
  for (i in seq_len(num.mod)) {
    starts[i] <- modules[[i]]$I
    mod.chroms[i] <- as.character(chromosomes[starts[i]])
  }
  sorted.is <- order(starts)
  mod.chroms <- ordered(mod.chroms, levels=levels(chromosomes))[,drop=T]
  list(modules=modules[sorted.is], chromosomes=mod.chroms[sorted.is])
}


#'
#' @param map.loc.agr The sum of the log ratio as returned by \code{sum.map.loc}.
#' @noRd
base.modules <- function(map.loc.agr, num.modules=1e4) {
  p.num <- NROW(map.loc.agr)
  num.modules <- min(num.modules, p.num - 1)
  cna.agr.diff.dist <- sort(abs(diff(map.loc.agr[,LogRatio])), decreasing=T)
  threshold <- cna.agr.diff.dist[num.modules]
  seg.starts <- vector('numeric')
  seg.ends <- vector('numeric')
  seg.chromosomes <- vector('character')
  chromosomes.levels <- levels(map.loc.agr[['Chromosome']])
  for (i in chromosomes.levels) {
    cna.agr <- map.loc.agr[Chromosome==i, LogRatio]
    chrom.start <- map.loc.agr[Chromosome==i, .(ChromosomeStart=.I),by=Chromosome][1,ChromosomeStart]
    ##########################
    # NOTE: magic number 1e9 #
    ##########################
    cna.agr.diff <- abs(diff(c(1e9, cna.agr))) > threshold
    curr.seg.starts <- which(cna.agr.diff)
    curr.seg.ends <- c(curr.seg.starts[2:length(curr.seg.starts)] - 1, length(cna.agr))
    curr.seg.starts <- curr.seg.starts + chrom.start - 1
    curr.seg.ends <- curr.seg.ends + chrom.start - 1
    seg.starts <- c(seg.starts, curr.seg.starts)
    seg.ends <- c(seg.ends, curr.seg.ends)
    seg.chromosomes <- c(seg.chromosomes, rep(i, length(curr.seg.starts)))
  }
  seg.chromosomes <- ordered(seg.chromosomes, levels=chromosomes.levels)
  mod.sort.i <- order(seg.starts)
  seg.starts <- seg.starts[mod.sort.i]
  seg.ends <- seg.ends[mod.sort.i]
  seg.chromosomes <- seg.chromosomes[mod.sort.i]

  modules <- Map(function(seg.start, seg.end) {
    list(I=seg.start,
         kw=seg.end - seg.start + 1,
         probes.ignore=vector('numeric'))
  }, seg.starts, seg.ends)
  # Ensure that extra factor levels are dropped
  list(modules=modules, chromosomes=seg.chromosomes[,drop=T])
}


#'
#' @noRd
###
# MATLAB MAPPING
# kws = params.kws
trunk.mods <- function (module1, module2, kws) {
  # Find best fit if (kw1 + kw2) is larger than the max kernel width allowed 
  mod1.probes <- module1$I:(module1$I + module1$kw - 1)
  # Strange try/catch structure here in the original MATLAB code
  # try
  mod2.probes <- module2$I:(module2$I + module2$kw - 1)
  # catch
  # mod2Probes = module2.I:(module2.I + module2.kw - 1)
  # end
  
  kw1 <- module1$kw
  kw2 <- module2$kw
  
  if (is.element('probes.ignore', names(module1))) {
    kw1 <- kw1 - length(module1$probes.ignore)
    mod1.probes <- setdiff(mod1.probes, module1$probes.ignore)
  }
  
  if (is.element('probes.ignore', names(module2))) {
    kw2 <- kw2 - length(module2$probes.ignore)
    mod2.probes <- setdiff(mod2.probes, module2$probes.ignore)
  }
  
  kw.all <- kw1 + kw2
  kw.max <- kws[length(kws)]
  if (kw.all <= kw.max) {
    return(list(mod1.probes=mod1.probes, mod2.probes=mod2.probes))
  }
  
  num.p.to.sub <- kw.all - kw.max
  kw.diff <- kw1 - kw2
  if (kw.diff >= num.p.to.sub) {
    kw1 <- kw1 - num.p.to.sub
  } else if (kw.diff > 0) {
    kw1 <- kw2
  } else if (kw.diff <= -num.p.to.sub) {
    kw2 <- kw2 - num.p.to.sub
  } else {
    kw2 <- kw1
  }
  kw.all <- kw1 + kw2
  if (kw.all > kw.max) {
    kw1 <- round(kw.max/2)
    kw2 <- kw.max - kw1
  }
  mod1.probes <- mod1.probes[(length(mod1.probes)-kw1+1):length(mod1.probes)]
  mod2.probes <- mod2.probes[1:kw2]
  
  list(mod1.probes=mod1.probes, mod2.probes=mod2.probes)
}


#'
#' @param kw1 Integer value for kw1.
#' @param kw2 Integer value for kw2.
#' @param kws The vector of kws as returned by \code{single.samp.param.estimation()$kws}.
#' @param mus The matrix of mus as returned by \code{all.samp.param.estimation()$mus}.
#' @param vars The matrix of vars as returned by \code{all.samp.param.estimation()$vars}.
#' @param int.sqrt.vs The matrix of int.sqrt.vs as returned by \code{all.samp.param.estimation()$int.sqrt.vs}.
#' @noRd
###
# MATLAB MAPPING
# kws = params.kws
# mus = params.mus
# vars = params.vars
# int.sqrt.vs = params.int.sqrt.vs
params.at.kws <- function(kw1, kw2, params) {
  kw.all <- kw1 + kw2
  
  kw.all.l.i <- tail(which(params$kws <= kw.all), n=1)
  kw.all.r.i <- which(params$kws >= kw.all)[1]
  kw.all.l <- params$kws[kw.all.l.i]
  kw.all.r <- params$kws[kw.all.r.i]
  
  # Left construction
  kw.break.l <- round(kw1/kw.all*kw.all.l)
  
  num.kws.half <- tail(which(params$kws <= kw.all.l/2), n=1)
  kw.breaks <- rep(0, 2 * num.kws.half)
  for (kw.j in seq_len(num.kws.half)) {
    kw.breaks[kw.j] <- params$kws[kw.j]
    kw.breaks[2*num.kws.half - kw.j + 1] <- kw.all.l - params$kws[kw.j]
  }
  
  unique.i <- which(!duplicated(kw.breaks))
  
  mus.l <- params$mus[kw.all.l.i, 1:(2 * num.kws.half)]
  vars.l <- params$vars[kw.all.l.i, 1:(2 * num.kws.half)]
  int.sqrt.vs.l <- params$int.sqrt.vs[kw.all.l.i, 1:(2 * num.kws.half)]
  
  if (length(unique.i) > 1) {
    kw.breaks.use <- kw.breaks[unique.i]
    mus.l.use <- mus.l[unique.i]
    vars.l.use <- vars.l[unique.i]
    int.sqrt.vs.l.use <- int.sqrt.vs.l[unique.i]
    break.i <- kw.breaks.use == kw.break.l
    
    if (sum(break.i) == 0) {
      mu.l <- interp1(kw.breaks.use, mus.l.use, kw.break.l, method='linear')
      var.l <- interp1(kw.breaks.use, vars.l.use, kw.break.l, method='linear')
      int.sqrt.v.l <- interp1(kw.breaks.use, int.sqrt.vs.l.use, kw.break.l, method='linear')
    } else {
      mu.l <- mus.l.use[break.i]
      var.l <- vars.l.use[break.i]
      int.sqrt.v.l <- int.sqrt.vs.l.use[break.i]
    }
  } else {
    mu.l <- mus.l[unique.i]
    var.l <- vars.l[unique.i]
    int.sqrt.v.l <- int.sqrt.vs.l[unique.i]
  }
  
  # Right construction
  kw.break.r <- round(kw1 / kw.all * kw.all.r)
  num.kws.half <- tail(which(params$kws <= kw.all.r/2), n=1)
  kw.breaks <- rep(0, 2 * num.kws.half)
  for (kw.j in seq_len(num.kws.half)) {
    kw.breaks[kw.j] <- params$kws[kw.j]
    kw.breaks[2 * num.kws.half - kw.j + 1] <- kw.all.r - params$kws[kw.j]
  }
  
  unique.i <- which(!duplicated(kw.breaks))
  
  mus.r <- params$mus[kw.all.r.i, 1:(2 * num.kws.half)]
  vars.r <- params$vars[kw.all.r.i, 1:(2 * num.kws.half)]
  int.sqrt.vs.r <- params$int.sqrt.vs[kw.all.r.i, 1:(2 * num.kws.half)]
  if (length(unique.i) > 1) {
    kw.breaks.use <- kw.breaks[unique.i]
    mus.r.use <- mus.r[unique.i]
    vars.r.use <- vars.r[unique.i]
    int.sqrt.vs.r.use <- int.sqrt.vs.r[unique.i]
    break.i <- kw.breaks.use == kw.break.r
    
    if (sum(break.i) == 0) {
      mu.r <- interp1(kw.breaks.use, mus.r.use, kw.break.r, method='linear')
      var.r <- interp1(kw.breaks.use, vars.r.use, kw.break.r, method='linear')
      int.sqrt.v.r <- interp1(kw.breaks.use, int.sqrt.vs.r.use, kw.break.r, method='linear')
    } else {
      mu.r <- mus.r.use[break.i]
      var.r <- vars.r.use[break.i]
      int.sqrt.v.r <- int.sqrt.vs.r.use[break.i]
    }
  } else {
    mu.r <- mus.r[unique.i]
    var.r <- vars.r[unique.i]
    int.sqrt.v.r <- int.sqrt.vs.r[unique.i]
  }
  
  a1 <- (kw.all - kw.all.l) / (kw.all.r - kw.all.l)
  a2 <- (kw.all.r - kw.all) / (kw.all.r - kw.all.l)
  
  if (kw.all.r != kw.all.l) {
    return(list(mu.est=a1*mu.r + a2*mu.l,
                var.est=a1*var.r + a2*var.l,
                int.sqrt.v.est=a1*int.sqrt.v.r + a2*int.sqrt.v.l))
  } else {
    return(list(mu.est=mu.l,
                var.est=var.l,
                int.sqrt.v.est=int.sqrt.v.l))
  }
}


#'
#' @noRd
###
# MATLAB MAPPING
# This is excatly same same function as getAnalNegLogEEuler()
anal.neg.log.e.euler <- function(int.sqrt.v, euler.num, t, num.tails) {
  log.b <- log(int.sqrt.v) - log(2*pi) - t^2/2
  log.a <- log(0.5*(euler.num)*erfc((t)/sqrt(2)))
  
  if (is.finite(log.a)) {
    e.euler <- -(log.a + log(1 + exp(log.b - log.a))) / log(10) -log10(num.tails)
  } else {
    e.euler <- -log.b/log(10) - log10(num.tails)
  }
  
  e.max <- -log10(int.sqrt.v / (2*pi) + 0.5*euler.num) - log10(num.tails)
  
  list(e.euler=e.euler, e.max=e.max)
}


#'
#' @param cna.matrix The matrix of log ratios as returned by \code{extract.matrix}.
#' @param kws The vector of kws as returned by \code{single.samp.param.estimation()$kws}.
#' @param mus The matrix of mus as returned by \code{all.samp.param.estimation()$mus}.
#' @param vars The matrix of vars as returned by \code{all.samp.param.estimation()$vars}.
#' @param int.sqrt.vs The matrix of int.sqrt.vs as returned by \code{all.samp.param.estimation()$int.sqrt.vs}.
#' @noRd
###
# MATLAB MAPPING
# kws = params.kws
# mus = params.mus
# vars = params.vars
# int.sqrt.vs = params.int.sqrt.vs
break.sig.v2 <- function(cna.agr, module1, module2, params) {
  
  trunk.results <- trunk.mods(module1, module2, params$kws)
  mod1.i <- trunk.results$mod1.probes
  mod2.i <- trunk.results$mod2.probes
  
  kw1 <- length(mod1.i)
  kw2 <- length(mod2.i)
  break.info <- list(reg.start=module1$I,
                     reg.end=module2$I + module2$kw - 1,
                     kw1=kw1,
                     kw2=kw2)
  
  params.results <- params.at.kws(kw1, kw2, params)
  dist.mu <- params.results$mu.est
  dist.var <- params.results$var.est
  int.sqrt.v <- params.results$int.sqrt.v.est
  
  f1 <- rowMeans(cna.agr[,mod1.i,drop=F])
  f2 <- rowMeans(cna.agr[,mod2.i,drop=F])
  f <- f2 - f1
  
  if (dist.var == 0) {
    break.info$t <- 0
    break.info$z <- 0
    ########################
    # NOTE: magic number p #
    ########################
    break.info$p <- -1e9
  } else {
    t <- sum(f)
    z <- (t - dist.mu) / sqrt(dist.var)
    if (abs(z) >= 0.5) {
      p <- anal.neg.log.e.euler(int.sqrt.v, 1, abs(z), 2)$e.euler
    } else {
      ########################
      # NOTE: magic number p #
      ########################
      p <- -1e9
    }
    break.info$t <- t
    break.info$z <- z
    break.info$p <- p
  }
  break.info
}


#'
#' @return The break.info list
#' @noRd
adj.sig <- function(cna.matrix, modules, check.i, params) {
  num.mods <- length(modules)
  if (check.i >= num.mods) {
    ########################
    # NOTE: magic number p #
    ########################
    # This is the break.info MATLAB structure
    return(list(kw1=1, kw2=1, t=0, z=0, p=-1e9, reg.start=0, reg.end=0))
  }
  break.sig.v2(cna.matrix, modules[[check.i]], modules[[check.i+1]], params)
}


adj.sigs <- function(cna.matrix, modules, params) {
  num.mods <- length(modules)
  break.info <- replicate(num.mods - 1, list())
  break.ps <- rep(0, num.mods - 1)
  for (mod.i in seq_len(num.mods - 1)) {
    break.info[[mod.i]] <- adj.sig(cna.matrix, modules, mod.i, params)
    break.ps[mod.i] <- break.info[[mod.i]]$p
  }
  list(break.info=break.info, break.ps=break.ps)
}


cna.cum.to.smooth <- function(cna.agr, p.num, kw1, kw2) {
  p.num <- min(p.num, ncol(cna.agr) - kw1 - kw2)
  lr.s.1 <- -(cna.agr[, (kw1 + 1):(p.num + kw1), drop=F] - cna.agr[,1:p.num, drop=F]) / kw1
  lr.s.2 <-  (cna.agr[, (kw2 + kw1 + 1):(p.num + kw2 + kw1), drop=F] - cna.agr[,(kw1 + 1):(p.num + kw1), drop=F]) / kw2
  cna.s <- lr.s.1 + lr.s.2
  cna.s
}


iter.at.t <- function(cna.agr, p.num, t, kw1, kw2) {
  t <- abs(t)
  cna.s <- cna.cum.to.smooth(cna.agr, p.num, kw1, kw2)
  num.perm <- nrow(cna.s)
  e.emp <- sum(rowSums(t(apply(cna.s >= t, 1, diff)) > 0)) + sum(rowSums(t(apply(cna.s <= -t, 1, diff)) > 0)) + sum(abs(cna.s[,1]) >= t)
  e.emp <- e.emp / num.perm
  e.emp
}


updated.adj.sigs <- function(cna.matrix, modules, merge.i, break.info, break.ps, params) {
  break.info[[merge.i]] <- NULL
  break.ps <- break.ps[-merge.i]
  if (merge.i > 1) {
    break.info[[merge.i-1]] <-  break.sig.v2(cna.matrix, modules[[merge.i-1]], modules[[merge.i]], params)
    break.ps[merge.i-1] <- break.info[[merge.i-1]]$p
  }
  if (merge.i < length(modules)) {
    break.info[[merge.i]] <- break.sig.v2(cna.matrix, modules[[merge.i]], modules[[merge.i+1]], params)
    break.ps[merge.i] <- break.info[[merge.i]]$p
  }
  list(break.info=break.info, break.ps=break.ps)
}


merge.adj.modules <- function(cna.matrix, modules, merge.i, break.info, break.ps, params) {
  module.l <- modules[[merge.i]]
  module.r <- modules[[merge.i+1]]
  merged.module <- module.l
  merged.module$kw <- module.l$kw + module.r$kw
  modules[[merge.i]] <- merged.module
  modules[[merge.i+1]] <- NULL
  c(list(modules=modules), updated.adj.sigs(cna.matrix, modules, merge.i, break.info, break.ps, params))
}


mod.sigs <- function(modules, break.info) {
  num.mods <- length(modules)
  modules[[1]]$l$amp <- 0
  ########################
  # NOTE: magic number p #
  ########################
  modules[[1]]$l$p <- -1e9
  modules[[length(modules)]]$r$amp <- 0
  ########################
  # NOTE: magic number p #
  ########################
  modules[[length(modules)]]$r$p <- -1e9
  modules[[length(modules)]]$amp <- 0
  ########################
  # NOTE: magic number p #
  ########################
  modules[[length(modules)]]$p  <- -1e9
  for (mod.i in seq_len(num.mods-1)) {
    modules[[mod.i]]$r$amp <- break.info[[mod.i]]$z
    modules[[mod.i]]$r$p <- break.info[[mod.i]]$p
    modules[[mod.i+1]]$l$amp <- break.info[[mod.i]]$z
    modules[[mod.i+1]]$l$p <- break.info[[mod.i]]$p
    modules[[mod.i]]$amp <- 0
    ########################
    # NOTE: magic number p #
    ########################
    modules[[mod.i]]$p <- -1e9
  }
  
  for (mod.i in seq_len(num.mods)) {
    if (modules[[mod.i]]$l$amp == 0) {
      modules[[mod.i]]$amp <- -modules[[mod.i]]$r$amp
      modules[[mod.i]]$p <- modules[[mod.i]]$r$p
    } else if (modules[[mod.i]]$r$amp == 0) {
      modules[[mod.i]]$amp <- modules[[mod.i]]$l$amp
      modules[[mod.i]]$p <- modules[[mod.i]]$l$p
    } else {
      modules[[mod.i]]$amp <- sign(modules[[mod.i]]$l$amp)*min(abs(modules[[mod.i]]$l$amp), abs(modules[[mod.i]]$r$amp))
      modules[[mod.i]]$p <- min(modules[[mod.i]]$l$p, modules[[mod.i]]$r$p)
    }
  }
  modules
}


#'
#' @param map.loc.agr The sum of the log ratio as returned by \code{sum.map.loc}.
#' @noRd
###
# MATLAB MAPPING
# p.num = data.P
# agr.cum.iter = params.agrCumIter
merge.modules.v2 <- function(cna.matrix, map.loc.agr, p.num, samp.is, modules,
                             sig.mods, agr.cum.iter, params, e.e,
                             use.perm=T, level=1) {
  if (isempty(samp.is)) {
    samp.is <- rep(T, length=nrow(cna.matrix))
  }
  
  if (!isempty(modules)) {
    modules.results <- sort.modules(modules, map.loc.agr[['Chromosome']])
  } else {
    modules.results <- base.modules(map.loc.agr)
  }
  
  modules <- modules.results$modules
  chromosomes <- modules.results$chromosomes
  unique.chromosomes <- levels(chromosomes)
  
  if (length(unique.chromosomes) > 1) {
    modules.use <- list()
    break.info <- list()
    for (chromosome in unique.chromosomes) {
      chromosome.i <- which(chromosomes == chromosome)
      merge.results <- merge.modules.v2(cna.matrix, map.loc.agr, p.num, samp.is,
                                        modules[chromosome.i], sig.mods, agr.cum.iter,
                                        params, e.e, level=level+1)
      curr.modules <- merge.results$modules
      curr.breaks <- merge.results$break.info
      sig.mods <- merge.results$sig.mods
      
      for (i in seq_len(length(curr.breaks))) {
        curr.breaks[[i]]$chromosome <- chromosome
      }
      modules.use <- c(modules.use, curr.modules)
      break.info <- c(break.info, curr.breaks)
    }
    modules <- modules.use
    return(list(modules=modules, break.info=break.info, sig.mods=sig.mods))
  }
  
  curr.cna.matrix <- t(colSums(cna.matrix[samp.is,]))
  
  adj.results <-  adj.sigs(curr.cna.matrix, modules, params)
  break.info <- adj.results$break.info
  break.ps <- adj.results$break.ps
  
  repeat {
    
    if (isempty(break.info))
      break
    
    min.p <- min(break.ps)
    ##########################
    # NOTE: magic number 1e9 #
    ##########################
    if (min.p == 1e9)
      break
    
    min.i <- which(break.ps == min.p)[1]
    
    if (min.p >= e.e) {
      
      if (use.perm) {
        do.perm <- F
        if (isempty(sig.mods)) {
          do.perm <- T
        } else {
          start.is <- (sig.mods$starts == break.info[[min.i]]$reg.start)
          end.is <- (sig.mods$ends == break.info[[min.i]]$reg.end)
          overlap.flag <- sum(start.is & end.is)
          if (overlap.flag == 0) {
            do.perm <- T
          }
        }
        if (do.perm) {
          min.p.perm <- iter.at.t(agr.cum.iter, p.num, break.info[[min.i]]$t, break.info[[min.i]]$kw1, break.info[[min.i]]$kw2)
          if (min.p.perm != 0) {
            min.p <- -log10(min.p.perm)
          }
          if (min.p >= e.e) {
            if (isempty(sig.mods)) {
              sig.mods$starts <- break.info[[min.i]]$reg.start
              sig.mods$ends <- break.info[[min.i]]$reg.end
              sig.mods$p <- min.p
            } else {
              sig.mods$starts <- c(sig.mods$starts, break.info[[min.i]]$reg.start)
              sig.mods$ends <- c(sig.mods$ends, break.info[[min.i]]$reg.end)
              sig.mods$p <- c(sig.mods$p, min.p)
            }
          }
        }
      }
      if (min.p >= e.e) {
        ##########################
        # NOTE: magic number 1e9 #
        ##########################
        break.ps[min.i] <- 1e9
        break.info[[min.i]]$p <- min.p
        next
      }
    }
    
    merge.results <- merge.adj.modules(curr.cna.matrix, modules, min.i, break.info, break.ps, params)
    modules <- merge.results$modules
    break.info <- merge.results$break.info
    break.ps <- merge.results$break.ps
  }
  
  modules <- mod.sigs(modules, break.info)
  
  if (level == 1) {
    for (i in seq_len(length(break.info))) {
      break.info[[i]]$chromosome <- 1
    }
  }
  
  return(list(modules=modules, break.info=break.info, sig.mods=sig.mods))
}