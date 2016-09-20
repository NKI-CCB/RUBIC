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


# NOTE: never called? not tested.
num.broad.peaks <- function(cna.matrix, amp.level, del.level, modules, params, fdr, e) {
  p.num <- ncol(cna.matrix)
  if (max(max(cna.matrix) > amp.level)) {
    t.sign <- +1
  } else {
    t.sign <- -1
  }
  all.events <- call.events(cna.matrix, amp.level, del.level, params$agr.cum.iter,
                            modules, params, fdr, e, t.sign)
  called.i <- rep(F, p.num)
  for (i in seq_len(length(all.events))) {
    called.i[all.events[[i]]$I:(all.events[[i]]$I + all.events[[i]]$kw - 1)] <- T
  }
  num.peaks <- sum(diff(c(F, called.i)) == 1)
}

compute.mod.is <- function (modules, sig.mods, sig.i) {
  mod.l.i <- 0
  mod.r.i <- 0
  for (i in seq_along(modules)) {
    if (modules[[i]]$I == sig.mods$starts[sig.i]) {
      mod.l.i <- i
    }
    if ((modules[[i]]$I + modules[[i]]$kw - 1) == sig.mods$ends[sig.i]) {
      mod.r.i <- i
    }
  }
  
  list(mod.l.i=mod.l.i, mod.r.i=mod.r.i)
}

analytical.to.perm.p <- function(modules, break.info, sig.mods) {
  num.breaks <- length(break.info)
  for (break.i in seq_len(num.breaks)) {
    start.is <- sig.mods$starts == break.info[[break.i]]$reg.start
    end.is <- sig.mods$ends == break.info[[break.i]]$reg.end
    sig.i <- which(start.is & end.is)[1]
    if (!is.na(sig.i)) {
      p <- sig.mods$p[sig.i]
      break.info[[break.i]]$p <- p
      mod.is <- compute.mod.is(modules, sig.mods, sig.i)
      modules[[mod.is$mod.l.i]]$r$p <- p
      modules[[mod.is$mod.r.i]]$l$p <- p
    }
  }
  
  num.mods <- length(modules)
  for (mod.i in seq_len(num.mods)) {
    if (modules[[mod.i]]$l$amp == 0) {
      modules[[mod.i]]$p <- modules[[mod.i]]$r$p
    } else if (modules[[mod.i]]$r$amp == 0) {
      modules[[mod.i]]$p <- modules[[mod.i]]$l$p
    } else {
      modules[[mod.i]]$p <- min(modules[[mod.i]]$l$p, modules[[mod.i]]$r$p)
    }
  }
  
  list(modules=modules, break.info=break.info)
}

segments.at.fdr.abs <- function(cna.matrix, map.loc.agr, amp.level, del.level, samps.use, params, fdr, t.sign, broad.fdr=F) {
  p.num <-ncol(cna.matrix)
  
  if (t.sign == 1) {
    cna.matrix[cna.matrix < amp.level] <- 0
  } else {
    cna.matrix[cna.matrix > del.level] <- 0
  }
  e <- fdr
  #########################
  # NOTE: magic number -3 #
  #########################
  merge.results <- merge.modules.v2(cna.matrix, map.loc.agr, p.num,
                                    samps.use, matrix(nrow=0, ncol=0), list(),
                                    params$agr.cum.iter, params, -3, use.perm=F)
  modules.base <- merge.results$modules
  
  modules.prev <- list()
  break.info.prev <- list()
  sig.mods.prev <- list()
  num.peaks.prev <- 0
  
  repeat {
    merge.results <- merge.modules.v2(cna.matrix, map.loc.agr, p.num,
                                      samps.use, modules.base, sig.mods.prev,
                                      params$agr.cum.iter, params, -log10(e))
    modules <- merge.results$modules
    break.info <- merge.results$break.info
    sig.mods <- merge.results$sig.mods
    
    if (broad.fdr) {
      num.peaks <- num.broad.peaks(cna.matrix, amp.level, del.level, modules, params, fdr, e)
    } else {
      num.peaks <- 2 * sum(diff(sapply(break.info, function(x) x$z) > 0) > 0) + 1
    }
    if (num.peaks <= num.peaks.prev) {
      break
    }
    modules.prev <- modules
    break.info.prev <- break.info
    sig.mods.prev <- sig.mods
    num.peaks.prev <- num.peaks
    e <- num.peaks * fdr
  }
  e <- num.peaks.prev * fdr
  modules <- modules.prev
  break.info <- break.info.prev
  
  c(analytical.to.perm.p(modules, break.info, sig.mods), e=e)
}

aggregate.segments <- function(cna.matrix, map.loc.agr, amp.level, del.level, params.p, params.n, fdr, samps.use=vector('numeric')) {

  segments.p <- segments.at.fdr.abs(cna.matrix, map.loc.agr, amp.level, del.level, samps.use, params.p, fdr, +1)
  segments.n <- segments.at.fdr.abs(cna.matrix, map.loc.agr, amp.level, del.level, samps.use, params.n, fdr, -1)
  
  list(segments.p=segments.p$modules, segments.n=segments.n$modules, e.p=segments.p$e, e.n=segments.n$e)
}