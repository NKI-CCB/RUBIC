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


min.kw.peak.index <- function(segments, t.sign, min.kw.peak=1e10) {
  num.segs <- length(segments)
  if (!is.element('l', names(segments[[1]]))) {
    return(NULL)
  }
  min.kw.peak.i <- 0
  
  for (i in 1:num.segs) {
    peak <- FALSE
    if (!is.element('recurr.peak', names(segments[[i]]))) {
      if (((segments[[i]]$l$amp == 0) && (segments[[i]]$r$amp == 0)) ||
          ((segments[[i]]$l$amp == 0) && (sign(segments[[i]]$r$amp) == -t.sign)) ||
          ((segments[[i]]$r$amp == 0) && (sign(segments[[i]]$l$amp) == t.sign)) ||
          ((sign(segments[[i]]$l$amp) == t.sign) && ((sign(segments[[i]]$l$amp) !=  sign(segments[[i]]$r$amp)))))
            {
              peak <- TRUE
            }
      if (peak && (segments[[i]]$kw < min.kw.peak)) {
        min.kw.peak.i <- i
        min.kw.peak <- segments[[i]]$kw
      }
    }
  }
  if (min.kw.peak.i == 0) {
    NULL
  } else {
    min.kw.peak.i
  }
}

percentile.v2 <- function(t, p.num, agr.cum.iter, kw, t.sign) {
  cna.s <- (agr.cum.iter[, (kw+1):(p.num+kw), drop=F] - agr.cum.iter[, 1:p.num, drop=F]) / kw
  mu.s <- mean(rowMeans(cna.s))
  var.s <- mean(apply(cna.s, 1, var))
  z <- (t - mu.s) / sqrt(var.s)
  
  if (t.sign == +1) {
    cna.b <- cna.s >= t
    z <- max(z, 0)
  }
  else {
    cna.b <- cna.s <= t
    z <- min(z, 0)
  }
  perc <- mean(cna.b)
  list(perc=perc, z=z)
}

is.peak.valid <- function(p.num, agr.cum.iter, segments, min.kw.peak.i, cna.agr, fdr, t.sign) {
  sig.info <- list(amp=0, p=-1e9)
  peak.valid <- FALSE
  
  i <- min.kw.peak.i
  
  seg.is <- segments[[i]]$I:(segments[[i]]$I + segments[[i]]$kw - 1)

  if (is.element('probes.ignore', names(segments[[i]]))) {
    seg.is <- setdiff(seg.is, segments[[i]]$probes.ignore)
  }
  
  seg.avg <- mean(cna.agr[seg.is])
  
  pz <- percentile.v2(seg.avg, p.num, agr.cum.iter, segments[[i]]$kw, t.sign)
  perc <- pz$perc
  z <- pz$z
  
  if (perc <= fdr) {
    sig.info$amp <- z
    sig.info$p <- -log10(perc)
    peak.valid <- TRUE
  } 
  
  list(peak.valid=peak.valid, sig.info=sig.info)
}

broadest.adj <- function(focal.events, seg.loc, type) {
  ind <- 0
  max.kw <- 0
  for (i in seq_len(length(focal.events))) {
    if (type == -1){
      curr.loc <- focal.events[[i]]$I + focal.events[[i]]$kw
    } else {
      curr.loc <- focal.events[[i]]$I - 1
    }
    if (curr.loc == seg.loc) {
      if (focal.events[[i]]$kw > max.kw) {
        max.kw <- focal.events[[i]]$kw
        ind <- i
      }
    }
  }
  
  ind
}

adjust.boundaries <- function(focal.events, segment) {
  seg.start <- segment$I
  seg.end <- segment$I + segment$kw - 1
  seg.sign <- sign(segment$amp)
  
  focal.l.I <- broadest.adj(focal.events, seg.start, -1)
  focal.r.I <- broadest.adj(focal.events, seg.end, +1)
  if (focal.l.I != 0) {
    if ((sign(focal.events[[focal.l.I]]$amp) == seg.sign) && (focal.events[[focal.l.I]]$r$amp != 0)) {
      segment$I <- focal.events[[focal.l.I]]$I
      segment$kw <- segment$kw + focal.events[[focal.l.I]]$kw
      focal.start <- focal.events[[focal.l.I]]$I
      focal.end <- focal.events[[focal.l.I]]$I + focal.events[[focal.l.I]]$kw - 1
      segment$probes.ignore <- unique(c(focal.start:focal.end, segment$probes.ignore))
    }
  }
  
  if (focal.r.I != 0) {
    if ((sign(focal.events[[focal.r.I]]$amp) == seg.sign) && (focal.events[[focal.r.I]]$l$amp != 0)) {
      segment$kw <- segment$kw + focal.events[[focal.r.I]]$kw
      focal.start <- focal.events[[focal.r.I]]$I
      focal.end <- focal.events[[focal.r.I]]$I + focal.events[[focal.r.I]]$kw - 1
      segment$probes.ignore <- unique(c(segment$probes.ignore, focal.start:focal.end))
    }
  }
  
  segment
}


join.adj.segs <- function(cna.matrix, segments, rel.peak.i, params, e) {
  
  if ((segments[[rel.peak.i]]$l$amp == 0) && (segments[[rel.peak.i]]$r$amp == 0)) {
    segments[rel.peak.i] <- NULL
    return(segments)
  }
  
  if (segments[[rel.peak.i]]$l$amp == 0) {
    segments[[rel.peak.i+1]]$l$amp <- 0
    segments[[rel.peak.i+1]]$l$p <- -1e9
    #if (!isempty(segments[[rel.peak.i+1]]$recurr.peak)) {
    if (is.element('recurr.peak', names(segments[[rel.peak.i+1]]))) {
      segments[[rel.peak.i+1]]$recurr.peak <- NULL
    }
    segments[rel.peak.i] <- NULL
    return(segments)
  }
  
  if (segments[[rel.peak.i]]$r$amp == 0) {
    segments[[rel.peak.i-1]]$r$amp <- 0
    segments[[rel.peak.i-1]]$r$p <- -1e9
    if (is.element('recurr.peak', names(segments[[rel.peak.i-1]]))) {
      segments[[rel.peak.i-1]]$recurr.peak <- NULL
    }
    segments[rel.peak.i] <- NULL
    return(segments)
  }

  seg.l <- segments[[rel.peak.i-1]]
  seg.r <- segments[[rel.peak.i+1]]
  break.info <- break.sig.v2(cna.matrix, seg.l, seg.r, params)
  
  segments[[rel.peak.i-1]]$r$amp <- break.info$z
  segments[[rel.peak.i-1]]$r$p <- break.info$p
  segments[[rel.peak.i+1]]$l$amp <- break.info$z
  segments[[rel.peak.i+1]]$l$p <- break.info$p
  
  if (break.info$p >= -log10(e)) {
    e.perm <- iter.at.t(params$agr.cum.iter, ncol(cna.matrix), break.info$t, break.info$kw1, break.info$kw2)
    if (e.perm == 0) {
      p.perm <- break.info$p
    } else {
      p.perm <- -log10(e.perm)
    }
    if (p.perm >= -log10(e)) {
      segments[[rel.peak.i-1]]$r$p <- p.perm
      segments[[rel.peak.i+1]]$l$p <- p.perm
      if (!isempty(segments[[rel.peak.i-1]]$recurr.peak)) {
        segments[[rel.peak.i-1]]$recurr.peak <- NULL
      }
      if (!isempty(segments[[rel.peak.i+1]]$recurr.peak)) {
        segments[[rel.peak.i+1]]$recurr.peak <- NULL
      }
      segments[rel.peak.i] <- NULL
      return(segments)
    }
  }  
  probes.ignore <- (segments[[rel.peak.i-1]]$I + segments[[rel.peak.i-1]]$kw):(segments[[rel.peak.i+1]]$I - 1)
  probes.ignore <- unique(c(segments[[rel.peak.i-1]]$probes.ignore, probes.ignore))
  probes.ignore <- unique(c(probes.ignore, segments[[rel.peak.i+1]]$probes.ignore))
  
  probe.start <- segments[[rel.peak.i-1]]$I
  probe.end <- segments[[rel.peak.i+1]]$I + segments[[rel.peak.i+1]]$kw - 1
  segments[[rel.peak.i-1]]$I <- probe.start
  segments[[rel.peak.i-1]]$kw <- probe.end - probe.start + 1
  segments[[rel.peak.i-1]]$r <- segments[[rel.peak.i+1]]$r
  segments[[rel.peak.i-1]]$probes.ignore <- probes.ignore
  if (is.element('recurr.peak', names(segments[[rel.peak.i-1]]))) { # why are we doing this?
    segments[[rel.peak.i-1]]$recurr.peak <- NULL
  }
  
  segments[rel.peak.i:(rel.peak.i+1)] <- NULL
  
  segments
}


call.events <- function(cna.matrix, amp.level, del.level,
                        segments, params, fdr, e, t.sign, samps.use=NULL) {
  p.num <- ncol(cna.matrix)
  
  if (isempty(samps.use)) {
    samps.use <- rep(TRUE, nrow(cna.matrix))
  }
  
  cna.p <- cna.matrix[samps.use,, drop=F]
  cna.p[cna.p < amp.level] <- 0
  cna.p.agr <- colSums(cna.p)
  cna.n <- cna.matrix[samps.use,, drop=F ]
  cna.n[cna.n > del.level] <- 0
  cna.n.agr <- colSums(cna.n)
  
  focal.events <- list()
  
  
  while (TRUE) {
    min.kw.peak.i <- min.kw.peak.index(segments, t.sign)
    
    if (is.null(min.kw.peak.i)) {
      break
    }
    
    if (t.sign == +1 ){
      peak.valid.results <- is.peak.valid(p.num, params$agr.cum.iter, segments, min.kw.peak.i, cna.p.agr, fdr, t.sign)
    } else {
      peak.valid.results <- is.peak.valid(p.num, params$agr.cum.iter, segments, min.kw.peak.i, cna.n.agr, fdr, t.sign)
    }
    peak.valid <- peak.valid.results$peak.valid
    sig.info <- peak.valid.results$sig.info
    
    if (peak.valid){
      segments[[min.kw.peak.i]]$amp <- sig.info$amp 
      segments[[min.kw.peak.i]]$p <- sig.info$p
      segments[[min.kw.peak.i]] <- adjust.boundaries(focal.events, segments[[min.kw.peak.i]])
      if (isempty(focal.events)) {
        focal.events <- list(segments[[min.kw.peak.i]])
      } else {
        focal.events <- c(focal.events, list(segments[[min.kw.peak.i]]))
      }
      if (t.sign == 1) {
        segments <- join.adj.segs(cna.p, segments, min.kw.peak.i, params, e)
      } else {
        segments <- join.adj.segs(cna.n, segments, min.kw.peak.i, params, e)
      }
      
    } else {
      segments[[min.kw.peak.i]]$recurr.peak <- FALSE
    }
  }
  focal.events
}