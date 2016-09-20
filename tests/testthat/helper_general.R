# Utility functions for the autotesters

# Convenience functions
as.r.per.sample <- function(mat.per.sample, nr, nc) {
  r.per.sample <- matrix(mat.per.sample, nrow=nr, ncol=nc)
  for(i in 1:nr) {
    for (j in 1:nc) {
      el <- r.per.sample[i,j]
      if (!is.null(el[[1]])) {
        r.el <- list()
        r.el[[1]] <- list(rapply(t(el[[1]][[1]]["mus",,]),c),
                          rapply(t(el[[1]][[1]]["vars",,]),c),
                          rapply(t(el[[1]][[1]]["diffVars",,]),c),
                          rapply(t(el[[1]][[1]]["covs",,]),c))
        names(r.el[[1]]) <- c("mus", "vars", "diff.vars", "covs")
        r.per.sample[i,j] <- r.el
      }
    }
  }
  r.per.sample
}


extract.chromosome.levels.test <- function(chromosomes) {
  chromosome.names <- toupper(unique(chromosomes))
  chromosome.numbers <- suppressWarnings(as.numeric(chromosome.names))
  chromosome.symbols <- chromosome.names[is.na(chromosome.numbers)]
  c(sort(chromosome.numbers),sort(chromosome.symbols))
}


# Convenience functions
as.r.chr <- function(mat.chr, levels=FALSE) {
  r.chr <- c()
  for (el in mat.chr) {
    if (el < 23)
      r.chr <- c(r.chr,as.character(el))
    else if (el == 23)
      r.chr <- c(r.chr,"X")
    else if (el == 24)
      r.chr <- c(r.chr,"Y")
  }
  if (levels == TRUE) {
    r.chr <- ordered(r.chr, levels=extract.chromosome.levels.test(r.chr))
  }
  r.chr
}


read.rubic.maploc <- function(mat.in, chr.levels=TRUE, probe.levels=TRUE) {
  mat.data <- mat.in$data[,,1]
  mat.numProbes <- as.numeric(mat.data$P)
  mat.numSamples <-  as.numeric(mat.data$S)
  
  map.loc <- data.table(
    Sample=rep(rapply(mat.data$sampleName,c), each=mat.numProbes),
    Chromosome=rep(as.r.chr(as.vector(mat.data$chrom),levels=FALSE),mat.numSamples),
    Probe=rep(rapply(mat.data$probes,c),mat.numSamples),
    Position=rep(as.vector(mat.data$relMapLoc),mat.numSamples),
    AbsPosition=rep(as.vector(mat.data$mapLoc),mat.numSamples),
    LogRatio=c(t(mat.data$cna)))
  
  if (chr.levels == TRUE) {
    map.loc[['Chromosome']] <- ordered(map.loc[['Chromosome']],extract.chromosome.levels.test(map.loc[['Chromosome']]))
  }
  if (probe.levels == TRUE) {
    map.loc[['Probe']] <- ordered(map.loc[['Probe']], levels=map.loc[['Probe']][1:mat.numProbes])
  }  
  map.loc
  
}


read.rubic.geneinfo <- function(mat.in) {
  mat.gene.info <- mat.in$geneInfo[,,1]
  gene.info <- data.table(ID=NA,
                          Name=as.vector(unlist(mat.gene.info$symbs)),
                          Start=as.vector(mat.gene.info$starts),
                          End=as.vector(mat.gene.info$ends),
                          Chromosome=as.r.chr(as.vector(mat.gene.info$chroms), levels=T))
  setkey(gene.info, Name)
  gene.info
}


as.r.module <- function(mat.module) {
  r.module <- list()
  r.module[["I"]] <- unlist(mat.module["I",,],use.names=FALSE)
  r.module[["kw"]] <- unlist(mat.module["kw",,],use.names=FALSE)
  r.module[["probes.ignore"]] <- unlist(mat.module["probesIgnore",,],use.names=FALSE)
  r.module
}


as.r.module.long <- function(mat.module) {
  r.module <- list()
  r.module[["I"]] <- unlist(mat.module[[1]]["I",,],use.names=FALSE)
  r.module[["kw"]] <- unlist(mat.module[[1]]["kw",,],use.names=FALSE)
  r.module[["probes.ignore"]] <- unlist(mat.module[[1]]["probesIgnore",,],use.names=FALSE)
  r.module[["l"]] <- list()
  r.module[["l"]][["amp"]] <- unlist(mat.module[[1]]["L",,]$L["amp",,],use.names=FALSE)
  r.module[["l"]][["p"]] <- unlist(mat.module[[1]]["L",,]$L["p",,],use.names=FALSE)  
  r.module[["r"]] <- list()
  r.module[["r"]][["amp"]] <- unlist(mat.module[[1]]["R",,]$R["amp",,],use.names=FALSE)
  r.module[["r"]][["p"]] <- unlist(mat.module[[1]]["R",,]$R["p",,],use.names=FALSE)  
  r.module[["amp"]] <- unlist(mat.module[[1]]["amp",,],use.names=FALSE)
  r.module[["p"]] <- unlist(mat.module[[1]]["p",,],use.names=FALSE)
  r.module
}


as.r.focal.event <- function(mat.module) {
  r.module <- as.r.module.long(mat.module)
  r.module[["chromosome"]] <- unlist(mat.module[[1]]["chromStr",,],use.names=FALSE)
  r.module[["loc.start"]] <- unlist(mat.module[[1]]["locStart",,],use.names=FALSE)
  r.module[["loc.end"]] <- unlist(mat.module[[1]]["locEnd",,],use.names=FALSE)
  r.module
}


as.r.focal.event.long <- function(mat.module) {
  r.module <- as.r.focal.event(mat.module)
  r.module[["gene.symbols"]] <- sort(unlist(mat.module[[1]]["geneSymbs",,],use.names=FALSE))
  r.module
}


as.r.focal.event.final <- function(mat.module) {
  r.module <- as.r.focal.event.long(mat.module)
  r.module$r$q <- unlist(mat.module[[1]]["R",,]$R["q",,], use.names=FALSE)
  r.module$l$q <- unlist(mat.module[[1]]["L",,]$L["q",,], use.names=FALSE)
  r.module
}

as.r.break.info.chrom <- function(mat.break.info) {
  r.bi <- lapply(mat.break.info, 
                 function(x) setNames(lapply(x[[1]][,,1], as.vector), c('reg.start', 'reg.end', 'kw1', 'kw2', 't', 'z', 'p','chromosome')))
  for (i in 1:length(r.bi)) {
    r.bi[[i]]$chromosome <- as.r.chr(r.bi[[i]]$chromosome, levels=FALSE)
  }
  r.bi
}


as.r.break.info <- function(mat.break.info) {
  lapply(mat.break.info, 
         function(x) setNames(lapply(x[[1]][,,1], as.vector), c('reg.start', 'reg.end', 'kw1', 'kw2', 't', 'z', 'p')))
}