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


#' Preprocess the gene.info \code{data.table} for further manipulations.
#' 
#' @param genes.info The \code{data.table} containing the gene informations.
#' @param filter.chr The \code{vector} of chromosomes to select. If NULL (default)
#'                   all chromosomes present in the file will be used.
#' @noRd
preprocess.genes.info <- function(genes.info, filter.chr=NULL) {
  expected.columns <- c('Ensembl Gene ID', 'Associated Gene Name',
                        'Chromosome Name', 'Gene Start (bp)', 'Gene End (bp)')
  column.names <- colnames(genes.info)
  if (!all(expected.columns %in% column.names)) {
    stop('Unexpected column names. Wrong file format?')
  }
  short.names <- c('ID', 'Name', 'Chromosome', 'Start', 'End')
  column.names <- colnames(genes.info)
  setnames(genes.info, expected.columns, short.names)
  # If specified only the chromosomes in filter.char will be retained
  if (!isempty(filter.chr)) {
    setkey(genes.info, Chromosome)
    genes.info <- genes.info[J(filter.chr)]
  }
  setkey(genes.info, Name)
  chromosome.levels <- extract.chromosome.levels(genes.info)
  genes.info[,Chromosome:=ordered(toupper(Chromosome), chromosome.levels)]
  genes.info
}


#' Read Biomart annotation file.
#'  
#' @param build.file The TSV file exported from Biomart.
#' @param filter.chr A character vector of wanted chromosomes.
#' 
#' @return The \code{data.table} containing the genome annotations.
#' 
#' @section Warning:
#' This function is expecting at least the following headers: 'Ensembl Gene ID',
#' 'Gene Start (bp)', 'Gene End (bp)', 'Chromosome Name', 'Associated Gene Name'.
#' 
#' @export
read.genes.info.tsv <- function(build.file, filter.chr=NULL) {
  #####################################################################
  # NOTE: there is a bug in the original code. Gene start positions   #
  #       are used for both start and end. The following code is      #
  #       generating the error:                                       #
  #       geneStarts = str2double(rawMatrix(2:end, 3));               #
  #       geneEnd = str2double(rawMatrix(2:end, 3));                  #
  #       -- in getGeneInfoFromTxt(), lines 7 and 8.                  #
  #####################################################################
  if (file.access(build.file, mode=4) == -1)
    stop(paste('The genes locations file', build.file, 'cannot be found or read'))
  
  preprocess.genes.info(fread(build.file), filter.chr)
}


#' Import Biomart annotations from the web service.
#' 
#' @param filter.chr A character vector of wanted chromosomes.
#' 
#' @return The \code{data.table} containing the genome annotations.
#' @export
read.genes.info.biomart <- function(filter.chr=NULL) {
  base.url <- 'http://www.biomart.org/biomart/martservice?query='
  query.xml <- c('<!DOCTYPE Query><Query client="true" processor="TSV" limit="-1" header="1">',
                 '<Dataset name="hsapiens_gene_ensembl" config="gene_ensembl_config">')
  if (!isempty(filter.chr)) {
    query.xml <- c(query.xml, paste0('<Filter name="chromosome_name" value="',
                                     paste(filter.chr, collapse=','),
                                     '"/>'))
  }
  query.xml <- c(query.xml, '<Attribute name="ensembl_gene_id"/>',
                 '<Attribute name="external_gene_name"/>',
                 '<Attribute name="chromosome_name"/>',
                 '<Attribute name="start_position"/>',
                 '<Attribute name="end_position"/>',
                 '</Dataset></Query>')

  url <- paste0(base.url, URLencode(paste0(query.xml, collapse='')))
  
  # From data.table > 1.9.5 fread passes showProgress=FALSE through
  # to download.file() as quiet=!showProgress. Thanks to a pull request from
  # Karl Broman and Richard Scriven for filing the issue, #741
  
  # For the moment this function display a very annoying red message.
  preprocess.genes.info(suppressWarnings(fread(url)))
}


filter.overlaps <- function(focal.events) {
  if (isempty(focal.events)) {
    return(focal.events)
  }
  
  num.events <- length(focal.events)
  focal.i <- rep(T, num.events)
  starts <- rep(0, num.events)
  ends <- rep(0, num.events)
  chroms <- rep(0, num.events)

  for (i in seq_len(num.events)) {
    starts[i] <- focal.events[[i]]$loc.start
    ends[i] <- focal.events[[i]]$loc.end
    chroms[i] <- focal.events[[i]]$chromosome
  }
  
  lengths <- ends - starts
  len.i <- order(lengths)
  for (i in seq_len(num.events)) {
    if (focal.i[len.i[i]] == F) {
      next
    }
    l.i <- starts <= starts[len.i[i]]
    r.i <- ends >= ends[len.i[i]]
    chrom.i <- chroms == chroms[len.i[i]]
    del.i <- setdiff(which(l.i & r.i & chrom.i), len.i[i])
    focal.i[del.i] <- F
  }
  focal.events[focal.i]
}


gene.locs <- function(map.loc.agr, called.events, markers=NULL) {
  lapply(called.events, function(event) {
    locs <- map.loc.agr[c(event$I, event$I + event$kw - 1), .(Chromosome, Position)]
    loc.start <- locs[1, Position]
    loc.end <- locs[2, Position]
    # Assert that locs[1, Chromosome] == locs[2, Chromosome]
    loc.chr <- as.character(locs[1, Chromosome])
    if (!is.null(markers)) {
      # If there are not marker found by the query a warning is displayed and
      # -/+ Inf is returned
      loc.adj.start <- suppressWarnings(markers[Chromosome==loc.chr & Position < loc.start,
                                               max(Position)])
      loc.adj.end <- suppressWarnings(markers[Chromosome==loc.chr & Position > loc.end,
                                             min(Position)])
      if (is.infinite(loc.adj.start)) {
        loc.adj.start <- loc.start
      }
      if (is.infinite(loc.adj.end)) {
        loc.adj.end <- loc.end
      }
    } else {
      loc.adj.start <- loc.start
      loc.adj.end <- loc.end
    }
    event$chromosome <- loc.chr
    event$loc.start <- loc.adj.start
    event$loc.end <- loc.adj.end
    event
  })
}


insert.genes <- function(genes.info, focal.events) {
  lapply(focal.events, function(event) {
    names <- genes.info[Chromosome==event$chromosome & Start < event$loc.end & End > event$loc.start, Name]
    event$gene.symbols <- sort(unique(names))
    event
  })
}


called.to.genes <- function(map.loc.agr, called.events, max.length, genes.info, markers=NULL) {

  focal.events <- gene.locs(map.loc.agr, called.events, markers)
  
  # Filter focal events
  if (max.length != 0) {
    focal.events <- Filter(function(event) {
      if ((event$loc.end - event$loc.start) < max.length) {
        return(T)
      }
      F
    }, focal.events)
  } else {
    focal.events <- called.events
  }
  
  # Filter geneless events
  focal.events <- Filter(function(event) {
    if (NROW(genes.info[Chromosome==event$chromosome & Start < event$loc.end & End > event$loc.start]) > 0) {
      return(T)
    }
    F
  }, focal.events)
  
  # Filter overlaps
  focal.events <- filter.overlaps(focal.events)
  
  # Inser genes
  focal.events <- insert.genes(genes.info, focal.events)

  focal.events
}


calc.break.qvalues <- function(focal.p.events,  focal.n.events) {
  len.p.events <- length(focal.p.events)
  len.n.events <- length(focal.n.events)
  p.l.p <- rep(0, len.p.events)
  p.r.p <- rep(0, len.p.events)
  p.l.n <- rep(0, len.n.events)
  p.r.n <- rep(0, len.n.events)

  for (i in seq_len(len.p.events)) {
    p.l.p[i] <- focal.p.events[[i]]$l$p
    p.r.p[i] <- focal.p.events[[i]]$r$p
    focal.p.events[[i]]$l$q <- focal.p.events[[i]]$l$p
    focal.p.events[[i]]$r$q <- focal.p.events[[i]]$r$p
  }

  for (i in seq_len(len.n.events)) {
    p.l.n[i] <- focal.n.events[[i]]$l$p
    p.r.n[i] <- focal.n.events[[i]]$r$p
    focal.n.events[[i]]$l$q <- focal.n.events[[i]]$l$p
    focal.n.events[[i]]$r$q <- focal.n.events[[i]]$r$p
  }

  num.tests <- 2 * (len.p.events + len.n.events)
  q.all <- rep(-1e9, num.tests)
  for (test.i in seq_len(num.tests)) {
    max.l.p <- max(p.l.p)
    max.r.p <- max(p.r.p)
    max.l.n <- max(p.l.n)
    max.r.n <- max(p.r.n)
    max.all <- max(max.l.p, max.r.p, max.l.n, max.r.n)
    if (max.all == -1e9) {
      break
    }
    if (max.l.p == max.all) {
      max.i <- which(p.l.p == max.l.p)[1]
      focal.p.events[[max.i]]$l$q <- max.l.p + log10(test.i)
      p.l.p[max.i] <- -1e9
    } else if (max.r.p == max.all) {
      max.i <- which(p.r.p == max.r.p)[1]
      focal.p.events[[max.i]]$r$q <- max.r.p + log10(test.i)
      p.r.p[max.i] <- -1e9
    } else if (max.l.n == max.all) {
      max.i <- which(p.l.n == max.l.n)[1]
      focal.n.events[[max.i]]$l$q <- max.l.n + log10(test.i)
      p.l.n[max.i] <- -1e9
    } else if (max.r.n == max.all) {
      max.i <- which(p.r.n == max.r.n)[1]
      focal.n.events[[max.i]]$r$q <- max.r.n + log10(test.i)
      p.r.n[max.i] <- -1e9
    }
    q.all[test.i] <- max.all + log10(test.i)
  }
  
  q.all <- sort(q.all[q.all != -1e9])
  list(focal.p.events=focal.p.events,  focal.n.events=focal.n.events,
       q.all=q.all)
}


sort.regions.on.genome <- function(focal.events) {
  focal.events[order(sapply(focal.events, function(x) x$I))]
}