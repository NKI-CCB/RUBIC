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

extract.chromosome.levels <- function(data) {
  chromosome.names <- toupper(unique(data[,Chromosome]))
  chromosome.numbers <- suppressWarnings(as.numeric(chromosome.names))
  chromosome.symbols <- chromosome.names[is.na(chromosome.numbers)]
  c(sort(chromosome.numbers),sort(chromosome.symbols))
}

#' Read the segmented CNA file in TSV format.
#' 
#' @param seg.file The name of the file containig the segmented CNA.
#' @param sample The index of the sample name column in the \code{seg.file}.
#' @param chromosome The index of the chromosome column in the \code{seg.file}.
#' @param start The index of the segment start column in the \code{seg.file}.
#' @param end The index of the segment end column in the \code{seg.file}.
#' @param log.ratio The index of the segment log ratio column in the \code{seg.file}.
#' @param header A boolean value that specifies whether the file contain a header line.
#' 
#' @return A \code{data.table} containing the segmented CNA.
#' @export
read.seg.file <- function(seg.file, sample=1, chromosome=2, start=3, end=4, log.ratio=6, header=T) {
  if (file.access(seg.file, mode=4) == -1)
    stop(paste('The segments file', seg.file, 'cannot be found or read'))
  # Order the columns based on their indices
  columns <- sort(c(sample=sample, chromosome=chromosome, start=start, end=end, log.ratio=log.ratio))
  seg.CNA <- fread(seg.file, sep='\t', select=columns, header=header)
  # Associate the new index to the columns
  columns <- setNames(1:5, names(columns))
  # Reorder the columns
  setcolorder(seg.CNA, c(columns['sample'], columns['chromosome'],
                         columns['start'], columns['end'], columns['log.ratio']))
  # Rename the columns
  setnames(seg.CNA, c('Sample', 'Chromosome', 'Start', 'End', 'LogRatio'))
  
  # Make sure chromosomes stay ordered independently of their alphabetical order
  chromosome.levels <- extract.chromosome.levels(seg.CNA)
  seg.CNA[,Chromosome:=ordered(toupper(Chromosome), chromosome.levels)]
  seg.CNA
}

#' Read the markers file in TSV format.
#' 
#' @param markers.file The name of the file containig the markers.
#' @param header A boolean value that specifies whether the file contain a header line.
#' 
#' @return A \code{data.table} containing the markers.
#' @export
read.markers <- function(markers.file, header=T) {
  if (file.access(markers.file, mode=4) == -1)
    stop(paste('The markers file', markers.file, 'cannot be found or read'))
  markers <- setNames(fread(markers.file, sep='\t', select=c(1,2,3), header=header), c('Name', 'Chromosome', 'Position'))
  # Make sure chromosomes stay ordered independently of their alphabetical order
  chromosome.levels <- extract.chromosome.levels(markers)
  unique(markers[,Chromosome:=ordered(toupper(Chromosome), chromosome.levels)], by=c('Chromosome', 'Position'))
}

#' Preprocess input file to crate a matrix containing the log ratio of each segment.
#' 
#' @param seg.CNA The \code{data.table} containing a segmented CNA as produced by \code{read.seg.file}.
#' @param markers The \code{data.table} containing the markers as produced by \code{read.markers}.
#' @param samples A vector containig the list of samples that awshould be analysed.
#' @param min.seg.markers The minimum number of markers per segment.
#' @param min.mean Log ratio means will be clipped at this mimimum value.
#' @param max.mean Log ratio means will be clipped at this maximum value.
#' @param ignore.chr A vector containing the list of chromosomes to ignore in the analysis.
#' @param spacing The spacing value to be added between chromosome in the absolute positioning.
#' 
#' @return A \code{data.table} containing all the samples preprocessed
#' @noRd
preprocess.map.loc <- function(seg.CNA, markers, samples=NULL,
                               min.seg.markers=1, min.mean=NA, max.mean=NA,
                               spacing=1e9) {

  setkey(seg.CNA, Sample)
  setkey(markers, Chromosome, Position)
  
  samples <- unique(samples)
  
  # Filter the matrix using the provided vector of samples
  if (length(samples) > 0) {
    seg.CNA <- seg.CNA[samples]
  }
  if (NROW(seg.CNA) == 0)
    stop('There are no samples left to analyse in seg.cna')
  
  # Create the location mapping
  # Counter necessary to number the probes
  p.counter <- 0
  # Join data.tables for all chromosomes adding an incremental probe name
  map.loc <- rbindlist(Map(function(x) {
    # Set all unique position pulling together all starts and ends in seg.CNA for all chromosomes
    positions <- unique(sort(c(seg.CNA[Chromosome==x]$Start, seg.CNA[Chromosome==x]$End)))
    # Generate the names of the probes as factors to keep an order which is not alphabetical
    p.num <- length(positions)
    probe.names <- paste0('P', p.counter+1:p.num)
    p.counter <<- p.counter + p.num
    probes <- ordered(probe.names, levels=probe.names)
    data.table(Probe=probes, Chromosome=x, Position=positions)
  }, unique(seg.CNA$Chromosome)))
  
  # NOTE: Specs show allow chromosome filtering, but chromosome filtering
  #       for CNA data type appears disabled in original sources
  
  # Generates a new column 'AbsPosition' containing an incremental
  # position over different chromosomes. Between each chromosome
  # an extra 'spacing' value is added  
  map.loc[,AbsPosition:=-spacing]
  # FIXME: Use data.table update rather than sapply
  sapply(unique(map.loc$Chromosome), function(x) {
    map.loc.start <- max(map.loc[,AbsPosition]) + spacing
    map.loc[Chromosome == x, AbsPosition:=Position + map.loc.start]
  })
  
  # Generate all segments assigning a unique segment id derived
  # by hashing Sample, Chromosome, Start, and End
  seg.CNA[, Segment:=unlist(Map(digest, Map(c, Sample, Chromosome, Start, End)))]
  
  # Get the number of markers for each segment
  setkey(seg.CNA, Chromosome, Start, End)
  # In order to use less memory the segment data set is grouped by samples first.
  # Then, each marker is joined with the segment where it belongs usig foverlaps.
  # Finally all the the markers beloning to the segments of the samples
  # are grouped together to count the number of markers per segment.
  seg.markers <- seg.CNA[, foverlaps(markers[, .(Name, Chromosome, Position, Position2=Position)],
                                     copy(.SD),
                                     by.x=c('Chromosome', 'Position','Position2'),
                                     type='within', nomatch=0L)[, .(ProbesNo=.N), by=Segment],
                         .SDcols=c('Chromosome', 'Start', 'End', 'Segment'),
                         by=Sample][, Sample:=NULL]

  # Join the sample segments with location map
  map.loc <- seg.CNA[,unique(foverlaps(map.loc[,.(Probe, Chromosome, AbsPosition, Position, Position2=Position)],
                                       copy(.SD),
                                       by.x=c('Chromosome','Position','Position2'),
                                       type='within'),
                             by='Probe', fromLast=T), by=Sample][, Position2:=NULL]
  
  # When samples have different segments replace the segments that are present
  # in some samples but missing in others with default values
  map.loc[is.na(LogRatio), c('Start','End','Segment'):=list(Position, Position, unlist(Map(digest, Map(c, Sample, Chromosome, Position, Position))))]
  # Update the segments group with the ones missing in some samples and
  # assign them 0 probes count to make sure they will always be marked
  # as CNV = 1 and therefore joined to other segments.
  seg.markers <- rbind(seg.markers, map.loc[is.na(LogRatio), .(Segment, ProbesNo=rep(0,.N))])
  map.loc[is.na(LogRatio), LogRatio:=rep(0,.N)]
  
  # Getting ready to join using the segment id
  setkey(map.loc, Segment)
  # Set CNV = 1 to each segment in which the number of markers is < min.seg.markers
  map.loc <- map.loc[seg.markers][, CNV:=unlist(Map(function(x) if(x < min.seg.markers) 1 else 0, ProbesNo))]
  setcolorder(map.loc, c('Sample', 'Segment', 'ProbesNo', 'Chromosome', 'Start', 'End', 'Probe', 'Position', 'AbsPosition', 'LogRatio', 'CNV'))
  
  # Join together segments which have less than min.seg.markers markers
  #####################################################################
  # NOTE: there is a bug in the original code when the last segment   #
  #       in a cromosome and the first segment in the next chromosome #
  #       having CNV == 1 they get joined together.                   #
  #       Furthermore, in the MATLAB code is only checked whether the #
  #       boundaries are on the chromosome boundaries, but because of #
  #       this bug, the chromosome boundaries can be inside the       #
  #       segment itself, and the check will fail to detect the jump  #
  #       to another choromosome.                                     #
  #####################################################################
  setkey(map.loc, Sample, AbsPosition)
  segment.grp <- map.loc[,T,by=.(Segment, Chromosome, CNV)][,.(Segment, SegmentGrp=cumsum(as.numeric(head(c(FALSE, c(tail(paste0(CNV, Chromosome),-1), NA) != paste0(CNV, Chromosome)), -1)) | !CNV) + CNV[1])]
  setkey(map.loc, Segment)
  map.loc <- map.loc[segment.grp]
  setkey(map.loc, SegmentGrp)
  # Update the LogRatio with the mean of the SegmentGrp
  map.loc[map.loc[,.(LogRatio=mean(LogRatio)), by=SegmentGrp], LogRatio:=i.LogRatio]
  setkey(map.loc, Sample, AbsPosition)
  
  plr.pointer <- NA_real_
  
  segment.lr <- map.loc[,.(CNV=CNV[1], Chromosome=Chromosome[1], LogRatio=LogRatio[1]), by=SegmentGrp][,.(SegmentGrp, LogRatio=unlist(Map(function(grp, cnv, lr, chr){
    # If the segment has to be joined
    if (cnv == 1) {
      pchr <- Chromosome[match(grp - 1, SegmentGrp)]
      nchr <- Chromosome[match(grp + 1, SegmentGrp)]
      # Use the previous pointer if the last log ratio belonged to a CNV==1 segment
      # otherwise use the log ratio of the previous segment
      plr <- if (is.na(plr.pointer)) LogRatio[match(grp - 1, SegmentGrp)] else plr.pointer
      # Get the log ratio of the next segment
      nlr <- LogRatio[match(grp + 1, SegmentGrp)]
      # Check if the segment is at the genome or chromosome boundaries
      if (!is.na(plr) & pchr != chr) plr <- NA
      if (!is.na(nlr) & nchr != chr) nlr <- NA
      if (is.na(plr) & is.na(nlr)) {
        lr <- 0
      } else {
        # Set the current log ratio to the closest value between the neighbors
        lrs <- c(nlr, plr)
        lrs <- lrs[!is.na(lrs)]
        lr <- lrs[which.min(abs(lrs - lr))]
      }
      plr.pointer <<- lr
    } else {
      plr.pointer <<- NA_real_
    }
    lr
  }, SegmentGrp, CNV, LogRatio, Chromosome)))]
  setkey(map.loc, SegmentGrp)
  map.loc[segment.lr, LogRatio:=i.LogRatio]
  map.loc[,c('Segment', 'ProbesNo', 'Start', 'End', 'CNV', 'SegmentGrp'):=NULL]
  setkey(map.loc, Sample, AbsPosition)
  
  # If min.mean and max.mean are set clip log ratio to these values
  if (!is.na(min.mean)) map.loc[,LogRatio:=sapply(LogRatio, max, min.mean)]
  if (!is.na(max.mean)) map.loc[,LogRatio:=sapply(LogRatio, min, max.mean)]
  
  map.loc
}

#' Returns the sum of the LogRatio of each probe across all samples.
#' 
#' @param map.loc The \code{data.table} containg segmented location mapping as returned by \code{preprocess.map.loc}.
#' 
#' @return The \code{data.table} containing the sum of the LogRatio of each probe across all samples.
#' @noRd
###
# MATLAB MAPPING
# sum.map.loc(map.loc)[,LogRatio] data.cnaAgr
sum.map.loc <- function(map.loc) {
  map.loc[,.(LogRatio=sum(LogRatio)), by=.(Probe, Chromosome, Position, AbsPosition)]
}


#' Extract a data matrix from the \code{map.loc data.table}.
#' 
#' @param map.loc The \code{data.table} containg segmented location mapping as returned by \code{preprocess.map.loc}.
#' @param column The \code{map.loc} column to be used as matrix values.
#' @noRd
extract.matrix <- function(map.loc, column='LogRatio') {
  data.matrix(dcast.data.table(map.loc, Sample ~ Probe, value.var=column)[,-1,with=F])
}


#' Returns the number of probes per sample.
#' 
#' @param map.loc The \code{data.table} containg segmented location mapping as returned by \code{preprocess.map.loc}.
#' 
#' @return The number of probes per sample.
#' @noRd
probes.per.sample <- function(map.loc) {
  map.loc[,.(.N),by=Sample][1]$N
}

#' Returns the number of samples.
#' 
#' @param map.loc The \code{data.table} containg segmented location mapping as returned by \code{preprocess.map.loc}.
#' 
#' @return The number of samples.
#' @noRd
samples.num <- function(map.loc) {
  length(data[,.GRP,by=Sample][,Sample])
}

#' Return the maximum number of probes per chromosome.
#' 
#' @param map.loc The \code{data.table} containg segmented location mapping as returned by \code{preprocess.map.loc}.
#' 
#' @return The maximum number of probes per chromosome.
#' @noRd
max.chrom.length <- function(map.loc) {
  map.loc[,.N,by=.(Sample,Chromosome)][,max(N)]
}

#' Ensures that there is a minimum number of probes before the analysis.
#' 
#' @param map.loc The \code{data.table} containg segmented location mapping as returned by \code{preprocess.map.loc}.
#' @param markers The \code{data.table} containing the markers as produced by \code{read.markers}.
#' @param min.probes The wante minimum number of probes. However, the minimum
#'                   number of probes that may result is \code{min(P, min.probes)}
#'                   (where P is the total number of probes in the dataset).
#'                   Setting \code{min.probes} lower than total number of probes
#'                   will result in improved computation time at the risk of
#'                   sacrificing statistical power.
#' @noRd
ensure.min.probes <- function(map.loc, markers, min.probes=2.6e5) {
  probes.num <- probes.per.sample(map.loc)
  # The data cannot have less than probes.num, therefore if min.probes has less than probes.num
  # return without changing anything
  if (min.probes < probes.num) {
    return(map.loc)
  }
  # NOTE: Original code here is bugged since it counts also the header line
  markers.num <- NROW(markers)
  # Reduce resolution removing some of the markers to fit the minimum of min.probes
  # NOTE: The original code do not order by chromosome position
  setkey(markers, Chromosome, Position)
  if (markers.num > min.probes) {
    probes.del <- markers.num / (min.probes - probes.num)
    probe.idx <- round(seq(max(probes.del, 1), markers.num, probes.del))
    probe.idx <- unique(probe.idx)
    markers <- markers[probe.idx]
  } else {
    markers <- copy(markers)
  }
  markers[,Name:=NULL]
  # Check that markers are within chromosome segments boundaries
  markers <- markers[setkey(map.loc[, .(lower=min(Position), upper=max(Position)), by=Chromosome], Chromosome)][between(Position, lower, upper)]
  # Create a new set of markers for each sample
  new.markers <- rbindlist(lapply(unique(map.loc[,Sample]), function(x) {
    data.table(Sample=x, markers)
  }))
  setkey(new.markers, Sample, Chromosome, Position)
  setkey(map.loc, Sample, Chromosome, Position)
  # Merge the new locations with the old ones
  map.loc <- merge(map.loc, new.markers, all=T)
  setkey(map.loc, Sample, Chromosome, Position)
  # Recompute probe names as factors to keep an order which is not alphabetical
  probe.names <- paste0('P', 1:map.loc[,.N,by=Sample][1]$N)
  map.loc[,Probe:=ordered(probe.names, levels=probe.names)]
  
  # Interpolate values based on the available positions
  interpolate.column <- function(positions, columnn.values) {
    # If there are no NAs there is no need for interpolation
    na.idx <- which(is.na(columnn.values))
    if (length(na.idx) > 0) {
      columnn.values[na.idx] <- interp1(positions[-na.idx], columnn.values[-na.idx], positions[na.idx], method='linear')
    }
    columnn.values
  }
  map.loc[,.(Probe,
             Position,
             AbsPosition=interpolate.column(Position, AbsPosition),
             LogRatio=interpolate.column(Position, LogRatio)), by=.(Sample, Chromosome)]
  }
