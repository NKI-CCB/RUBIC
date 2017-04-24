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


#' Return the aggregate profile of amplifications.
#'
#' @param map.loc
#' @param amp.level
#'
#' @noRd
amp.profile <- function(map.loc, amp.level) {
    map.loc[,.(LogRatio=sum(LogRatio[LogRatio >= amp.level])), by=.(Probe, Chromosome, Position, AbsPosition)]
}


#' Return the aggregate profile of deletions.
#'
#' @param map.loc
#' @param del.level
#'
#' @noRd
del.profile <- function(map.loc, del.level) {
  map.loc[,.(LogRatio=sum(LogRatio[LogRatio <= del.level])), by=.(Probe, Chromosome, Position, AbsPosition)]
}


focal.events.as.data.table <- function(focal.events) {
  rbindlist(lapply(focal.events, function(event) {
    list(ID=event$I, Chromosome=event$chromosome,
         Start=event$loc.start, End=event$loc.end)
  }))
}


#' Transforms the list of segments in a data.table.
#'
#' @param map.loc.agr
#' @param segments
#' @param focal.events
#' @param markers
#'
#' @return A \code{data.table} containing segment information and whether they are
#'         a focal event.
#'
#' @noRd
segments.as.data.table <- function(map.loc.agr, segments, markers) {
  dt <- rbindlist(lapply(segments, function(segment) {
    locs <- map.loc.agr[segment$I:(segment$I + segment$kw - 1),
                        c(.SD[c(1,NROW(.SD)),.(Chromosome, Position, AbsPosition)],
                          LogRatio=mean(LogRatio))]

    loc.chr <- as.character(locs[1, Chromosome])
    loc.start <- locs[1, Position]
    loc.end <- locs[2, Position]
    loc.abs.start <- locs[1, AbsPosition]
    loc.abs.end <- locs[2, AbsPosition]

    if (!is.null(markers)) {
      # If there are not marker found by the query a warning is displayed and
      # -/+ Inf is returned
      loc.adj.start <- suppressWarnings(markers[Chromosome==loc.chr & Position < loc.start,
                                                max(Position)])
      loc.adj.end <- suppressWarnings(markers[Chromosome==loc.chr & Position > loc.end,
                                              min(Position)])
      if (is.infinite(loc.adj.start)) {
        loc.adj.start <- loc.start
        loc.abs.adj.start <- loc.abs.start
      } else {
        # NOTE: For plotting purposes we are using the mean, while in the focal.events
        #       gene.locs we use the maximum extension
        loc.adj.start <- mean(c(loc.start, loc.adj.start))
        loc.abs.adj.start <- mean(c(loc.abs.start, map.loc.agr[Position==loc.adj.start, AbsPosition]))
      }
      if (is.infinite(loc.adj.end)) {
        loc.adj.end <- loc.end
        loc.abs.adj.end <- loc.abs.end
      } else {
        # NOTE: For plotting purposes we are using the mean, while in the focal.events
        #       gene.locs we use the maximum extension
        loc.adj.end <- mean(c(loc.end, loc.adj.end))
        loc.abs.adj.end <- mean(c(loc.abs.end, map.loc.agr[Position==loc.adj.end, AbsPosition]))
      }
    } else {
      loc.adj.start <- loc.start
      loc.adj.end <- loc.end
      loc.abs.adj.start <- loc.abs.start
      loc.abs.adj.end <- loc.abs.end
    }

    list(ID=segment$I, Chromosome=loc.chr,
         Start=loc.adj.start, End=loc.adj.end,
         StartAbs=loc.abs.adj.start, EndAbs=loc.abs.adj.end,
         Mu=locs[1, LogRatio])
  }))
  dt[,Chromosome:=ordered(toupper(Chromosome), levels(map.loc.agr$Chromosome))]
}


#' Format number without using scientific notation.
#'
#' @param x A vector of numbers.
#' @return A vector of strings of formatted numbers.
#' 
#' @noRd
location.format <- function(x) {
  format(x, scientific=F, trim=T)
}

#' Draws a striped bar of gene locations.
#'
#' @param genes A \code{data.table} contaning at least Start and End columns
#'              which identify gene locations.
#' @param min.x The lower limit of the x axis.
#' @param max.x The upper limit of the x axis.
#' @param focal.p An optional \code{data.table} as returned by
#'                \code{focal.events.as.data.table} for the positive focal events.
#' @param focal.n An optional \code{data.table} as returned by
#'                \code{focal.events.as.data.table} for the negative focal events.
#'
#' @noRd
plot.genes.bar <- function(genes, min.x, max.x, focal.p=NULL, focal.n=NULL) {
  plot <- ggplot(genes)
  if (NROW(focal.p) > 0)
    plot <- plot + geom_rect(data=focal.p,
                             aes(xmin=Start, ymin=0, xmax=End, ymax=1),
                             alpha=.7, size=NA, fill='red')
  if (NROW(focal.n) > 0)
    plot <- plot + geom_rect(data=focal.n,
                             aes(xmin=Start, ymin=0, xmax=End, ymax=1),
                             alpha=.7, size=NA, fill='blue')
  plot <- plot + geom_rect(aes(xmin=Start, xmax=End, ymin=0, ymax=1), fill='grey45') +
    labs(x = NULL, y = NULL) +
    scale_x_continuous(expand=c(0, 0), limits=c(min.x, max.x)) +
    scale_y_continuous(expand=c(0, 0)) +
    theme(panel.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.line=element_blank(),
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          legend.text.align=0,
          plot.margin=unit(c(0,0,0,0), "lines"))
  plot
}


#' Plots aggreate log ratios, segmented profile, focal abberrations, and optionally
#' locations of selected genes.
#'
#' @param amp.profile A \code{data.table} as returned by \code{amp.profile}.
#' @param del.profile A \code{data.table} as returned by \code{del.profile}.
#' @param segments.p An optional \code{data.table} as returned by
#'                   \code{segments.as.data.table} for the positive segments.
#' @param segments.n An optional \code{data.table} as returned by
#'                   \code{segments.as.data.table} for the negative segments.
#' @param focal.p An optional \code{data.table} as returned by
#'                   \code{focal.events.as.data.table} for the positive focal events.
#' @param focal.n An optional \code{data.table} as returned by
#'                   \code{focal.events.as.data.table} for the negative focal events.
#' @param genes An optional \code{data.table} contaning at least two colums named
#'              'Start' and 'End' corresponding to the gene positions.
#' @param steps A boolean value indicating whether to plot the segments as steps
#'              or horizontal segments.
#' @param chromosome The chromosome to be plotted.
#' 
#' @noRd
plot.cna <- function(amp.profile=NULL, del.profile=NULL,
                     segments.p=NULL, segments.n=NULL,
                     focal.p=NULL, focal.n=NULL,
                     genes=NULL, steps=T, chromosome=NULL, offset.table=NULL) {

  min.x <- +Inf
  max.x <- -Inf

  #TODO: if (is.null(chromosome))
  #      transform local coordinates in absolute
  #      and plot the whole genome


  if (NROW(amp.profile) > 0) {
    if (!is.null(chromosome)) {
      amp.profile <- amp.profile[Chromosome==chromosome]
    }
    min.x <- suppressWarnings(min(min.x, amp.profile[,min(Position)]))
    max.x <- suppressWarnings(max(max.x, amp.profile[,max(Position)]))
  }

  if (NROW(del.profile) > 0) {
    if (!is.null(chromosome)) {
      del.profile <- del.profile[Chromosome==chromosome]
    }
    min.x <- suppressWarnings(min(min.x, del.profile[,min(Position)]))
    max.x <- suppressWarnings(max(max.x, del.profile[,max(Position)]))
  }

  if (NROW(segments.p) > 0) {
    if (!is.null(chromosome)) {
      segments.p <- segments.p[Chromosome==chromosome]
    }
    min.x <- suppressWarnings(min(min.x, segments.p[,min(Start)]))
    max.x <- suppressWarnings(max(max.x, segments.p[,max(End)]))
  }

  if (NROW(segments.n) > 0) {
    if (!is.null(chromosome)) {
      segments.n <- segments.n[Chromosome==chromosome]
    }
    min.x <- suppressWarnings(min(min.x, segments.n[,min(Start)]))
    max.x <- suppressWarnings(max(max.x, segments.n[,max(End)]))
  }

  if (NROW(focal.p) > 0) {
    if (!is.null(chromosome)) {
      focal.p <- focal.p[Chromosome==chromosome]
    }
    min.x <- suppressWarnings(min(min.x, focal.p[,min(Start)]))
    max.x <- suppressWarnings(max(max.x, focal.p[,max(End)]))
  }

  if (NROW(focal.n) > 0) {
    if (!is.null(chromosome)) {
      focal.n <- focal.n[Chromosome==chromosome]
    }
    min.x <- suppressWarnings(min(min.x, focal.n[,min(Start)]))
    max.x <- suppressWarnings(max(max.x, focal.n[,max(End)]))
  }

  if (is.infinite(min.x) || is.infinite(max.x))
    stop('There is no data to plot.')

  if (NROW(genes) > 0) {
    if (!is.null(chromosome)) {
      genes <- genes[Chromosome==chromosome]
    }
    min.x <- suppressWarnings(min(min.x, genes[,min(Start)]))
    max.x <- suppressWarnings(max(max.x, genes[,max(End)]))
  }

  if (NROW(genes) > 0) {
    top.plot <- plot.genes.bar(genes, min.x, max.x, focal.p, focal.n)
  }

  if (is.null(chromosome)){ 
    if(!is.null(amp.profile) && !is.null(del.profile)) {
      if (NROW(amp.profile) > 0) {
        amp.profile[,Position := as.numeric(Position)]
        amp.profile[,Position := .SD[,Position]+offset.table[.(.SD[,'Chromosome']),'max',mult='first',with=FALSE][,max], by=Chromosome, .SDcols=c('Chromosome', 'Position')]
        min.x <- suppressWarnings(min(min.x, amp.profile[,min(Position)]))
        max.x <- suppressWarnings(max(max.x, amp.profile[,max(Position)]))
      }
      if (NROW(segments.p) > 0) {
        segments.p[,Start := .SD[,Start]+offset.table[.(.SD[,'Chromosome']),'max',mult='first',with=FALSE][,max], by=Chromosome, .SDcols=c('Chromosome', 'Start')]
        segments.p[,End := .SD[,End]+offset.table[.(.SD[,'Chromosome']),'max',mult='first',with=FALSE][,max], by=Chromosome, .SDcols=c('Chromosome', 'End')]
        min.x <- suppressWarnings(min(min.x, segments.p[,min(Start)]))
        max.x <- suppressWarnings(max(max.x, segments.p[,max(End)]))
      }
      if (NROW(focal.p) > 0) {
        focal.p[,Start := as.numeric(Start)]
        focal.p[,Start := .SD[,Start]+offset.table[.(.SD[,'Chromosome']),'max',mult='first',with=FALSE][,max], by=Chromosome, .SDcols=c('Chromosome', 'Start')]
        focal.p[,End := as.numeric(End)]
        focal.p[,End := .SD[,End]+offset.table[.(.SD[,'Chromosome']),'max',mult='first',with=FALSE][,max], by=Chromosome, .SDcols=c('Chromosome', 'End')]
        min.x <- suppressWarnings(min(min.x, focal.p[,min(Start)]))
        max.x <- suppressWarnings(max(max.x, focal.p[,max(End)]))
      }
      
      if (NROW(del.profile) > 0) {
        del.profile[,Position := as.numeric(Position)]
        del.profile[,Position := .SD[,Position]+offset.table[.(.SD[,'Chromosome']),'max',mult='first',with=FALSE][,max], by=Chromosome, .SDcols=c('Chromosome', 'Position')]
        min.x <- suppressWarnings(min(min.x, del.profile[,min(Position)]))
        max.x <- suppressWarnings(max(max.x, del.profile[,max(Position)]))
        
      }
      if (NROW(segments.n) > 0) {
        segments.n[,Start := .SD[,Start]+offset.table[.(.SD[,'Chromosome']),'max',mult='first',with=FALSE][,max], by=Chromosome, .SDcols=c('Chromosome', 'Start')]
        segments.n[,End := .SD[,End]+offset.table[.(.SD[,'Chromosome']),'max',mult='first',with=FALSE][,max], by=Chromosome, .SDcols=c('Chromosome', 'End')]
        min.x <- suppressWarnings(min(min.x, segments.n[,min(Start)]))
        max.x <- suppressWarnings(max(max.x, segments.n[,max(End)]))
        
      }
      if (NROW(focal.n) > 0) {
        focal.n[,Start := as.numeric(Start)]
        focal.n[,Start := .SD[,Start]+offset.table[.(.SD[,'Chromosome']),'max',mult='first',with=FALSE][,max], by=Chromosome, .SDcols=c('Chromosome', 'Start')]
        focal.n[,End := as.numeric(End)]
        focal.n[,End := .SD[,End]+offset.table[.(.SD[,'Chromosome']),'max',mult='first',with=FALSE][,max], by=Chromosome, .SDcols=c('Chromosome', 'End')]
        min.x <- suppressWarnings(min(min.x, focal.n[,min(Start)]))
        max.x <- suppressWarnings(max(max.x, focal.n[,max(End)]))
        
      }
    }  
 
  }
  
  main.plot <- ggplot(NULL, aes(x=Position, y=LogRatio))

  if(NROW(amp.profile) > 0)
    main.plot <- main.plot + geom_point(data=amp.profile, color='#d6a9a9')
  if(NROW(del.profile))
    main.plot <- main.plot + geom_point(data=del.profile, color='#a9a9d6')

  if ((NROW(del.profile) == 0 && NROW(segments.n) == 0 && NROW(focal.n) == 0) ||
       (NROW(amp.profile) == 0 && NROW(segments.p) == 0 && NROW(focal.p) == 0)) {
    rect.min.y <- 0
    rect.max.y <- +Inf
    if (NROW(amp.profile) == 0) {
      main.plot <- main.plot + scale_y_reverse(name='Aggregate log ratios')
    }

  } else {
    main.plot <- main.plot + geom_hline(aes(yintercept=0), size=.2, color='grey', linetype="dashed") +
      scale_y_continuous(name='Aggregate log ratios')
    rect.min.y <- -Inf
    rect.max.y <- +Inf
  }

  if(NROW(focal.p) > 0)
    main.plot <- main.plot + geom_rect(data=focal.p, aes(xmin=Start, xmax=End), ymin=rect.min.y, ymax=rect.max.y, alpha=.2, size=NA, fill='red', inherit.aes=F)
  if(NROW(focal.n) > 0)
    main.plot <- main.plot + geom_rect(data=focal.n, aes(xmin=Start, xmax=End), ymin=rect.min.y, ymax=rect.max.y, alpha=.2, size=NA, fill='blue', inherit.aes=F)
  if(NROW(segments.p) > 0) {
    if (steps) {
      setkey(segments.p, End)
      main.plot <- main.plot +
        # Add one extra row to the data.table that contain the end of the segment
        # as its start in order to complete the last step of the ladder.
        # It assumes that the last row is the the last segment
        geom_step(data=rbind(segments.p,
                             segments.p[NROW(segments.p),
                                        modifyList(.SD, list(Start=End))]),
                  aes(x=Start, y=Mu), size=.6, color='red')
    }else{
      main.plot <- main.plot + geom_segment(data=segments.p, aes(x=Start, y=Mu, xend=End, yend=Mu), size=.6, color='red')
    }
  }
  if(NROW(segments.n) > 0) {
    if (steps) {
      setkey(segments.n, End)
      main.plot <- main.plot +
        # Add one extra row to the data.table that contain the end of the segment
        # as its start in order to complete the last step of the ladder.
        # It assumes that the last row is the the last segment
        geom_step(data=rbind(segments.n,
                             segments.n[NROW(segments.n),
                                        modifyList(.SD, list(Start=End))]),
                  aes(x=Start, y=Mu), size=.6, color='blue')
    } else {
      main.plot <- main.plot + geom_segment(data=segments.n, aes(x=Start, y=Mu, xend=End, yend=Mu), size=.6, color='blue')
    }
  }

  if (is.null(chromosome)) {
    main.plot <- main.plot + geom_vline(xintercept=offset.table[,max][-c(1)])
  }  
  
  main.plot <- main.plot + scale_x_continuous(labels=location.format, limits=c(min.x, max.x), expand=c(0, 0)) +
    theme(panel.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text=element_text(color='black'),
          axis.ticks=element_line(color='black'),
          axis.title.x=element_blank()
          )
  
  main.grob <- ggplotGrob(main.plot)

  if (NROW(genes) > 0) {

    top.grob <- ggplotGrob(top.plot)
    main.grob <- gtable_add_rows(main.grob, unit(.025, "npc"), 1)
    main.grob <- gtable_add_grob(main.grob, top.grob$grobs[[which(top.grob$layout$name == "panel")]], 2, 4, 2, 4)
  }
  # FIXME: The current version of ggplot on CRAN (1.0.1) don't support
  #        yet the possibility to use ggsave to save grobs,
  #        therefore we need to display it.
  grid.newpage()
  grid.draw(main.grob)
}


plot.cna.chromosome <- function(chromosome, amp.profile=NULL, del.profile=NULL,
                                segments.p=NULL, segments.n=NULL,
                                focal.p=NULL, focal.n=NULL,
                                genes=NULL, steps=T) {
  plot.cna.call <- as.list(match.call(expand.dots=FALSE))
  #steps <- plot.cna.call$steps
  #plot.cna.call$steps <- NULL
  # Removes function name and chromosome, and process the other arguments
  new.arguments <- lapply(plot.cna.call[c(-1,-2)], function(a) {
    paste0(a, '[Chromosome=="', chromosome ,'"]')
  })
  #plot.cna.call[[1]] <- 'plot.cna'
  plot.cna.call <- modifyList(plot.cna.call, new.arguments)
  plot.cna.call$steps <- steps
  return(eval(plot.cna, plot.cna.call, parent.frame()))
}


generate.all.plots.eps <- function(dir,
                                   map.loc, amp.level, del.level,
                                   segments.p, segments.n,
                                   focal.p.events, focal.n.events,
                                   markers, steps=T, genes=NULL,
                                   width=11, height=5) {

  dir.create(dir)
  
  a.prof <- amp.profile(map.loc, amp.level)
  d.prof <- del.profile(map.loc, del.level)
  focal.p <- focal.events.as.data.table(focal.p.events)
  focal.n <- focal.events.as.data.table(focal.n.events)
  segs.p <- segments.as.data.table(a.prof, segments.p, markers)
  segs.n <- segments.as.data.table(d.prof, segments.n, markers)

  extension <- 'eps'

  setEPS()
  for (chromosome in levels(map.loc[,Chromosome])) {
    a.file.name <- file.path(dir, paste0('chromosome_', chromosome, '_all.', extension))
    p.file.name <- file.path(dir, paste0('chromosome_', chromosome, '_gains.', extension))
    n.file.name <- file.path(dir, paste0('chromosome_', chromosome, '_losses.', extension))

    cairo_ps(a.file.name, width=width, height=height)
    plot.cna(amp.profile=a.prof, del.profile=d.prof,
             segments.p=segs.p, segments.n=segs.n,
             focal.p=focal.p, focal.n=focal.n,
             genes=genes, steps=steps, chromosome=chromosome)
    dev.off()

    cairo_ps(p.file.name, width=width, height=height)
    plot.cna(amp.profile=a.prof,
             segments.p=segs.p,
             focal.p=focal.p,
             genes=genes, steps=steps, chromosome=chromosome)
    dev.off()

    cairo_ps(n.file.name, width=width, height=height)
    plot.cna(del.profile=d.prof,
             segments.n=segs.n,
             focal.n=focal.n,
             genes=genes, steps=steps, chromosome=chromosome)
    dev.off()
  }
  
  offset.table <- markers[, .(max=max(Position)), by=Chromosome]
  offset.table[, max:=c(0,max[0:(length(max)-1)])]
  offset.table[,max := cumsum(max)]
  
  a.file.name <- file.path(dir, paste0('chromsome_all_all.', extension))  
  p.file.name <- file.path(dir, paste0('chromsome_all_gains.', extension))
  n.file.name <- file.path(dir, paste0('chromosome_all_losses.', extension))
  
  cairo_ps(a.file.name, width=width, height=height)
  plot.cna(amp.profile=a.prof, del.profile=d.prof,
           segments.p=segs.p, segments.n=segs.n,
           focal.p=focal.p, focal.n=focal.n,
           genes=NULL, steps=steps, offset.table=offset.table) # genes NULL so always whole genome
  dev.off()
  
  cairo_ps(p.file.name, width=width, height=height)
  plot.cna(amp.profile=a.prof,
           segments.p=segs.p,
           focal.p=focal.p,
           genes=NULL, steps=steps, offset.table=offset.table) # genes NULL so always whole genome
  dev.off()     
  
  cairo_ps(n.file.name, width=width, height=height)
  plot.cna(del.profile=d.prof,
           segments.n=segs.n,
           focal.n=focal.n,
           genes=NULL, steps=steps, offset.table=offset.table) # genes NULL so always whole genome
  dev.off()  
  
}
