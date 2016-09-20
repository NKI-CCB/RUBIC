#' Segmented CNA
#'
#' Segmented CNA simulated data used by RUBIC example.
#'
#' @format \code{\link[data.table]{data.table}} with 5 columns:
#' \describe{
#' \item{Sample}{Sample anme}
#' \item{Chromosome}{Chromosome name}
#' \item{Start}{Segment start}
#' \item{End}{Segment end}
#' \item{LogRatio}{Log ratios for this segment}
#' }
#' @examples
#'   seg.cna
"seg.cna"

#' Markers locations
#'
#' The markers locations used by RUBIC example that indicates the exact
#' locations of measurement probes (markers) for the given platform.
#' For sequencing data, copy number values are often estimated with fixed
#' bin sizes (prior to segmentations). In this case each marker should be
#' associated with a bin and the center genomic position of the bin.
#'
#' @format \code{\link[data.table]{data.table}} with 3 columns:
#' \describe{
#' \item{Name}{Probe name}
#' \item{Chromosome}{Chromosome name}
#' \item{Position}{Probe position on the chromosome}
#' }
#' @examples
#'   markers
"markers"

#' Samples names
#'
#' Samples names used to filter samples in the \code{\link{seg.cna}} data for
#' the RUBIC example.
#'
#' @format \code{\link[data.table]{data.table}} with 1 column:
#' \describe{
#' \item{V1}{Sample name}
#' }
#' @examples
#'   samples
"samples"

#' Gene locations
#'
#' Gene locations used by RUBIC example to call focal events.
#'
#' @format \code{\link[data.table]{data.table}} with 5 columns:
#' \describe{
#' \item{ID}{Gene ID}
#' \item{Name}{Gene symbol}
#' \item{Chromosome}{Chromosome name}
#' \item{Start}{Gene start}
#' \item{End}{Gene end}
#' }
#' @examples
#'   genes
"genes"