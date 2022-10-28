#' A universal parameter settings object for Rsubread alignment
#'
#' This object is a named vector of multiple options that can be chosen for
#' functions that involve alignment with Rsubread, namely \code{align_target()}
#' and \code{filter_host()}. Both functions take an object for the parameter
#' \code{settings}, which are provided by \code{align_details} by default, or
#' may be given by a user-created object containing the same information.
#'
#' The default options included in \code{align_details} are \code{type = "dna"},
#' \code{maxMismatches = 3}, \code{nsubreads = 10}, \code{phredOffset = 33},
#' \code{unique = FALSE}, and \code{nBestLocations = 16}. Full descriptions of
#' these parameters can be read by accessing \code{?Rsubread::align}.
#' @name align_details
#' @docType data
#' @format list
#' @usage data(align_details)
#' @keywords datasets
#' @examples
#' data("align_details")
#'
"align_details"

#' A universal parameter object for Bowtie 2 16S alignment
#'
#' This character string provides several Bowtie 2 options to provide an
#' optimized alignment specifically optimized for 16S amplicon sequencing data.
#' This object can be used with functions that use the Bowtie 2 aligner through
#' the \code{Rbowtie2} package, namely \code{align_target_bowtie()} and
#' \code{filter_host_bowtie}. These settings can be substituted for default
#' settings by passing \code{} to the \code{bowtie2_options} argument.
#'
#' The default parameters listed in this object are "-local -R 2 -N 0 -L 25 -i
#' S,1,0.75 -k 10 --score-min L,100,1.28".
#'
#' Further delineation of Bowtie 2 parameters is provided in the
#' \href{https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml}{Bowtie 2
#' manual}.
#'
#' These parameters were chosen from the article
#' \href{https://doi.org/10.1101/2022.07.27.501757}{Metagenomic profiling
#' pipelines improve taxonomic classification for 16S amplicon sequencing data},
#' Faits et al. [Preprint].
#'
#' @name bt2_16S_params
#' @docType data
#' @format list
#' @usage data(bt2_16S_params)
#' @keywords datasets
#' @examples
#' data("bt2_16S_params")
#' 
"bt2_16S_params"

#' A universal parameter object for Bowtie 2 loose alignment
#'
#' This character string provides several Bowtie 2 options to provide a loose
#' alignment useful for metagenomes. This object can be used with functions that
#' use the Bowtie 2 aligner through the \code{Rbowtie2} package, namely
#' \code{align_target_bowtie()} and \code{filter_host_bowtie}. These settings
#' can be substituted for default settings by passing \code{} to the
#' \code{bowtie2_options} argument.
#'
#' The default parameters listed in this object are "--local -k 100 -D 20 -R 3
#' -L 3 -N 1 -p 8 --gbar 1 --mp 3".
#'
#' Further delineation of Bowtie 2 parameters is provided in the
#' \href{https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml}{Bowtie 2
#' manual}.
#'
#' These parameters were chosen from the article
#' \href{https://doi.org/10.1101/2022.07.27.501757}{bowtie2: Relaxed Parameters
#' for Generous Alignments to Metagenomes}.
#'
#' @name bt2_loose_params
#' @docType data
#' @format list
#' @usage data(bt2_loose_params)
#' @keywords datasets
#' @examples
#' data("bt2_loose_params")
#' 
"bt2_loose_params"


#' Taxonomy table
#'
#' This is a taxonomy table to be used in \code{download_refseq()}. Last updated
#' 10/18/22.
#'
#' @keywords datasets
#' 
"taxonomy_table"