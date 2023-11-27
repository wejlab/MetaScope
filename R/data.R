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
#' The default parameters listed in this object are
#' "--local -R 2 -N 0 -L 25 -i S,1,0.75 -k 5 --score-min L,0,1.88"
#' 
#' Note that k is actually 10 and is doubled internally from 5.
#' The \code{score-min} function was chosen such that the minimum alignment score
#' allowed requires 98% identity.
#'
#' Further delineation of Bowtie 2 parameters is provided in the
#' \href{https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml}{Bowtie 2
#' manual}.
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

#' A universal parameter object for Bowtie 2 metagenomic or non-16S alignment
#'
#' This character string provides several Bowtie 2 options to provide a 95% identity
#' alignment useful for metagenomes. This object can be used with functions that
#' use the Bowtie 2 aligner through the \code{Rbowtie2} package, namely
#' \code{align_target_bowtie()} and \code{filter_host_bowtie}. These settings
#' can be substituted for default settings by passing \code{} to the
#' \code{bowtie2_options} argument.
#'
#' The default parameters listed in this object are
#' "--local -R 2 -N 0 -L 25 -i S,1,0.75 -k 5 --score-min L,0,1.7".
#'
#' Further delineation of Bowtie 2 parameters is provided in the
#' \href{https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml}{Bowtie 2
#' manual}.
#'
#' @name bt2_regular_params
#' @docType data
#' @format list
#' @usage data(bt2_regular_params)
#' @keywords datasets
#' @examples
#' data("bt2_regular_params")
#'
"bt2_regular_params"

#' A universal parameter object for Bowtie 2 metagenomic alignment where the host genome is thought to be absent from the reference database
#'
#' This character string provides several Bowtie 2 options to conduct an 
#' alignment useful for metagenomes, especially in the case where a genome may
#' not be present in the reference database. This object can be used with functions that
#' use the Bowtie 2 aligner through the \code{Rbowtie2} package, namely
#' \code{align_target_bowtie()} and \code{filter_host_bowtie}. These settings
#' can be substituted for default settings by passing \code{} to the
#' \code{bowtie2_options} argument.
#'
#' The default parameters listed in this object are
#' "--local -R 2 -N 0 -L 25 -i S,1,0.75 -k 5 --score-min L,0,1.4".
#'
#' Further delineation of Bowtie 2 parameters is provided in the
#' \href{https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml}{Bowtie 2
#' manual}.
#'
#' @name bt2_missing_params
#' @docType data
#' @format list
#' @usage data(bt2_missing_params)
#' @keywords datasets
#' @examples
#' data("bt2_missing_params")
#'
"bt2_missing_params"
