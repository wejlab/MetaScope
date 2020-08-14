#' A universal parameter settings object for Rsubread alignment
#'
#' This object is a named vector of multiple options that can be chosen for
#' functions that involve alignment with Rsubread, namely \code{align_target()}
#' and \code{filter_host()}. Both functions take an object for the parameter
#' \code{settings}, which are provided by \code{align_details} by default, or may
#' be given by a user-created object containing the same information.
#'
#' The default options included in \code{align_details} are \code{type = "dna"},
#' \code{nthreads = 8}, \code{maxMismatches = 5}, \code{nsubreads = 10},
#' \code{phredOffset = 33}, \code{unique = FALSE}, and
#' \code{nBestLocations = 16}. Full descriptions of these parameters can be
#' read by acessing \code{?Rsubread::align}.
#'
#' @keywords datasets
#'
#' @examples
#' data("align_details")
#'
"align_details"


