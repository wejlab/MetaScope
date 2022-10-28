#' @name check_samtools_exists
#' @title Check if samtools exists on the system
#' @description This is an internal function that is not meant to be used
#'   outside of the package. It checks whether samtools exists on the system.
#' @return Returns TRUE if samtools exists on the system, else FALSE.
#'

check_samtools_exists <- function() {
  if (file.exists(Sys.which("samtools"))) return(TRUE)
  return(FALSE)
}
