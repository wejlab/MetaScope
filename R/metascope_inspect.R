#' Inspect a Bowtie index
#'
#' This function can be use to call the \code{bowtie2-inspect} wrapper which
#' wraps the \code{bowtie2-inspect-s} and \code{bowtie2-inspect-l} binaries.
#'
#' @param bt2_base The path to the Bowtie2 index files minus minus trailing
#' .1.bt2/.2.bt2 extension names.
#' @param ... Required additional arguments to be passed on to the \code{bowtie2-inspect}
#' binaries. See below for details.
#'
#' @details All additional arguments in ... are intepreted as
#' parameters to be passed on to \code{bowtie2-inspect} wrapper. All of them
#' \code{Character} or \code{Numeric} scalar. You can put all additional
#' arguments in one \code{Character} (e.g. "--across 60 --names") with white
#' space separation, or put them in different \code{Character}
#' (e.g. "--across","60","--names"). See the output of
#' \code{metascope_inspect_usage()} for details about available parameters.
#'
#' @references
#' {
#' Langmead B, Wilks C, Antonescu V, Charles R. Scaling
#' read aligners to hundreds of threads on general-purpose processors.
#' Bioinformatics. 2018 Jul 18. doi: 10.1093/bioinformatics/bty648.
#'
#' Langmead B, Salzberg SL. Fast gapped-read alignment with
#' Bowtie 2. Nature Methods. 2012 Mar 4;9(4):357-9. doi: 10.1038/nmeth.1923.
#' }
#' @export metascope_inspect
#' @examples
#' \dontrun{
#' metascope_inspect(bt2_base)
#' }



metascope_inspect <- function(bt2_base, ...) {

    # Handle additional arguments
    arg_options <- c(...)
    if (is.null(arg_options)) return(NULL)
    arg_options <- trimws(arg_options)
    arg_options <- arg_options[nzchar(arg_options)]
    if (length(arg_options) == 0) return(NULL)

    # Combine arguments together
    bt2_base <- shQuote(bt2_base)
    args <- c(arg_options, bt2_base)

    # Call bowtie2-inspect wrapper
    script_path <- system.file("bowtie2-inspect", package = "MetaScope")
    if (script_path == "" || !file.exists(script_path)) {
        stop("bowtie2-inspect script not found in package")
    }

    tryCatch({
        cmd <- paste("python3", shQuote(script_path), paste(args, collapse = " "))
        system(cmd)
    },
    error = function(e) {
        cmd <- paste("python", shQuote(script_path), paste(args, collapse = " "))
        system(cmd)
    })


}

#' @name metascope_inspect_usage
#' @title Print available arguments that can be passed to metascope_inspect()
#' @description Calling metascope_inspect_usage() prints the available arguments
#' that can be passed to the ... argument of the metascope_inspect() function of
#' the package. Note that some arguments are invalid if they are already
#' handled as explicit function arguments.
#' @return Information about available arguments that can be passed to
#' metascope_inspect()
#' @references
#' {
#' Langmead B, Wilks C, Antonescu V, Charles R. Scaling
#' read aligners to hundreds of threads on general-purpose processors.
#' Bioinformatics. 2018 Jul 18. doi: 10.1093/bioinformatics/bty648.
#'
#' Langmead B, Salzberg SL. Fast gapped-read alignment with
#' Bowtie 2. Nature Methods. 2012 Mar 4;9(4):357-9. doi: 10.1038/nmeth.1923.
#' }
#' @export metascope_inspect_usage
#' @examples
#' metascope_inspect_usage()


metascope_inspect_usage <- function() {
    if(R.Version()$arch=="i386"){
        return("bowtie2 is not available for 32bit, please use 64bit R instead")
    }

    script_path <- system.file("bowtie2-inspect", package = "MetaScope")
    if (script_path == "" || !file.exists(script_path)) {
        stop("bowtie2-inspect script not found in package")
    }

    tryCatch({
        cmd <- paste("python3", shQuote(script_path), "-h")
        system(cmd)
    },
    error = function(e) {
        cmd <- paste("python", shQuote(script_path), "-h")
        system(cmd)
    })
}