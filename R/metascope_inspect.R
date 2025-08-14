#' Inspect a Bowtie index
#'
#' This function can be use to call the \code{bowtie2-inspect} wrapper which
#' wraps the \code{bowtie2-inspect-s} and \code{bowtie2-inspect-l} binaries.
#'
#'
#'
#'
#'
#'
#'

metascope_inspect <- function(bt2_base, ...) {

    # Handle additional arguments
    arg_options <- c(...)
    if (is.null(arg_options)) return(NULL)
    arg_options <- trimws(arg_options)
    arg_options <- arg_options[nzchar(arg_options)]
    if (length(arg_options) == 0) return(NULL)

    Rbowtie2::checkPathExist(bt2_base, "bt2_base")


    # Combine arguments together
    bt2_base <- shQuote(bt2_base)
    args <- c(arg_options, bt2_base)

    # Call bowtie2-inspect wrapper using Rbowtie .callbinary
    tryCatch(
    {
        invisible(Rbowtie2::.callbinary(lang = "python3",
                                        bin1 = "bowtie2-inspect",
                                        args1 = paste(args, collapse = " ")))
    },
    error = function(e) {
        invisible(Rbowtie2::.callbinary(lang = "python",
                                        bin1 = "bowtie2-inspect",
                                        args1 = paste(args, collapse = " ")))
    },
    warning = function(w) {
    },
    finally = {
    }
    )
}