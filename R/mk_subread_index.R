#' Make a Subread index
#'
#' This function is a wrapper for the \code{Rsubread::buildindex} function. It
#' will generate one or more Subread indexes from a .fasta file. If the library
#' is too large (default >4GB) it will automatically be split into multiple
#' indexes, with _1, _2, etc at the end of the ref_lib basename.
#'
#' @param ref_lib The name/location of the reference library file, in
#'   (uncompressed) .fasta format.
#' @param split The maximum allowed size of the genome file (in GB). If the
#'   \code{ref_lib} file is larger than this, the function will split the
#'   library into multiple parts.
#' @param mem The maximum amount of memory (in MB) that can be used by the index
#'   generation process (used by the Rsubread::buildindex function).
#'
#' @return Creates one or more Subread indexes for the supplied reference .fasta
#'   file. If multiple indexes are created, the libraries will be named the
#'   \code{ref_lib} basename + "_1", "_2", etc. The function returns the names of the
#'   folders holding these files.
#'
#' @export
#'
#' @examples
#' #### Create a subread index from the example reference library
#' now <- Sys.time()
#' ## Create a temporary directory to store the reference library
#' ref_temp <- tempfile()
#' dir.create(ref_temp)
#'
#' ## Download reference genome
#' out_fasta <- download_refseq('Duck circovirus', reference = FALSE,
#'                              representative = FALSE, out_dir = ref_temp,
#'                              compress = TRUE, patho_out = FALSE,
#'                              caching = TRUE)
#'
#' ## Make subread index of reference library
#' mk_subread_index(out_fasta)
#' unlink(ref_temp)
#' now - Sys.time()
#'

mk_subread_index <- function(ref_lib, split = 4, mem = 8000) {
  gb <- 1073741824
  ref_size <- file.info(ref_lib)$size
  split_libs <- ceiling(ref_size / (split * gb))
  if (ref_size > (split * gb)) {
    message("Library size is ", round(ref_size / gb, 1),
    " GBs. Splitting the library into ", split_libs,
      " separate components.")
    # Making connections to the .fasta library
    # generating the new split .fasta files
    con <- file(ref_lib, open = "r")
    out_cons <- paste("outcon_", seq_len(split_libs), sep = "")
    out_files <- paste(tools::file_path_sans_ext(ref_lib),
      "_", seq_len(split_libs), ".", tools::file_ext(ref_lib),
      sep = "")
    for (i in seq_len(split_libs)) {
      assign(out_cons[i], file(out_files[i], open = "w"))
    }
    ## Reading ref library and splitting to split libraries
    ngenomes <- -1
    while ((length(one_line <-
                   readLines(con, n = 1, warn = FALSE)) > 0)) {
      if (substr(one_line, 1, 1) == ">") {
        ngenomes <- ngenomes + 1
      }
      writeLines(one_line, get(out_cons[(ngenomes %% split_libs) + 1]))
    }
    message("Printed ", ngenomes + 1, " sequences to ", split_libs,
            " genome files")
    close(con)
    for (ocon in out_cons) close(get(ocon))
    for (lib in out_files) Rsubread::buildindex(
      basename = tools::file_path_sans_ext(lib),
      reference = lib, memory = mem)
  } else {
    message("Library size is ", round(ref_size / gb, 1),
            " GBs. Building the Rsubread index")
    Rsubread::buildindex(
      basename = tools::file_path_sans_ext(ref_lib, compression = TRUE),
      reference = ref_lib, memory = mem)
  }
  return(
    paste(tools::file_path_sans_ext(ref_lib, compression = TRUE),
    "_", seq_len(split_libs), sep = ""))
}
