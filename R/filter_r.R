#' Helper function to remove reads matched to filter libraries
remove_matches <- function(){
  
  # Note that output of align_R is already sorted by chromosome
  # Maybe then we sort the BAM files from filtering and do it that way?
  
  # Now use Rsamtools; steal code from align_target??
  # Go into original BAM file and remove anything with same read name
  # To do this - check helper function in align_target
  
  # rename bam file?
}


#' Align reads against one or more filter libraries and subsequently
#' remove mapped reads
#'
#' Although we have already filtered out any unmapped reads, which may belong
#' to one or more host species or otherwise, there may still remain some sort
#' of unwelcome contamination in the data from the filter species, which we
#' wish to remove. To do this, we employ `filter_r()`, which takes as input
#' the location of the BAM file created from `align_target()`, and produces a
#' sorted BAM file with any reads that match the filter libraries removed. We
#' will then use this final BAM file downstream for further analysis.
#'
#' @param reads_bam Location of the merged and sorted .bam file to align that
#' was previously output by the `align_target()` function
#' @inheritParams align_target
#' @param project_name
#' @param threads
#' @param mismatch
#'
#' @return
#' 
#' @examples
#' # How to get previous BAM file??? Rerun code? Or include it?
#' # Create an index for filter libraries
#' # filter
#' 

filter_r <- function(reads_bam, libs,
                     project_name = tools::file_path_sans_ext(reads_bam),
                     threads = 8, mismatch = 5) {
  # Initialize list of names
  header_names <- vector(mode = "list", length(libs))
  
  for (i in seq_along(libs)) {
    # Create output file name for BAM
    lib_file <- paste(project_name, ".", libs[i], ".bam", sep = "")
    # Align BAM to the lib & generate new file
    Rsubread::align(index = libs[i], readfile1 = reads,
                    input_format = "bam",
                    output_file = lib_file,
                    type = "dna", nthreads = threads,
                    unique = FALSE, nBestLocations = 16,
                    maxMismatches = mismatch)
    # Remove umapped reads
    filter_unmapped_reads(lib_file)
    # Save into a running list the BAM headers specifying the read names that
    # positively hit the current filter library
    header_names[[i]] <- scanBamHeader(lib_file)

    # throw away BAM file
    file.remove(lib_file)
    # throw away other files...
  }
  
  # Helper function to sort headers and filter BAM file
  remove_matches(header_names)
  
  # output that filtered BAM file
  message(paste("DONE! Alignments written to ", project_name, ".bam",
                sep = ""))
  return(paste(project_name, ".bam", sep = ""))
}
