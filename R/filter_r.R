#' Helper function to remove reads matched to filter libraries
#' 
#' Within the `filter_r()` function, we align our sequencing sample to all
#' filter libraries of interest. The `remove_matches()` function allows
#' use the header 
#' 
#' @param reads_bam
#' @param header_names
#' 
#' @return
#' 
#' @examples
#' 
remove_matches <- function(reads_bam, header_names){
  
  # Note: reads_BAM and filter-aligned files are already sorted by chromosome
  this_header <- Rsamtools::scanBamHeader(reads_bam,
                                            what = c("text"))[[1]]$text
  # Extract all header info for reference sequence dictionary
  SQ_ind <- names(this_header) %in% "@SQ"
  target_header <- sapply(this_header[SQ_ind], function(x) x[1])
  
  # header_names is a list potentially, need to do sapply or loop
  for (k in seq_along(header_names)) {
    header_names[[k]] %in% target_header
  }

  # Go into original BAM file and remove anything with same read name
  # To do this - check helper function in align_target
  
  # rename bam file?
  
  # To do:
  # See if there are "SQ" that are present in both, then filter?
  
}


#' Align reads against one or more filter libraries and subsequently
#' remove mapped reads
#'
#' Although we have already filtered out any unmapped reads, which may belong
#' to one or more host species or otherwise, there may still remain some sort
#' of unwelcome contamination in the data from the filter species which we
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
#' @return This function overwrites the inputted BAM file and removes any
#' contamination from species in the filter genome library.
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
    Rsubread::align(index = libs[i], readfile1 = reads_bam,
                    input_format = "bam",
                    output_file = lib_file,
                    type = "dna", nthreads = threads,
                    unique = FALSE, nBestLocations = 16,
                    maxMismatches = mismatch)
    # sort BAM file and remove umapped reads (package helper function)
    filter_unmapped_reads(lib_file)
    
    # Function - filter out read names
    
   
    # throw away BAM file
    file.remove(lib_file)
    # throw away vcf file
    file.remove(paste(lib_file, ".indel.vcf", sep = ""))
    # keep summary file for now
  }
  
  # Helper function to sort headers and filter BAM file
  remove_matches(header_names)
  
  # output that filtered BAM file
  message(paste("DONE! Alignments written to ", project_name, ".bam",
                sep = ""))
  return(paste(project_name, ".bam", sep = ""))
}


filtered_bam <- Rsamtools::filterBam(sorted_bamfile, destination = bamfile,
                                     index = bam_index, 
                                     indexDestination = FALSE,
                                     param = Rsamtools::ScanBamParam(
                                       flag = Rsamtools::scanBamFlag(
                                         isUnmappedQuery = F)))
# flag: is query in list?
# remember to index file
# cross-check lists, insert T/F flag?
# How to delete a set of reads from a BAM file