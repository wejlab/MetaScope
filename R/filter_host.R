globalVariables(c("align_details"))
#' Helper function to remove reads matched to filter libraries
#' 
#' Within the \code{filter_host()} function, we align our sequencing sample to all
#' filter libraries of interest. The \code{remove_matches()} function allows
#' for removal of any target reads that are also aligned to filter libraries.
#' It is not intended for use by users.
#' 
#' @param reads_bam The name of a merged, sorted .bam file that has previously
#' been aligned to a reference library. Likely, the output from running an
#' instance of \code{align_target()}.
#' @param read_names A \code{list} of target query names from \code{reads_bam}
#' that have also aligned to a filter reference library. Each \code{list}
#' element should be a vector of read names.
#' @param name_out The name of the .bam file that to which the filtered alignments
#' will be written.
#' 
#' @return The name of a filtered, sorted .bam file written to the user's
#' current working directory.
#' 

remove_matches <- function(reads_bam, read_names, name_out) {
  # Note: reads_BAM and filter-aligned files are already sorted by chromosome
  # index bam file
  bam_index <- Rsamtools::indexBam(reads_bam)

  # obtain vector of target query names from .bam file
  target_reads <- Rsamtools::scanBam(reads_bam)[[1]]$qname

  # some aligned reads may be duplicated; remove these, and unlist filter names
  filter_reads <- unique(unlist(read_names))

  # define logical vector of which reads to keep, based on query names (qnames)
  filter_which <- !(target_reads %in% filter_reads)
  
  # create BamFile instance to set yieldSize
  bf <- Rsamtools::BamFile(reads_bam, yieldSize = length(filter_which))
  
  filtered_bam <- Rsamtools::filterBam(bf, destination = name_out,
                                       index = bam_index,
                                       indexDestination = FALSE,
                                       filter = filter_which,
                                       param = Rsamtools::ScanBamParam(
                                         what = "qname"))

  # clean up
  file.remove(bam_index)

  return(reads_bam)
}


#' Align reads against one or more filter libraries and subsequently
#' remove mapped reads
#'
#' After a sample is aligned to a target library with \code{align_target()},
#' we may use \code{filter_host()} to remove unwelcome host contamination using
#' filter reference libraries. This function takes as input the name
#' of the .bam file produced via \code{align_target()}, and produces a
#' sorted .bam file with any reads that match the filter libraries removed.
#' This resulting .bam file may be used upstream for further analysis.
#' It is not intended for use by users.
#'
#' @param reads_bam The name of a merged, sorted .bam file that has previously
#' been aligned to a reference library. Likely, the output from running an
#' instance of \code{align_target()}.
#' @param output The desired name of the output .bam file. Default is
#' the basename of \code{reads_bam} + \code{.filtered.bam}.
#' @inheritParams align_target
#'
#' @return The name of a filtered, sorted .bam file written to the user's
#' current working directory.
#' 
#' @export
#' 
#' @examples
#' # How to get previous BAM file??? Rerun code? Or include it?
#' # Create an index for filter libraries
#' # filter
#' 
#' #' ## Get reference and filter genome libraries
#' download_refseq('viral', compress = FALSE)
#' download_refseq('bacterial', compress = FALSE)
#'
#' ## Make and align to a single reference genome library
#' mk_subread_index('viral.fasta')
#' readPath <- system.file("extdata", "virus_example.fastq",
#'                          package = "MetaScope")
#' viral_map <- align_target(readPath, "viral", project_name = "virus_example")
#'
#' # Filter with separate reference genome library
#' mk_subread_index('bacterial.fasta')
#' filter_host(viral_map, libs = "bacterial")
#' 

filter_host <- function(reads_bam, libs, lib_dir=NULL,
                     output = paste(tools::file_path_sans_ext(reads_bam),
                                    "filtered", "bam", sep = "."),
                     settings = align_details) {
  # Initialize list of names
  read_names <- vector(mode = "list", length(libs))
  
  for (i in seq_along(libs)) {
    # Create output file name for BAM
    lib_file <- paste(tools::file_path_sans_ext(reads_bam),
                      ".", libs[i], ".bam", sep = "")
    # Align BAM to the lib & generate new file
    Rsubread::align(index = paste(lib_dir,libs[i],sep=""), 
                    readfile1 = reads_bam,
                    input_format = "bam",
                    output_file = lib_file,
                    type = settings[["type"]],
                    nthreads = settings[["nthreads"]],
                    maxMismatches = settings[["maxMismatches"]],
                    nsubreads = settings[["nsubreads"]],
                    phredOffset = settings[["phredOffset"]],
                    unique = settings[["unique"]],
                    nBestLocations = settings[["nBestLocations"]])
    # sort BAM file and remove umapped reads (package helper function)
    filter_unmapped_reads(lib_file)
    
    # Extract target query names from mapped BAM file
    read_names[[i]] <- Rsamtools::scanBam(lib_file)[[1]]$qname
    
    # throw away BAM, vcf file
    file.remove(lib_file)
    file.remove(paste(lib_file, ".indel.vcf", sep = ""))
    file.remove(paste(lib_file, ".summary", sep = ""))
  }
  
  # helper function to sort headers and filter BAM file
  remove_matches(reads_bam, read_names, output)
  
  # output final filtered BAM file
  message(paste("DONE! Alignments written to", output, sep = " "))
  
  return(output)
}
