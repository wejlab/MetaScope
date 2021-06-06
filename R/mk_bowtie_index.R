#' Make a Bowtie2 index
#' 
#' This function is a wrapper for the \code{Rbowtie2::bowtie2_build} function.
#' It will create either small (.bt2) or large Bowtie2 indexes (.bt2l) depending
#' on the combined size of the reference fasta files.
#' 
#' @param ref_dir \code{Character} scalar. The location of the directory that contains the reference
#' files in uncompressed fasta format. This directory should contain only the reference fasta files to be indexed.
#' @param index_dir \code{Character} scalar. The location of the directory where the Bowtie2 index files 
#' should be created.
#' @param index_name \code{Character} scalar. The basename of the Bowtie2 index files to be created.
#' @param ... \code{Character} scalar. Optional arguments that can be passed to the mk_bowtie_index() function. 
#' All options should be passed as one string. To see all the available options that can be passed to the function 
#' use Rbowtie2::bowtie2_build_usage().
#' @param overwrite \code{Logical}. Whether existing files should be overwritten. Default is FALSE.
#' 
#' @return Creates multiple Bowtie2 indexes (.bt2 or .bt2l) of the supplied 
#' reference .fasta files. Returns the file path of the folder containing these
#' files.
#' 
#' @export
#' 
#' @examples  
#' # Code not run
#' \dontrun{
#' ## Download all RefSeq reference viral genomes and make an index
#' 
#' # Create a temporary directory to store reference fasta files
#' ref_temp <- tempdir()
#' 
#' # Download all the RefSeq reference viral genomes
#' download_refseq('viruses', compress = FALSE)
#' 
#' # Move downloaded fasta file to temporary directory 
#' file.rename(from = file.path(".", "Viruses.fasta"), to = file.path(ref_temp, "Viruses.fasta"))
#' 
#' # Create bowtie index files in current directory
#' mk_bowtie_index(ref_dir = ref_temp, index_dir = ".", index_name = "virus", "--threads 4", overwrite=FALSE)
#' }


mk_bowtie_index <- function(ref_dir, index_dir, index_name, ..., overwrite=FALSE){
  
  ref_dir <- dir(ref_dir, full.names = TRUE)
  index_dir <- tools::file_path_as_absolute(index_dir)
  
  
  Rbowtie2::bowtie2_build(references = ref_dir, 
                          bt2Index=file.path(index_dir, index_name),
                          ...,
                          overwrite = overwrite
                          )
  
  return(tools::file_path_as_absolute(index_dir))
}