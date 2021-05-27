#' Make a Bowtie2 index
#' 
#' This function is a wrapper for the \code{Rbowtie2::bowtie2_build} function.
#' It will create large Bowtie2 indexes (.bt2l) from a reference fasta files.
#' 
#' @param ref_dir The name/location of the directory that contains the reference
#' libraries in uncompressed fasta format
#' @param index_dir The name/location of the directory where the index files 
#' should be created
#' @param index_name The basename of the index files. Default is the name of 
#' the last directory specified in index_dir 
#' @param threads The number of threads to be used. Default is 4.
#' @param overwrite Whether existing files should be overwritten. Default is false.
#' 
#' @return Creates multiple large Bowtie2 indexes (.bt2l) for the supplied 
#' reference .fasta files. Returns the file path of the folder holding these
#' files.
#' 
#' @export


mk_bowtie_index <- function(ref_dir, index_dir, index_name = basename(index_dir), threads = 4, overwrite=FALSE){
  
  ref_dir <- dir(ref_dir, full.names = TRUE)
  index_dir <- tools::file_path_as_absolute(index_dir)
  
  Rbowtie2::bowtie2_build(references = ref_dir, 
                          bt2Index=file.path(index_dir, index_name),
                          paste("--threads ",threads),
                          overwrite = overwrite
                          )
  
  return(paste("Bowtie2 index files located at ", tools::file_path_as_absolute(index_dir)))
}