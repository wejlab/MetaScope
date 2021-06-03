#' Make a Bowtie2 index
#' 
#' This function is a wrapper for the \code{Rbowtie2::bowtie2_build} function.
#' It will create either small (.bt2) or large Bowtie2 indexes (.bt2l) of the 
#' given reference fasta files depending on their combined file size.
#' 
#' @param ref_dir \code{Character} scalar. The name/location of the directory that contains the reference
#' libraries in uncompressed fasta format
#' @param index_dir \code{Character} scalar. The name/location of the directory where the index files 
#' should be created
#' @param index_name \code{Character} scalar. The basename of the index files to be created. Default is the
#' basename of the index_dir
#' @param threads \code{Integer}. The number of threads to be used. Default is 4.
#' @param overwrite \code{Logical}. Whether existing files should be overwritten. Default is FALSE.
#' 
#' @return Creates multiple Bowtie2 indexes (.bt2 or .bt2l) of the supplied 
#' reference .fasta files. Returns the file path of the folder containing these
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