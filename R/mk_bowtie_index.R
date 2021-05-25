#' Make a Bowtie2 index
#' 
#' This function is a wrapper for the \code{Rbowtie2::bowtie2_build} function.
#' It will create large Bowtie2 indexes (.bt2l) from a fasta file.
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
#' @export


mk_bowtie_index <- function(ref_dir, index_dir, index_name = basename(index_dir), threads = 4, overwrite=FALSE){
  Rbowtie2::bowtie2_build(references = ref_dir, 
                          bt2lIndex=file.path(index_dir, index_name),
                          paste("--threads ",threads),
                          overwrite = overwrite
                          )
}