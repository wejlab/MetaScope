#' Align microbiome reads to set of indexed Bowtie2 libraries
#' 
#' @param read1 Path to the .fastq file to align. 
#' @param read2 Optional: Location of the mate pair .fastq file to align.
#' @param lib_dir Path to the directory that contains the Bowtie2 indexes.
#' @param libs The basename of the Bowtie2 indexes to align against 
#' (without trailing .bt2 or .bt2l extensions).
#' @param align_dir Path to the directory where the output alignment file 
#' should be created.
#' @param align_file The basename of the output alignment file file 
#' (without trailing .sam or .bam extensions). 
#' @param align_format The format of the alignment file. Default is "bam" 
#' but can also pass "sam" to the function. NOTE: To create the bam file, first
#' the sam file must be created but afterwards the sam file will be deleted. 
#' @param  bowtie2_options Optional: Additional parameters that can be passed to
#' the align_target_bowtie() function. To see all the available parameters
#' use Rbowtie2::bowtie2_usage(). Default parameters are the parameters are the 
#' default parameters that PathoScope 2.0 uses. NOTE: Users should pass all their
#' parameters as one string and if optional parameters are given then the user 
#' is responsible for entering all the parameters to be used by Bowtie2. NOTE:
#' The only parameters that should NOT be specified here is the threads.
#' @param threads The number of threads that can be utilized by the function.
#' Default is 8 threads.
#' @param overwrite Whether existing files should be overwritten. 
#' Default is FALSE.
#' 
#' @return Returns the path to the directory where the output alignment 
#' file is stored. 
#' 
#' @export
#' 
#' @examples
#' # Code not run
#' \dontrun{
#' 
#' ## Create alignment bam file produced by Bowtie2
#' 
#' # Create a temporary directory to store reference fasta files
#' ref_temp <- tempfile()
#' dir.create(ref_temp)
#' 
#' # Create a temporary directory to store the index files 
#' lib_temp <- tempfile()
#' dir.create(lib_temp)
#' 
#' # Create a temporary directory to store the index files 
#' align_temp <- tempfile()
#' dir.create(align_temp)
#' 
#' # Download all the RefSeq reference viral genomes to current directory
#' download_refseq('viruses', compress = TRUE)
#' 
#' # Move downloaded fasta file to temporary reference directory 
#' file.rename(from = file.path(".", "Viruses.fasta.gz"), 
#' to = file.path(ref_temp, "Viruses.fasta.gz"))
#' 
#' # Create bowtie index files in temporary index directory
#' mk_bowtie_index(ref_dir = ref_temp, lib_dir = lib_temp, lib_name = "virus", 
#' overwrite=FALSE)
#' 
#' # Get path to example reads
#' readPath <- system.file("extdata", "virus_example.fastq", package = "MetaScope")
#' 
#' # Create alignment file 
#' align_target_bowtie(read1 = readPath, lib_dir = lib_temp, 
#' libs = "virus", align_dir = align_temp, align_file = "example", 
#' align_format = "bam", overwrite = TRUE)
#' }

align_target_bowtie <- function(read1, read2 = NULL, 
                                lib_dir, libs,
                                align_dir, align_file, 
                                align_format = "bam",
                                bowtie2_options = NULL, threads = 8,
                                overwrite = FALSE){
  
  
  # Convert user specified paths to absolute paths for debugging purposes
  lib_dir <- tools::file_path_as_absolute(lib_dir)
  align_dir <- tools::file_path_as_absolute(align_dir)
  
  # If user does not specify their own parameters then use default 
  # PathoScope parameters
  if (missing(bowtie2_options)){
    bowtie2_options <- paste("--very-sensitive-local -k 100 --score-min L,20,1.0",
                             "--threads",threads)
  }
  else
    bowtie2_options <- paste(bowtie2_options,"--threads",threads)
  
  
  # If user does not specify mate pair reads then align using Bowtie2 unpaired 
  # alignment  
  if (missing(read2)){
    
    message("Attempting to perform Bowtie2 unpaired alignment")
    
    Rbowtie2::bowtie2(bt2Index = file.path(lib_dir, libs), 
                      output = file.path(align_dir, align_file),
                      outputType = align_format,
                      seq1 = read1,
                      overwrite = overwrite, 
                      bowtie2_options)
    
    message(paste("Successfully created",paste0(align_file,".",align_format)))
    
    # Sort .bam file and remove unmapped reads using filter_unmapped_reads()
    # if format is set to "bam"
    if (align_format == "bam"){
      
      message(paste("Attempting to filter unmapped reads from",
                    paste0(align_file,".",align_format)))
      
      bam_location <- file.path(align_dir,paste0(align_file,".bam"))
      filter_unmapped_reads(bam_location)
      
      message(paste("Successfully filtered unmapped reads from",
                    paste0(align_file,".",align_format)))
    }
    
    
    return(align_dir)
    
  }
  
  # If user specifies mate pair reads then use Bowtie2 pair-end alignment 
  else{
    
    message("Attempting to perform Bowtie2 paired-end alignment")
    
    Rbowtie2::bowtie2(bt2Index = file.path(lib_dir, libs), 
                      output = file.path(align_dir, align_file),
                      outputType = align_format,
                      seq1 = read1,
                      seq2 = read2,
                      overwrite = overwrite, 
                      bowtie2_options)
    
    message(paste("Successfully created",paste0(align_file,".",align_format)))
    
    # Sort .bam file and remove unmapped reads using filter_unmapped_reads()
    # if format is set to "bam"
    if (align_format == "bam"){
      
      message(paste("Attempting to filter unmapped reads from",
                    paste0(align_file,".",align_format)))
      
      bam_location <- file.path(align_dir,paste0(align_file,".bam"))
      filter_unmapped_reads(bam_location)
      
      message(paste("Successfully filtered unmapped reads from",
                    paste0(align_file,".",align_format)))
    }
    
    
    return(align_dir)
    
  }
}