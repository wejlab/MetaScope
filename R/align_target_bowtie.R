#' Align microbiome reads to set of indexed Bowtie2 libraries
#' 
#' This is the main MetaScope target library mapping function, using RBowtie2
#' and multiple libraries. Aligns to each library separately, filters
#' unmapped reads from each file, and then merges and sorts the .bam files
#' from each library into one output file. If desired, output can be
#' passed to `filter_host_bowtie()` to remove reads that also map to filter library
#' genomes.
#'
#' @param read1 Path to the .fastq file to align.
#' @param read2 Optional: Location of the mate pair .fastq file to align.
#' @param lib_dir Path to the directory that contains the Bowtie2 indexes.
#' @param libs The basename of the Bowtie2 indexes to align against
#' (without trailing .bt2 or .bt2l extensions).
#' @param align_dir Path to the directory where the output alignment file
#' should be created.
#' @param align_file The basename of the output alignment file
#' (without trailing .bam extension).
#' @param  bowtie2_options Optional: Additional parameters that can be passed to
#' the align_target_bowtie() function. To see all the available parameters
#' use Rbowtie2::bowtie2_usage(). Default parameters are the parameters are the
#' default parameters that PathoScope 2.0 uses. NOTE: Users should pass all their
#' parameters as one string and if optional parameters are given then the user
#' is responsible for entering all the parameters to be used by Bowtie2. NOTE:
#' The only parameter that should NOT be specified here is the threads.
#' @param threads The number of threads that can be utilized by the function.
#' Default is 8 threads.
#' @param overwrite Whether existing files should be overwritten.
#' Default is FALSE.
#'
#' @return Returns the path to where the output alignment file is stored.
#'
#' @export
#'
#' @examples
#' #### Align example reads to an example reference library using Rbowtie2
#'
#' ## Create a temporary directory to store reference library 
#' ref_temp <- tempfile()
#' dir.create(ref_temp)
#'
#' ## Create a temporary directory to store the reference library index files
#' lib_temp <- tempfile()
#' dir.create(lib_temp)
#'
#' ## Create a temporary directory to store the bam file
#' align_temp <- tempfile()
#' dir.create(align_temp)
#'
#' ## Create object with path to example reference library
#' refPath <- system.file("extdata","target.fasta", package = "MetaScope")
#'
#' ## Copy the reference library to the temporary directory 
#' file.copy(from = refPath, to = file.path(ref_temp, "target.fasta"))
#'
#' ## Create the bowtie index files in the temporary index library directory
#' mk_bowtie_index(ref_dir = ref_temp, lib_dir = lib_temp, lib_name = "target",
#' overwrite=FALSE)
#'
#' ## Create object with path to the example reads
#' readPath <- system.file("extdata", "reads.fastq", package = "MetaScope")
#'
#' ## Align reads to the reference library
#' align_target_bowtie(read1 = readPath, lib_dir = lib_temp,
#' libs = "target", align_dir = align_temp, align_file = "bowtie_target",
#' overwrite = TRUE)


align_target_bowtie <- function(read1, read2 = NULL, lib_dir, libs, align_dir,
                                align_file, bowtie2_options = NULL,
                                threads = 8, overwrite = FALSE) {

    # Convert user specified paths to absolute paths for debugging purposes
    lib_dir <- tools::file_path_as_absolute(lib_dir)
    align_dir <- tools::file_path_as_absolute(align_dir)
    # If user does not specify their own parameters then use default
    # PathoScope parameters
    if (missing(bowtie2_options)) {
        bowtie2_options <- paste(
            "--very-sensitive-local -k 100 --score-min L,20,1.0", "--threads",
            threads)
    } else bowtie2_options <- paste(bowtie2_options, "--threads", threads)

    bam_files <- numeric(length(libs))
    for (i in seq_along(libs)) {
        # Do not attach the .bam extension because Rbowtie2 does this already
        bam_files[i] <- 
            file.path(align_dir,
                      paste(basename(tools::file_path_sans_ext(read1)),
                            ".", libs[i], sep = ""))

        message("Attempting to perform Bowtie2 alignment on ",
                libs[i], " index")
        Rbowtie2::bowtie2_samtools(bt2Index = file.path(lib_dir, libs[i]),
                                   output = bam_files[i], outputType = "bam",
                                   seq1 = read1, seq2 = read2,
                                   overwrite = overwrite,
                                   ... = bowtie2_options)
        # Attach .bam extension to bam files in order to call this function
        filter_unmapped_reads(paste0(bam_files[i], ".bam"))
    }

    # Create variable names for files
    outputFile <- file.path(align_dir, paste0(align_file, ".bam"))
    bam_files <- paste0(bam_files, ".bam")
    # If more than one libraries were aligned to then combine the bam files
    if (length(bam_files) > 1) {
        message("Merging the bam files into ",align_file,".bam")
        merge_bam_files(bam_files, tools::file_path_sans_ext(outputFile))
    } else file.rename(bam_files, outputFile)
    
    message("DONE! Alignments written to ", outputFile)
    return(outputFile)
}
