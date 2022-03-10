globalVariables(c("align_details"))

#' Helper function for filter_host_bowtie to write interim fastq file
mk_interim_fastq <- function(bf, read_loc, maxMemory) {
    message("Creating Intermediate Fastq file")
    # Obtain read names
    open(bf)
    innames <- Rsamtools::scanBam(bf, param = Rsamtools::ScanBamParam(
        what = c("qname")), maxMemory = maxMemory) %>% .[[1]] %>%
        unlist(.$qnames) %>% unname()
    # Remove repeated reads from subsequent pulled info
    ind <- !duplicated(innames)
    # Read in quality scores
    open(bf)
    inqual <- Rsamtools::scanBam(bf, param = Rsamtools::ScanBamParam(
        what = c("qual")), maxMemory = maxMemory) %>% .[[1]] %>%
        .$qual %>% as.character() %>% .[ind]
    # Appending sequences without saving (to conserve memory)
    open(bf)
    Rsamtools::scanBam(bf, param = Rsamtools::ScanBamParam(
        what = c("seq")), maxMemory = maxMemory) %>% .[[1]] %>% .$seq %>%
        as.character() %>% .[ind] %>% dplyr::tibble(
            Header = innames[ind], Sequence = ., Quality = inqual) %>%
        # Keep first occurrences of reads
        dplyr::mutate(Plus = "+",
                      Header = stringr::str_c("@", .data$Header)) %>%
        dplyr::select(.data$Header, .data$Sequence, .data$Plus,
                      .data$Quality) %>%
        as.matrix() %>% t() %>% as.character() %>% dplyr::as_tibble() %>%
        data.table::fwrite(., file = read_loc, compress = "gzip",
                           col.names = F, quote = F)
    message("Finished Creating Intermediate Fastq File")
}

#' Helper function to remove reads matched to filter libraries
#'
#' Within the \code{filter_host()} function, we align our sequencing sample to all
#' filter libraries of interest. The \code{remove_matches()} function allows
#' for removal of any target reads that are also aligned to filter libraries.
#'
#' This function is not intended for use by users.
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
#' @examples
#' 
#' #readPath <- system.file("extdata", "subread_target.bam",
#' #                        package = "MetaScope")
#'
#' ## Assume that the first 10 query names aligned to first filter library
#' ## And another 10 aligned to second filter library
#' # qnames <- Rsamtools::scanBam(readPath)[[1]]$qname
#' # read_names <- list(qnames[1:10], qnames[30:40])
#' # out <- "subread_target.filtered.bam"
#' 
#' # remove_matches(readPath, read_names, out)
#'

remove_matches <- function(reads_bam, read_names, name_out) {
    # Note: reads_BAM and filter-aligned files are already sorted by chromosome
    # index bam file
    bam_index <- Rsamtools::indexBam(reads_bam)
    # obtain vector of target query names from .bam file
    target_reads <- Rsamtools::scanBam(reads_bam)[[1]]$qname
    # some aligned reads may be duplicated; remove these, and unlist
    filter_reads <- unique(unlist(read_names))
    # define logical vector of which reads to keep, based on query names
    filter_which <- !(target_reads %in% filter_reads)
    # create BamFile instance to set yieldSize
    bf <- Rsamtools::BamFile(reads_bam, yieldSize = length(target_reads))
    filtered_bam <- Rsamtools::filterBam(
        bf, destination = name_out, index = bam_index,
        indexDestination = FALSE, filter = filter_which,
        param = Rsamtools::ScanBamParam(what = "qname"))
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
#' #### Filter reads from bam file that align to any of the filter libraries
#' 
#' ## Assuming a bam file has been created previously with align_target()
#'
#' ## Create object with path to the example filter library
#' refPath <- system.file("extdata","filter.fasta", package = "MetaScope")
#' 
#' ## Copy the example filter library to the current directory
#' file.copy(from = refPath, to = file.path(".", "filter.fasta"))
#' 
#' ## Make subread index of filter library
#' mk_subread_index('filter.fasta')
#'
#' ## Create object with path to the previously aligned bam file
#' readPath <- system.file("extdata", "subread_target.bam", package = "MetaScope")
#'
#' ## Filter bam file 
#' filter_host(readPath, libs = "filter")
#'

filter_host <- function(reads_bam, libs, lib_dir = NULL,
                        output = paste(tools::file_path_sans_ext(reads_bam),
                                       "filtered", "bam", sep = "."),
                        settings = align_details) {
    # Initialize list of names
    read_names <- vector(mode = "list", length(libs))

    for (i in seq_along(libs)) {
        # Create output file name for BAM
        lib_file <- paste(tools::file_path_sans_ext(reads_bam), ".",
                          libs[i], ".bam", sep = "")
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
        # Throw away BAM, vcf file
        file.remove(lib_file)
        file.remove(paste(lib_file, ".indel.vcf", sep = ""))
        file.remove(paste(lib_file, ".summary", sep = ""))
    }

    # helper function to sort headers and filter BAM file
    remove_matches(reads_bam, read_names, output)
    # output final filtered BAM file
    message("DONE! Alignments written to ", output)
    return(output)
}

#' Align reads against one or more filter libraries and subsequently
#' remove mapped reads
#' 
#' After a sample is aligned to a target library with \code{align_target_bowtie()},
#' we may use \code{filter_host_bowtie()} to remove unwelcome host contamination using
#' filter reference libraries. This function takes as input the name of the .bam 
#' file produced via \code{align_target_bowtie()}, and produces a
#' sorted .bam file with any reads that match the filter libraries removed.
#' This resulting .bam file may be used downstream for further analysis.
#' 
#' @param reads_bam The name of a merged, sorted .bam file that has previously
#' been aligned to a reference library. Likely, the output from running an
#' instance of \code{align_target_bowtie()}.
#' @param lib_dir Path to the directory that contains the filter Bowtie2 index
#' files.
#' @param libs The basename of the filter libraries 
#' (without .bt2 or .bt2l extension)
#' @param output The desired name of the output .bam file. Default is
#' the basename of \code{unfiltered_bam} + \code{.filtered.bam}.
#' @param bowtie2_options Optional: Additional parameters that can be passed to
#' the filter_host_bowtie() function. To see all the available parameters
#' use Rbowtie2::bowtie2_usage(). Default parameters are the parameters are the 
#' default parameters that PathoScope 2.0 uses. NOTE: Users should pass all their
#' parameters as one string and if optional parameters are given then the user 
#' is responsible for entering all the parameters to be used by Bowtie2. NOTE:
#' The only parameters that should NOT be specified here is the threads.
#' @param threads The amount of threads available for the function.
#' Default is 8 threads.
#' @param overwrite Whether existing files should be overwritten. 
#' Default is FALSE.
#' @param maxMemory The number of MB of RAM that the sortBAM function can use.
#' The smaller the number, the more temporary files will be produced, which will
#' hopefully save memory usage. The default is 512 MB.
#' 
#' @return The name of a filtered, sorted .bam file written to the user's
#' current working directory.
#' 
#' @export
#' 
#' @examples
#' #### Filter reads from bam file that align to any of the filter libraries
#' 
#' ## Assuming a bam file has already been created with align_target_bowtie()
#' 
#' ## Create a temporary directory to store the filter library
#' ref_temp <- tempfile()
#' dir.create(ref_temp)
#' 
#' ## Create a temporary directory to store the filter library index files
#' lib_temp <- tempfile()
#' dir.create(lib_temp)
#' 
#' ## Create a temporary directory to store the filtered bam file
#' align_temp <- tempfile()
#' dir.create(align_temp)
#' 
#' ## Create object with path to previously created bam file 
#' bamPath <- system.file("extdata", "bowtie_target.bam", package = "MetaScope")
#' 
#' ## Create object with path to the filter library
#' refPath <- system.file("extdata","filter.fasta", package = "MetaScope")
#' 
#' ## Move the filter library to the temporary reference directory 
#' file.copy(from = refPath, to = file.path(ref_temp, "filter.fasta"))
#' 
#' ## Create the bowtie index files in the temporary index directory
#' mk_bowtie_index(ref_dir = ref_temp, lib_dir = lib_temp, lib_name = "filter",
#' overwrite=FALSE)
#' 
#' ## Filter reads from the bam file that align to the filter library
#' filter_host_bowtie(reads_bam = bamPath, lib_dir = lib_temp, libs = "filter")
#'

filter_host_bowtie <- function(reads_bam, lib_dir, libs,
                               output = paste(
                                   tools::file_path_sans_ext(reads_bam),
                                   "filtered", "bam", sep = "."),
                               bowtie2_options = NULL, threads = 8,
                               overwrite = FALSE, maxMemory = 512) {
    # If no optional parameters are passed then use default parameters else use user parameters 
    if (missing(bowtie2_options)) {
        bowtie2_options <- paste("--very-sensitive-local -k 100",
                                 "--score-min L,20,1.0", "--threads", threads)
    } else bowtie2_options <- paste(bowtie2_options,"--threads", threads)
    # Convert reads_bam into fastq (Default 1,000,000,000 read chunks) 
    bf <- Rsamtools::BamFile(reads_bam, yieldSize = 1000000000)
    read_loc <- file.path(dirname(reads_bam), "intermediate.fastq")
    # Make Intermediate Fastq file
    mk_interim_fastq(bf, read_loc, maxMemory)
    read_names <- vector(mode = "list", length(libs)) # Init list of names
    for (i in seq_along(libs)) {
        # Create output file name for BAM
        lib_file <- paste(tools::file_path_sans_ext(reads_bam),".",
                          libs[i], ".bam", sep = "")
        # Align reads to lib and generate new filter BAM file
        Rbowtie2::bowtie2_samtools(
            bt2Index = file.path(lib_dir,libs[i]),
            output = tools::file_path_sans_ext(lib_file),
            outputType = "bam", seq1 = read_loc, ... = bowtie2_options,
            overwrite = overwrite)
        # sort BAM file and remove umapped reads (package helper function)
        filter_unmapped_reads(lib_file)
        # Extract target query names from mapped BAM file
        read_names[[i]] <- Rsamtools::scanBam(lib_file)[[1]]$qname
        # Throw away BAM file
        file.remove(lib_file)
    }
    # remove intermediate fastq file
    file.remove(read_loc)
    # helper function to sort headers and filter BAM file
    remove_matches(reads_bam, read_names, output)
    # output final filtered BAM file
    message("DONE! Alignments written to ", output)
    return(output)
}
