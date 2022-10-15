globalVariables(c("align_details"))

# Helper function for mk_interim_fastq()
reduceByYield_iterate <- function(X, YIELD, MAP, REDUCE, init) {
  if (!isOpen(X)) {
    open(X)
    on.exit(close(X))
  }
  # Result must be a data frame
  result <- if (missing(init)) {
    data <- YIELD(X)
    if (DONE(data)) return(list())
    MAP(data)
  } else init
  pb <- utils::txtProgressBar(
    min = 0,  max = mx <- nrow(result), style = 3)
  repeat {
    if(nrow(result) == 0) break
    data <- YIELD(X)
    value <- MAP(data)
    result <- REDUCE(result, value)
    setTxtProgressBar(pb, mx - nrow(result))
  }
  close(pb)
}

# Helper function for filter_host_bowtie to write interim fastq file
mk_interim_fastq <- function(reads_bam, read_loc, YS = 10000) {
  message("Writing Intermediate FASTQ file")
  # Pull out all query names
  bf_init <- Rsamtools::BamFile(reads_bam, yieldSize = 100000000)
  allqname <- Rsamtools::scanBam(bf_init, param = Rsamtools::ScanBamParam(
    what = "qname"))[[1]]$qname %>% unique()
  init_df <- dplyr::tibble("Header" = unname(unlist(allqname))) %>%
    dplyr::distinct(.data$Header, .keep_all = TRUE)
  # Define BAM file with smaller yield size
  bf <- Rsamtools::BamFile(reads_bam, yieldSize = YS)
  # Define functions
  YIELD <- function(bf) {
    to_pull <- c("qname", "qual", "seq")
    value <- Rsamtools::scanBam(bf, param = Rsamtools::ScanBamParam(
      what = to_pull))[[1]]
    return(value)
  }
  MAP <- function(value) {
    y <- dplyr::tibble("Header" = value$qname |> unlist(),
                       "Sequence" = as.character(value$seq),
                       "Quality" = as.character(value$qual)) %>%
      dplyr::distinct(.data$Header, .keep_all = TRUE)
    return(y)
  }
  REDUCE <- function(x, y) {
    all_join <- dplyr::left_join(x, y, by = "Header")
    # Print all current reads to file
    all_join %>% dplyr::filter(!is.na(Sequence)) %>%
      dplyr::mutate("Plus" = "+",
                    "Header" = stringr::str_c("@", .data$Header)) %>%
      dplyr::select(.data$Header, .data$Sequence, .data$Plus,
                    .data$Quality) %>% as.matrix() %>% t() %>%
      as.character() %>% dplyr::as_tibble() %>%
      data.table::fwrite(file = paste0(read_loc, ".gz"), compress = "gzip",
                         col.names = FALSE, quote = FALSE,
                         append = TRUE)
    empty_vals <- all_join %>% dplyr::filter(is.na(.data$Sequence)) %>%
      dplyr::select(.data$Header)
    return(empty_vals)
  }
  reduceByYield_iterate(bf, YIELD, MAP, REDUCE, init = init_df)
  message("Intermediate FASTQ file written to", read_loc)
}

#' Helper function to remove reads matched to filter libraries
#'
#' Within the \code{filter_host()} function, we align our sequencing sample to all
#' filter libraries of interest. The \code{remove_matches()} function allows
#' for removal of any target reads that are also aligned to filter libraries.
#'
#' This function is not intended for direct use.
#'
#' @param reads_bam The name of a merged, sorted .bam file that has previously
#' been aligned to a reference library. Likely, the output from running an
#' instance of \code{align_target()}.
#' @param read_names A \code{list} of target query names from \code{reads_bam}
#' that have also aligned to a filter reference library. Each \code{list}
#' element should be a vector of read names.
#' @param output The name of the .bam or .rds file that to which the filtered
#' alignments will be written.
#' @inheritParams filter_host_bowtie
#' @inheritParams filter_host
#' @param threads The number of threads to be used in filtering the bam file.
#' @param aligner The aligner which was used to create the bam file.
#'
#' @return Depending on input \code{make_bam}, either the name of a filtered,
#' sorted .bam file written to the user's current working directory, or
#' an RDS file containing a data frame of only requisite information to run
#' \code{metascope_id()}.
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
#' # remove_matches_bam(readPath, read_names, out, YS = 1000, threads = 1,
#' #                    aligner = "subread", make_bam = FALSE)
#'

remove_matches <- function(reads_bam, read_names, output, YS, threads,
                           aligner, make_bam) {
    message("Removing reads mapped to host indices")
    # some aligned reads may be duplicated; remove these, and unlist
    filter_names <- sort(unique(unlist(read_names)))
    if (make_bam) {
        name_out <- paste0(output, ".bam")
        # obtain vector of target query names from .bam file
        target_reads <- Rsamtools::scanBam(reads_bam)[[1]]$qname
        # define logical vector of which reads to keep, based on query names
        filter_which <- !(target_reads %in% filter_names)
        # index target bam file (reads_BAM already sorted by chromosome)
        bam_index <- Rsamtools::indexBam(reads_bam)
        # create BamFile instance to set yieldSize
        bf <- Rsamtools::BamFile(reads_bam, yieldSize = length(filter_which))
        filtered_bam <- Rsamtools::filterBam(
            bf, destination = name_out, index = bam_index,
            indexDestination = FALSE, filter = filter_which,
            param = Rsamtools::ScanBamParam(what = "qname"))
        file.remove(bam_index)
    } else if (!make_bam) {
        name_out <- paste0(output, ".rds")
        bf <- Rsamtools::BamFile(reads_bam, yieldSize = YS)
        YIELD <- function(x) {
            to_pull <- c("qname", "rname", "cigar", "qwidth", "pos")
            if (identical(aligner, "bowtie")) {
                Rsamtools::scanBam(x, param = Rsamtools::ScanBamParam(
                    what = to_pull, tag = c("AS")))[[1]] |> as.data.frame()
            } else if (identical(aligner, "subread")) {
                Rsamtools::scanBam(x, param = Rsamtools::ScanBamParam(
                    what = to_pull, tag = c("NM")))[[1]] |> as.data.frame()
            }
        }
        MAP <- function(value) value[!(value$qname %in% filter_names), ]
        REDUCE <- function(x, y, iterate = TRUE) dplyr::bind_rows(x, y)
        DONE <- function(value) nrow(value) == 0
        BiocParallel::register(BiocParallel::MulticoreParam(threads))
        suppressWarnings(GenomicFiles::reduceByYield(
            bf, YIELD, MAP, REDUCE, DONE, parallel = TRUE)) %>%
            saveRDS(.data, name_out)
    }
    message("DONE! Alignments written to ", name_out)
    return(name_out)
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
#' @param YS yieldSize, an integer. The number of alignments to be read in from
#' the bam file at once for the creation of an intermediate fastq file.
#' Default is 1000000.
#' @inheritParams filter_host_bowtie
#' @inheritParams align_target
#' @return The name of a filtered, sorted .bam file written to the user's
#' current working directory. Or, if \code{make_bam = FALSE}, an RDS file
#' containing a data frame of only requisite information to run
#' \code{metascope_id()}.
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
#' filter_host(readPath, libs = "filter", threads = 1)
#'

filter_host <- function(reads_bam, libs, lib_dir = NULL, make_bam = FALSE,
                        output = paste(tools::file_path_sans_ext(reads_bam),
                            "filtered", sep = "."), settings = align_details,
                        YS = 500000, threads = 8) {
    # Initialize list of names
    read_names <- vector(mode = "list", length(libs))
    for (i in seq_along(libs)) {
        # Create output file name for BAM
        lib_file <- paste(tools::file_path_sans_ext(reads_bam), ".",
                          libs[i], ".bam", sep = "")
        # Align BAM to the lib & generate new file
        Rsubread::align(index = paste(lib_dir, libs[i], sep = ""),
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
        # Extract target query names from mapped BAM file
        read_names[[i]] <- Rsamtools::scanBam(
            lib_file, param = Rsamtools::ScanBamParam(
                what = c("qname"),
                flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE)))
        # Throw away BAM, vcf file
        file.remove(lib_file)
        file.remove(paste(lib_file, ".indel.vcf", sep = ""))
        file.remove(paste(lib_file, ".summary", sep = ""))
    }
    # helper function to obtain final host-filtered read names
    # remove intermediate fastq file
    name_out <- remove_matches(reads_bam, read_names, output, YS, threads,
                               "subread", make_bam)
    return(name_out)
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
#' Alternatively, an RDS data frame can be output for a smaller output file
#' that is created more efficiently (through parallelization) and is still
#' compatible with \code{metascope_id()}.
#'
#' @param reads_bam The name of a merged, sorted .bam file that has previously
#' been aligned to a reference library. Likely, the output from running an
#' instance of \code{align_target_bowtie()}.
#' @param lib_dir Path to the directory that contains the filter Bowtie2 index
#' files.
#' @param libs The basename of the filter libraries
#' (without .bt2 or .bt2l extension)
#' @param make_bam Logical, whether to also output a bam file with host reads
#' filtered out. An rds file will be created instead if \code{FALSE}.
#' Default is \code{FALSE}.
#' @param output The desired name of the output .bam or .rds file. Extension
#' is automatically defined by whether \code{make_bam = TRUE.} Default is
#' the basename of \code{unfiltered_bam} + \code{.filtered} + extension.
#' @param bowtie2_options Optional: Additional parameters that can be passed to
#' the filter_host_bowtie() function. To see all the available parameters
#' use Rbowtie2::bowtie2_usage(). Default parameters are the parameters are the
#' default parameters that PathoScope 2.0 uses. NOTE: Users should pass all their
#' parameters as one string and if optional parameters are given then the user
#' is responsible for entering all the parameters to be used by Bowtie2. NOTE:
#' The only parameters that should NOT be specified here is the threads.
#' @param YS_1 yieldSize, an integer. The number of alignments to be read in from
#' the bam file at once for the creation of an intermediate fastq file.
#' Default is 1000000.
#' @param YS_2 yieldSize, am integer. The number of alignments to be read in from
#' the bam file at once for the creation of a filtered bam file. Smaller
#' chunks are generally needed for this step, which is why it is a better idea to
#' keep `YS_2` smaller than `YS_1` to conserve memory.
#' Default is 100000.
#' @param threads The amount of threads available for the function.
#' Default is 8 threads.
#' @param overwrite Whether existing files should be overwritten.
#' Default is FALSE.
#'
#' @return The name of a filtered, sorted .bam file written to the user's
#' current working directory. Or, if \code{make_bam = FALSE}, an RDS file
#' containing a data frame of only requisite information to run
#' \code{metascope_id()}.
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
#'                 overwrite=FALSE)
#'
#' ## Filter reads from the bam file that align to the filter library
#' filter_host_bowtie(reads_bam = bamPath, lib_dir = lib_temp,
#'                    libs = "filter", threads = 1)
#'

filter_host_bowtie <- function(reads_bam, lib_dir, libs, make_bam = FALSE,
                               output = paste(
                                   tools::file_path_sans_ext(reads_bam),
                                   "filtered", sep = "."),
                               bowtie2_options = NULL, YS_1 = 1000000,
                               YS_2 = 100000,
                               threads = 8, overwrite = FALSE) {
    # If user does not specify parameters, specify for them
    if (missing(bowtie2_options)) {
      # Relaxed Parameters for Generous Alignments to Metagenomes
      bowtie2_options <- paste("--local -k 100 -D 20 -R 3 -L 3",
                               "-N 1 -p 8 --gbar 1 --mp 3",
                               "--threads", threads)
    } else bowtie2_options <- paste(bowtie2_options, "--threads", threads)
    # Convert reads_bam into fastq (parallelized)
    read_loc <- file.path(dirname(reads_bam), "intermediate.fastq")
    mk_interim_fastq(reads_bam, read_loc, YS_1)
    # Init list of names
    read_names <- vector(mode = "list", length(libs))
    for (i in seq_along(libs)) {
        # Create output file name for BAM
        lib_file <- paste(tools::file_path_sans_ext(reads_bam), ".",
                          libs[i], ".bam", sep = "")
        message("Attempting to perform Bowtie2 alignment on ", libs[i],
                " index")
        # Align reads to lib and generate new filter BAM file
        Rbowtie2::bowtie2_samtools(
            bt2Index = file.path(lib_dir, libs[i]),
            output = tools::file_path_sans_ext(lib_file),
            outputType = "bam", seq1 = read_loc, ... = bowtie2_options,
            overwrite = overwrite)
        # Extract target query names from mapped BAM file
        read_names[[i]] <- Rsamtools::scanBam(
            lib_file, param = Rsamtools::ScanBamParam(
                what = c("qname"),
                flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE)))
        # Throw away BAM file
        file.remove(lib_file)
    }
    # remove intermediate fastq file
    file.remove(read_loc)
    name_out <- remove_matches(reads_bam, read_names, output, YS_2, threads,
                               "bowtie", make_bam)
    return(name_out)
}
