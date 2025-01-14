# Helper function for remove_matches
reduceByYield_RM <- function(X, YIELD, MAP, DONE, filter_names, num_reads) {
  if (!Rsamtools::isOpen(X)) {
    Rsamtools::open.BamFile(X)
    on.exit(Rsamtools::close.BamFile(X))
  }
  pb <- utils::txtProgressBar(min = 0,  max = num_reads, style = 3)
  k <- 0
  repeat {
    if (DONE(data <- YIELD(X))) break
    k <- k + nrow(data)
    MAP(data, filter_names)
    utils::setTxtProgressBar(pb, k)
  }
  base::close(pb)
}

# Helper function for mk_interim_fastq()
reduceByYield_iterate <- function(X, YIELD, MAP, REDUCE, init) {
  if (!Rsamtools::isOpen(X)) {
    Rsamtools::open.BamFile(X)
    on.exit(Rsamtools::close.BamFile(X))
  }
  # Result must be a data frame
  result <- init
  pb <- utils::txtProgressBar(min = 0,
                              max = mx <- nrow(result),
                              style = 3)
  repeat {
    if (nrow(result) == 0)
      break
    data <- YIELD(X)
    value <- MAP(data)
    result <- REDUCE(result, value)
    utils::setTxtProgressBar(pb, mx - nrow(result))
  }
  close(pb)
}

# Helper function for filter_host_bowtie to write interim fastq file
mk_interim_fastq <- function(reads_bam, read_loc, YS, quiet) {
  unlink(read_loc, recursive = FALSE)
  if (!quiet) message("Writing Intermediate FASTQ file to ", read_loc)
  # Pull out all query names
  bf_init <- Rsamtools::BamFile(reads_bam, yieldSize = 100000000)
  allqname <-
    Rsamtools::scanBam(bf_init, param = Rsamtools::ScanBamParam(
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
    all_join %>% dplyr::filter(!is.na(.data$Sequence)) %>%
      dplyr::mutate("Plus" = "+",
                    "Header" = stringr::str_c("@", .data$Header)) %>%
      dplyr::select(.data$Header, .data$Sequence, .data$Plus,
                    .data$Quality) %>% as.matrix() %>% t() %>%
      as.character() %>% dplyr::as_tibble() %>%
      data.table::fwrite(
        file = read_loc, compress = "gzip", col.names = FALSE,
        quote = FALSE, append = TRUE)
    empty_vals <- all_join %>% dplyr::filter(is.na(.data$Sequence)) %>%
      dplyr::select(.data$Header)
    return(empty_vals)
  }
  reduceByYield_iterate(bf, YIELD, MAP, REDUCE, init = init_df)
  if (!quiet) message("Intermediate FASTQ file written to", read_loc)
}

#' Helper function to remove reads matched to filter libraries
#'
#' Using the \code{filter_host()} function, we align our sequencing sample to
#' all filter libraries of interest. The \code{remove_matches()} function allows
#' for removal of any target reads that are also aligned to filter libraries.
#'
#' This function is not intended for direct use.
#'
#' @param reads_bam The name of a merged, sorted .bam file that has previously
#'   been aligned to a reference library. Likely, the output from running an
#'   instance of \code{align_target()}.
#' @param read_names A \code{list} of target query names from \code{reads_bam}
#'   that have also aligned to a filter reference library. Each \code{list}
#'   element should be a vector of read names.
#' @param output The name of the .bam or .csv.gz file that to which the filtered
#'   alignments will be written.
#' @inheritParams filter_host_bowtie
#' @inheritParams filter_host
#' @param threads The number of threads to be used in filtering the bam file.
#'   Default is 1.
#' @param aligner The aligner which was used to create the bam file.
#'
#' @return Depending on input \code{make_bam}, either the name of a filtered,
#'   sorted .bam file written to the user's current working directory, or an RDS
#'   file containing a data frame of only requisite information to run
#'   \code{metascope_id()}.
#'

remove_matches <- function(reads_bam, read_names, output, YS, threads,
                           aligner, make_bam, quiet) {
  if (!quiet) message("Removing reads mapped to host indices")
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
    Rsamtools::filterBam(
      bf, destination = name_out, index = bam_index,
      indexDestination = FALSE, filter = filter_which,
      param = Rsamtools::ScanBamParam(what = "qname"))
    # file.remove(bam_index) keep bam index
  } else if (!make_bam) {
    name_out <- paste0(output, ".csv.gz")
    numread <- Rsamtools::BamFile(reads_bam, yieldSize = 100000000) %>%
      Rsamtools::scanBam(param = Rsamtools::ScanBamParam(
        what = "pos")) %>% magrittr::extract2(1) %>% unlist() %>% length()
    # Define BAM file with smaller yield size
    bf <- Rsamtools::BamFile(reads_bam, yieldSize = YS)
    YIELD <- function(X) {
      to_pull <- c("qname", "rname", "cigar", "qwidth", "pos")
      if (identical(aligner, "bowtie2")) {
        out <- Rsamtools::scanBam(X, param = Rsamtools::ScanBamParam(
          what = to_pull, tag = c("AS")))[[1]]
        out[to_pull] %>% dplyr::as_tibble() %>%
          dplyr::mutate(tag = out$tag$AS) %>% return()
      } else if (identical(aligner, "subread")) {
        out <- Rsamtools::scanBam(X, param = Rsamtools::ScanBamParam(
          what = to_pull, tag = c("NM")))[[1]]
        out[to_pull] %>% dplyr::as_tibble() %>%
          dplyr::mutate(tag = out$tag$NM) %>% return()
      }
    }
    MAP <- function(value, filter_names) {
      value %>%
        dplyr::filter(!(.data$qname %in% filter_names)) %>%
        data.table::fwrite(file = name_out, compress = "gzip",
                           col.names = FALSE, quote = TRUE, append = TRUE,
                           nThread = threads)
      return("")
    }
    DONE <- function(data) nrow(data) == 0
    reduceByYield_RM(bf, YIELD, MAP, DONE, filter_names, numread)
  }
  if (!quiet) message("DONE! Alignments written to ", name_out)
  return(name_out)
}

#' Use Rsubread to align reads against one or more filter libraries and
#' subsequently remove mapped reads
#'
#' After aligning your sample to a target library with \code{align_target()},
#' use \code{filter_host()} to remove unwelcome host contamination using filter
#' reference libraries. This function takes as input the name of the .bam file
#' produced via \code{align_target()}, and produces a sorted .bam file with any
#' reads that match the filter libraries removed. This resulting .bam file may
#' be used upstream for further analysis. This function uses Rsubread. For the
#' Rbowtie2 equivalent of this function, see \code{filter_host_bowtie}.
#'
#' A compressed .csv can be created to produce a smaller output file that is
#' created more efficiently and is still compatible with \code{metascope_id()}.
#'
#' @param reads_bam The name of a merged, sorted .bam file that has previously
#'   been aligned to a reference library. Likely, the output from running an
#'   instance of \code{align_target()}.
#' @param lib_dir Path to the directory that contains the filter Subread index
#'   files.
#' @param libs The basename of the filter libraries (without index
#'   extension).
#' @inheritParams filter_host_bowtie
#' @inheritParams align_target
#' @return The name of a filtered, sorted .bam file written to the user's
#'   current working directory. Or, if \code{make_bam = FALSE}, a .csv.gz file
#'   containing a data frame of only requisite information to run
#'   \code{metascope_id()}.
#'
#' @export
#'
#' @examples
#' #### Filter reads from bam file that align to any of the filter libraries
#'
#' ## Assuming a bam file has been created previously with align_target()
#' \donttest{
#' ## Create temporary directory
#' filter_ref_temp <- tempfile()
#' dir.create(filter_ref_temp)
#'
#'#' # Create temporary taxonomizr accession
#' namesText<-c(
#'   "1\t|\troot\t|\t\t|\tscientific name\t|",
#'   "2\t|\tBacteria\t|\tBacteria <prokaryotes>\t|\tscientific name\t|",
#'   "418127\t|\tStaphylococcus aureus subsp. aureus Mu3\t|\t|\tscientific name\t|",
#'   "273036\t|\tStaphylococcus aureus RF122\t|\t|\tscientific name\t|",
#'   "426430\t|\tStaphylococcus aureus subsp. aureus str. Newman\t|\t|\tscientific name\t|",
#'   "1280\t|\tStaphylococcus aureus\t|\t|\tscientific name\t|",
#'   "1279\t|\tStaphylococcus\t|\t|\tscientific name\t|",
#'   "90964\t|\tStaphylococcaceae\t|\t|\tscientific name\t|",
#'   "1385\t|\tBacillales\t|\t|\tscientific name\t|",
#'   "91061\t|\tBacilli\t|\t|\tscientific name\t|",
#'   "1239\t|\tBacillota\t|\t|\tscientific name\t|",
#'   "1783272\t|\tBacillati\t|\t|\tscientific name\t|"
#' )
#'
#' nodesText<-c(
#'   "1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|",
#'   "2\t|\t131567\t|\tsuperkingdom\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|",
#'   "418127\t|\t1280\t|\tstrain",
#'   "273036\t|\t1280\t|\tstrain",
#'   "426430\t|\t1280\t|\tstrain",
#'   "1280\t|\t1279\t|\tspecies",
#'   "1279\t|\t90964\t|\tgenus",
#'   "90964\t|\t1385\t|\tfamily",
#'   "1385\t|\t91061\t|\torder",
#'   "91061\t|\t1239\t|\tclass",
#'   "1239\t|\t1783272\t|\tphylum",
#'   "1783272\t|\t2\t|\tkingdom",
#'   "131567\t|\t1\t|\tno rank"
#' )
#'
#' tmp_accession<-tempfile()
#' taxonomizr::read.names.sql(textConnection(namesText),tmp_accession, overwrite = TRUE)
#' taxonomizr::read.nodes.sql(textConnection(nodesText),tmp_accession, overwrite = TRUE)
#'
#'
#' ## Download filter genome
#' all_species <- c("Staphylococcus aureus subsp. aureus str. Newman")
#' all_ref <- vapply(all_species, MetaScope::download_refseq,
#'                   reference = FALSE, representative = FALSE, compress = TRUE,
#'                   out_dir = filter_ref_temp, caching = FALSE,
#'                   accession_path = tmp_accession,
#'                   FUN.VALUE = character(1))
#' ind_out <- vapply(all_ref, mk_subread_index, FUN.VALUE = character(1))
#'
#' ## Get path to example reads
#' readPath <- system.file("extdata", "subread_target.bam",
#'                         package = "MetaScope")
#'
#' ## Copy the example reads to the temp directory
#' refPath <- file.path(filter_ref_temp, "subread_target.bam")
#' file.copy(from = readPath, to = refPath)
#' utils::data("align_details")
#' align_details[["type"]] <- "rna"
#' align_details[["phredOffset"]] <- 10
#' # Just to get it to align - toy example!
#' align_details[["maxMismatches"]] <- 10
#'
#' ## Align and filter reads
#' filtered_map <- filter_host(
#'   refPath, lib_dir = filter_ref_temp,
#'   libs = stringr::str_replace_all(all_species, " ", "_"),
#'   threads = 1, subread_options = align_details)
#'
#' ## Remove temporary directory
#' unlink(filter_ref_temp, recursive = TRUE)
#' unlink(tmp_accession, recursive = TRUE)
#' }
#'

filter_host <- function(reads_bam, lib_dir = NULL, libs, make_bam = FALSE,
                        output = paste(tools::file_path_sans_ext(reads_bam),
                                       "filtered", sep = "."),
                        subread_options = align_details, YS = 100000,
                        threads = 1, quiet = TRUE) {
  data_env <- new.env(parent = emptyenv())
  utils::data("align_details", envir = data_env, package = "MetaScope")
  align_details <- data_env[["align_details"]]
  # Convert reads_bam into fastq
  read_loc <- file.path(dirname(output), "intermediate.fastq.gz")
  mk_interim_fastq(reads_bam, read_loc, YS, quiet = quiet)
  # Align
  read_names <- vector(mode = "list", length(libs))
  for (i in seq_along(libs)) {
    # Create output file name for BAM
    lib_file <- paste(tools::file_path_sans_ext(reads_bam), ".",
                      libs[i], ".bam", sep = "")
    # Align BAM to the lib & generate new file
    Rsubread::align(index = file.path(lib_dir, libs[i]),
                    readfile1 = read_loc,
                    input_format = "fastq",
                    output_file = lib_file,
                    type = subread_options[["type"]],
                    nthreads = threads,
                    maxMismatches = subread_options[["maxMismatches"]],
                    nsubreads = subread_options[["nsubreads"]],
                    phredOffset = subread_options[["phredOffset"]],
                    unique = subread_options[["unique"]],
                    nBestLocations = subread_options[["nBestLocations"]])
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
  # remove intermediate fastq file
  file.remove(read_loc)
  # Filter out host-aligned reads
  name_out <- remove_matches(reads_bam, read_names, output, YS, threads,
                             "subread", make_bam, quiet = quiet)
  return(name_out)
}

#' Use Rbowtie2 to align reads against one or more filter libraries and
#' subsequently remove mapped reads
#'
#' After a sample is aligned to a target library with
#' \code{align_target_bowtie()}, we may use \code{filter_host_bowtie()} to
#' remove unwelcome host contamination using filter reference libraries. This
#' function takes as input the name of the .bam file produced via
#' \code{align_target_bowtie()}, and produces a sorted .bam or .csv.gz file with
#' any reads that match the filter libraries removed. This resulting .bam file
#' may be used downstream for further analysis. This function uses Rbowtie2 For
#' the Rsubread equivalent of this function, see \code{filter_host}.
#'
#' A compressed .csv can be created to produce a smaller output file that is
#' created more efficiently and is still compatible with \code{metascope_id()}.
#'
#' The default parameters are the same that PathoScope 2.0 uses.
#' "--very-sensitive-local -k 100 --score-min L,20,1.0"
#'
#' @param reads_bam The name of a merged, sorted .bam file that has previously
#'   been aligned to a reference library. Likely, the output from running an
#'   instance of \code{align_target_bowtie()}.
#' @param lib_dir Path to the directory that contains the filter Bowtie2 index
#'   files.
#' @param libs The basename of the filter libraries (without .bt2 or .bt2l
#'   extension).
#' @param make_bam Logical, whether to also output a bam file with host reads
#'   filtered out. A .csv.gz file will be created instead if \code{FALSE}.
#'   Creating a bam file is costly on resources over creating a compressed
#'   csv file with only relevant information, so default is \code{FALSE}.
#' @param output The desired name of the output .bam or .csv.gz file. Extension is
#'   automatically defined by whether \code{make_bam = TRUE.} Default is the
#'   basename of \code{unfiltered_bam} + \code{.filtered} + extension.
#' @param bowtie2_options Optional: Additional parameters that can be passed to
#'   the filter_host_bowtie() function. To see all the available parameters use
#'   Rbowtie2::bowtie2_usage(). See Details for default parameters.
#'   NOTE: Users should pass all their parameters as
#'   one string and if optional parameters are given then the user is
#'   responsible for entering all the parameters to be used by Bowtie2.
#'   The only parameters that should NOT be specified here is the threads.
#' @param YS yieldSize, an integer. The number of alignments to be read in from
#'   the bam file at once for chunked functions. Default is 100000.
#' @param threads The amount of threads available for the function. Default is 1
#'   thread.
#' @param overwrite Whether existing files should be overwritten. Default is
#'   FALSE.
#' @param quiet Turns off most messages. Default is \code{TRUE}.
#'
#' @return The name of a filtered, sorted .bam file written to the user's
#'   current working directory. Or, if \code{make_bam = FALSE}, a .csv.gz file
#'   containing a data frame of only requisite information to run
#'   \code{metascope_id()}.
#'
#' @export
#'
#' @examples
#' #### Filter reads from bam file that align to any of the filter libraries
#'
#' ## Assuming a bam file has already been created with align_target_bowtie()
#' # Create temporary filter library
#' filter_ref_temp <- tempfile()
#' dir.create(filter_ref_temp)
#'
#' #' ## Create temporary taxonomizr accession
#' namesText<-c(
#' "3052462\t|\tOrthoebolavirus zairense\t|\t|\tscientific name\t|",
#' "3044781\t|\tOrthoebolavirus\t|\t|\tscientific name\t|",
#' "11266\t|\tFiloviridae\t|\t|\tscientific name\t|",
#' "11157\t|\tMononegavirales\t|\t|\tscientific name\t|",
#' "2497574\t|\tMonjiviricetes\t|\t|\tscientific name\t|",
#' "2497570\t|\tHaploviricotina\t|\t|\tscientific name\t|",
#' "2497569\t|\tNegarnaviricota\t|\t|\tscientific name\t|",
#' "2732396\t|\tOrthornavirae\t|\t|\tscientific name\t|",
#' "2559587\t|\tRiboviria\t|\t|\tscientific name\t|",
#' "10239\t|\tViruses\t|\t|\tscientific name\t|")
#'
#' nodesText<-c(
#'   "3052462\t|\t3044781\t|\tspecies",
#'   "3044781\t|\t11266\t|\tgenus",
#'   "11266\t|\t11157\t|\tfamily",
#'   "11157\t|\t2497574\t|\torder",
#'   "2497574\t|\t2497570\t|\tclass",
#'   "2497570\t|\t2497569\t|\tsubphylum",
#'   "2497569\t|\t2732396\t|\tsphylum",
#'   "2732396\t|\t2559587\t|\tkingdom",
#'   "2559587\t|\t10239\t|\tclade",
#'   "10239\t|\t1\t|\tsuperkingdom",
#'   "1\t|\t1\t|\tno rank")
#'
#' tmp_accession<-tempfile()
#' taxonomizr::read.names.sql(textConnection(namesText),tmp_accession, overwrite = TRUE)
#' taxonomizr::read.nodes.sql(textConnection(nodesText),tmp_accession, overwrite = TRUE)
#'
#'
#' ## Download reference genome
#' MetaScope::download_refseq("Orthoebolavirus zairense",
#'                            reference = FALSE,
#'                            representative = FALSE,
#'                            compress = TRUE,
#'                            out_dir = filter_ref_temp,
#'                            caching = TRUE,
#'                            accession_path = tmp_accession)
#'
#' ## Create temp directory to store the indices
#' index_temp <- tempfile()
#' dir.create(index_temp)
#'
#' ## Create filter index
#' MetaScope::mk_bowtie_index(
#'   ref_dir = filter_ref_temp,
#'   lib_dir = index_temp,
#'   lib_name = "filter",
#'   overwrite = TRUE
#' )
#'
#' ## Create temporary folder to hold final output file
#' output_temp <- tempfile()
#' dir.create(output_temp)
#'
#' ## Get path to example bam
#' bamPath <- system.file("extdata", "bowtie_target.bam",
#'                        package = "MetaScope")
#' target_copied <- file.path(output_temp, "bowtie_target.bam")
#' file.copy(bamPath, target_copied)
#'
#' ## Align and filter reads
#' filter_out <-
#'   filter_host_bowtie(
#'     reads_bam = target_copied,
#'     lib_dir = index_temp,
#'     libs = "filter",
#'     threads = 1
#'   )
#'
#' ## Remove temporary directories
#' unlink(filter_ref_temp, recursive = TRUE)
#' unlink(index_temp, recursive = TRUE)
#' unlink(output_temp, recursive = TRUE)
#' unlink(tmp_accession, recursive = TRUE)

filter_host_bowtie <- function(reads_bam, lib_dir, libs, make_bam = FALSE,
                               output = paste(
                                 tools::file_path_sans_ext(reads_bam),
                                 "filtered", sep = "."),
                               bowtie2_options = NULL, YS = 100000,
                               threads = 1, overwrite = FALSE, quiet = TRUE) {
  # If user does not specify parameters, specify for them
  if (is.null(bowtie2_options)) {
    bowtie2_options <- paste("--local -R 2 -N 0 -L 25 -i S,1,0.75 -k 5 --score-min L,0,1.7",
                             "--threads", threads)
  } else bowtie2_options <- paste(bowtie2_options, "--threads", threads)
  # Convert reads_bam into fastq
  read_loc <- file.path(dirname(output), "intermediate.fastq.gz")
  mk_interim_fastq(reads_bam, read_loc, YS, quiet = quiet)
  # Init list of names
  read_names <- vector(mode = "list", length(libs))
  for (i in seq_along(libs)) {
    # Create output file name for BAM
    lib_file <- paste(tools::file_path_sans_ext(reads_bam), ".",
                      libs[i], ".bam", sep = "")
    if (!quiet) message("Attempting to perform Bowtie2 alignment on ",
                        libs[i], " index")
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
    unlink(".bowtie2.cerr.txt")
  }
  # remove intermediate fastq file
  file.remove(read_loc)
  # Filter out host-aligned reads
  name_out <- remove_matches(reads_bam, read_names, output, YS, threads,
                             "bowtie2", make_bam, quiet = quiet)
  return(name_out)
}
