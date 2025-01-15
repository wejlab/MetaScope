#' Filter unmapped reads
#'
#' This function will remove all unmapped reads or lines in a .bam file
#' (warning: overwrites the original file!). This function is needed because
#' combining multiple .bam files from different microbial libraries may lead to
#' some reads that mapped to one library and have unmapped entries from another
#' library. This will remove any unmapped entries and leave all reference mapped
#' lines in the .bam file.
#'
#' It is not intended for direct use.
#'
#' @param bamfile Location for the .bam file to filter & remove all unmapped
#'   reads
#'
#' @return This function will overwrite the existing .bam file with a new .bam
#' file in the same location that has only mapped lines. The function itself
#' returns the output .bam file name.
#'

filter_unmapped_reads <- function(bamfile) {
  sorted_bamfile <- Rsamtools::sortBam(bamfile,
                                       paste(tools::file_path_sans_ext(bamfile),
                                             ".sorted", sep = ""))
  bam_index <- Rsamtools::indexBam(sorted_bamfile)
  filtered_bam <- Rsamtools::filterBam(
    sorted_bamfile,
    destination = bamfile,
    index = bam_index,
    indexDestination = FALSE,
    param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(
      isUnmappedQuery = FALSE))
  )
  # clean up
  file.remove(sorted_bamfile)
  file.remove(bam_index)
  # return filtered file name
  return(filtered_bam)
}

#' Create a combined .bam header
#'
#' This function generates a combined header from multiple .bam files from
#' different reference libraries (e.g. a split bacterial library). It is not
#' intended for use by users.
#'
#' @param bam_files A character vector of the locations/file names of .bam files
#'   from which to combine the headers.
#' @param header_file A file name and location for the output file for the
#'   combined header. This will be a .sam format file without any reads.
#'   Defaults to 'header_tmp.sam'.
#'
#' @return This function will return a combined header from all the supplied
#' .bam files.
#'

combined_header <-
  function(bam_files, header_file = "header_tmp.sam") {
    # get first and last line of header
    bam_head <- Rsamtools::scanBamHeader(bam_files[1])
    n <- length(bam_head[[1]]$text)
    last <- c(names(bam_head[[1]]$text)[n], bam_head[[1]]$text[[n]])
    # open and print the first line of header
    head_con <- file(header_file, open = "w")
    writeLines(c(names(bam_head[[1]]$text)[1], bam_head[[1]]$text[[1]]),
               con = head_con,
               sep = "\t")
    writeLines("\n", con = head_con, sep = "")
    # print genomes from all .bam files
    for (bfile in bam_files) {
      bam_head <- Rsamtools::scanBamHeader(bfile)
      for (j in 2:(length(bam_head[[1]]$text) - 1)) {
        writeLines(c(names(bam_head[[1]]$text)[j], bam_head[[1]]$text[[j]]),
                   con = head_con,
                   sep = "\t")
        writeLines("\n", con = head_con, sep = "")
      }
    }
    writeLines(last, con = head_con, sep = "\t")
    close(head_con)
    return(header_file)
  }

#' Replace the header from a .bam file
#'
#' This function replaces the header from one .bam file with a header from a
#' different .sam file. This function mimics the function of the 'reheader'
#' function in samtools. It is not intended for use by users.
#'
#' @param head A file name and location for the .sam file with the new header.
#' @param old_bam A file name and location for the .bam file which you would
#' @param new_bam A file name for the new .bam file with a replaced header.
#'   Defaults to the same base filename plus 'h.bam'. For example, 'example.bam'
#'   will be written as 'exampleh.bam'.
#'
#' @return This function will return a new .bam file with a replaced header. The
#' function also outputs the new .bam filename.
#'

bam_reheader_R <- function(head, old_bam, new_bam = paste(
  tools::file_path_sans_ext(old_bam), "h.bam", sep = "")) {
  new_sam <- paste(tools::file_path_sans_ext(new_bam), ".sam", sep = "")
  new_sam_con <- file(new_sam, open = "w")
  head_con <- file(head, open = "r")
  while (length(oneLine <- readLines(head_con, n = 1, warn = FALSE)) > 0) {
    writeLines(oneLine, new_sam_con)
  }
  close(head_con)
  old_sam <- Rsamtools::asSam(old_bam, overwrite = TRUE)
  old_sam_con <- file(old_sam, open = "r")
  while (length(oneLine <- readLines(old_sam_con, n = 1, warn = FALSE)) > 0) {
    if (substr(oneLine, 1, 1) != "@") {
      writeLines(oneLine, new_sam_con)
    }
  }
  close(new_sam_con)
  close(old_sam_con)
  file.remove(old_sam)
  new_bam <- Rsamtools::asBam(new_sam, overwrite = TRUE)
  file.remove(new_sam)
  file.remove(paste(new_bam, ".bai", sep = ""))
  return(new_bam)
}

#' Merge multiple .bam files
#'
#' This function merges .bam files. It first used the combined_header function
#' to generate a combined header for all the files, reheaders the files, and
#' then merges and sorts the .bam files. It is similar to the 'samtools merge'
#' function, but it allows the .bam files to have different headers. It is not
#' intended for direct use.
#'
#' @param bam_files A list of file names for the .bam files to be merged.
#' @param destination A file name and location for the merged .bam file.
#' @param head_file A file name and location for the combined header file.
#'   Defaults to the destination. For example, 'example.bam' will be written as
#'   'example.bam'.
#' @param quiet Turns off most messages. Default is \code{TRUE}.
#'
#' @return This function merges .bam files and combines them into a single file.
#' The function also outputs the new .bam filename.
#'

merge_bam_files <- function(bam_files, destination,
                            head_file = paste(destination, "_header.sam",
                                              sep = ""), quiet = TRUE) {
  if (!quiet) message("Merging .bam files")
  com_head <- combined_header(bam_files, header_file = head_file)
  # Paths for merged bam files
  unsortbam <- paste0(destination, "_unsorted.bam")
  if (check_samtools_exists()) {
    if (!quiet) message("samtools found on system. Merging with samtools.")
    lapply(bam_files, function(bf) sys::exec_wait(
      "samtools", c("reheader", head_file, bf), std_out = FALSE,
      std_err = TRUE))
    sys::exec_wait("samtools", c("merge", unsortbam, bam_files),
                   std_out = FALSE, std_err = TRUE)
    lapply(c(bam_files, head_file), unlink, force = TRUE, recursive = TRUE)
    if (!quiet) message("Sorting merged bam file")
    merged_bam_sorted <- paste0(destination, ".bam")
    sys::exec_wait("samtools", c("sort", "-o", merged_bam_sorted, unsortbam),
      std_out = TRUE, std_err = TRUE)
    lapply(unsortbam, unlink, force = TRUE, recursive = TRUE)
  } else {
    if (!quiet) message("samtools not found. Merging with Rsamtools instead.")
    bam_files_h <- NULL
    for (i in seq_along(bam_files)) {
      if (!quiet) message("Reheading bam file")
      new_bam_h <- bam_reheader_R(com_head, bam_files[i])
      bam_files_h <- c(bam_files_h, new_bam_h)
      file.remove(bam_files[i])
      # remove .bam and .vcf and .bam.summary files for each alignment
      file.remove(paste(bam_files[i], ".indel.vcf", sep = ""))
      file.remove(paste(bam_files[i], ".summary", sep = ""))
    }
    merged_bam <- Rsamtools::mergeBam(bam_files_h, unsortbam, overwrite = TRUE)
    # clean up
    file.remove(com_head)
    lapply(bam_files_h, file.remove)
    if (!quiet) message("Sorting merged bam file")
    # sort merged bam file
    merged_bam_sorted <- Rsamtools::sortBam(merged_bam, destination)
    file.remove(merged_bam)
  }
  # return merged and sorted bam
  return(merged_bam_sorted)
}

#' Align microbiome reads to a set of reference libraries
#'
#' This is the main MetaScope target library mapping function, using Rsubread
#' and multiple libraries. Aligns to each library separately, filters unmapped
#' reads from each file, and then merges and sorts the .bam files from each
#' library into one output file. If desired, output can be passed to
#' `filter_host()` to remove reads that also map to filter library genomes.
#'
#' @inheritParams align_target_bowtie
#' @param lib_dir Path to the index files for all libraries.
#' @param libs A vector of character strings giving the basenames of the Subread
#'   index files for alignment. If ALL indices to be used are located in the
#'   current working directory, set \code{lib_dir = NULL}. Default is
#'   \code{lib_dir = NULL}.
#' @param subread_options A named \code{list} specifying alignment parameters
#'   for the \code{Rsubread::align()} function, which is called inside
#'   \code{align_target()}. Elements should include type, nthreads,
#'   maxMismatches, nsubreads, phredOffset, unique, and nBestLocations.
#'   Descriptions of these parameters are available under
#'   \code{?Rsubread::align}. Defaults to the global \code{align_details}
#'   object.
#'
#' @return This function writes a merged and sorted .bam file after aligning to
#' all reference libraries given, along with a summary report file, to the
#' user's working directory. The function also outputs the new .bam filename.
#'
#' @export
#'
#' @examples
#' #### Align example reads to an example reference library using Rsubread
#' \donttest{
#' ## Create temporary directory
#' target_ref_temp <- tempfile()
#' dir.create(target_ref_temp)
#'
#' tax <- "Ovine atadenovirus D"
#'
#' #' ## Create temporary taxonomizr accession
#' namesText<-c(
#' "130327\t|\tOvine atadenovirus D\t|\t|\tscientific name\t|",
#' "100953\t|\tAtadenovirus\t|\t|\tscientific name\t|",
#' "10508\t|\tAdenoviridae\t|\t|\tscientific name\t|",
#' "2732559\t|\tRowavirales\t|\t|\tscientific name\t|",
#' "2732529\t|\tTectiliviricetes\t|\t|\tscientific name\t|",
#' "2732008\t|\tPreplasmiviricota\t|\t|\tscientific name\t|",
#' "2732005\t|\tBamfordvirae\t|\t|\tscientific name\t|",
#' "2732004\t|\tVaridnaviria\t|\t|\tscientific name\t|",
#' "10239\t|\tViruses\t|\t|\tscientific name\t|")
#'
#' nodesText<-c(
#'   "130327\t|\t100953\t|\tspecies",
#'   "100953\t|\t10508\t|\tgenus",
#'   "10508\t|\t2732559\t|\tfamily",
#'   "2732559\t|\t2732529\t|\torder",
#'   "2732529\t|\t2732008\t|\tclass",
#'   "2732008\t|\t2732005\t|\tphylum",
#'   "2732005\t|\t2732004\t|\tkingdom",
#'   "2732004\t|\t10239\t|\tclade",
#'   "10239\t|\t1\t|\tsuperkingdom",
#'   "1\t|\t1\t|\tno rank")
#'
#' tmp_accession<-tempfile()
#' taxonomizr::read.names.sql(textConnection(namesText),tmp_accession, overwrite = TRUE)
#' taxonomizr::read.nodes.sql(textConnection(nodesText),tmp_accession, overwrite = TRUE)
#'
#'
#' ## Download genome
#' all_ref <- MetaScope::download_refseq(tax,
#'                                       reference = FALSE,
#'                                       representative = FALSE,
#'                                       compress = TRUE,
#'                                       out_dir = target_ref_temp,
#'                                       caching = TRUE,
#'                                       accession_path = tmp_accession)
#'
#' ## Create subread index
#' ind_out <- mk_subread_index(all_ref)
#'
#' ## Get path to example reads
#' readPath <- system.file("extdata", "reads.fastq",
#'                         package = "MetaScope")
#' ## Copy the example reads to the temp directory
#' refPath <- file.path(target_ref_temp, "reads.fastq")
#' file.copy(from = readPath, to = refPath)
#'
#' ## Modify alignment parameters object
#' data("align_details")
#' align_details[["type"]] <- "rna"
#' align_details[["phredOffset"]] <- 50
#' # Just to get it to align - toy example!
#' align_details[["maxMismatches"]] <- 100
#'
#' ## Run alignment
#' target_map <- align_target(refPath,
#'                            libs = stringr::str_replace_all(tax, " ", "_"),
#'                            lib_dir = target_ref_temp,
#'                            subread_options = align_details)
#'
#' ## Remove temporary folder
#' unlink(target_ref_temp, recursive = TRUE)
#' unlink(tmp_accession, recursive = TRUE)
#' }

align_target <- function(read1, read2 = NULL, lib_dir = NULL, libs,
                         threads = 1,
                         align_file = tools::file_path_sans_ext(read1),
                         subread_options = align_details, quiet = TRUE) {
  data_env <- new.env(parent = emptyenv())
  utils::data("align_details", envir = data_env, package = "MetaScope")
  align_details <- data_env[["align_details"]]
  if (!is.null(lib_dir)) lib_dir <- tools::file_path_as_absolute(lib_dir)
  bam_files <- numeric(length(libs))
  for (i in seq_along(libs)) {
    if (!quiet) message("Attempting to perform subread alignment on ", libs[i],
            " index")
    bam_files[i] <- paste(tools::file_path_sans_ext(read1), ".", libs[i],
                          ".bam", sep = "")
    Rsubread::align(
      index = file.path(lib_dir, libs[i]),
      readfile1 = read1,
      readfile2 = read2,
      input_format = "fastq",
      output_format = "BAM",
      output_file = bam_files[i],
      type = subread_options[["type"]],
      nthreads = threads,
      maxMismatches = subread_options[["maxMismatches"]],
      nsubreads = subread_options[["nsubreads"]],
      phredOffset = subread_options[["phredOffset"]],
      unique = subread_options[["unique"]],
      nBestLocations = subread_options[["nBestLocations"]]
    )
    ## remove umapped reads
    if (!quiet) message("Filtering unmapped reads")
    filter_unmapped_reads(bam_files[i])
  }
  if (!quiet) message("Library alignment complete")
  # If more than one library was aligned, then combine bam files
  if (length(bam_files) > 1) {
    if (quiet) message("Merging the bam files into ", align_file, ".bam")
    merge_bam_files(bam_files, align_file, quiet=quiet)
  } else {
    file.rename(bam_files, paste(align_file, ".bam", sep = ""))
    # remove Rsubread .vcf and .bam.summary files for now
    file.remove(paste(bam_files, ".indel.vcf", sep = ""))
    file.remove(paste(bam_files, ".summary", sep = ""))
  }
  if (!quiet) message("DONE! Alignments written to ", align_file, ".bam")
  return(paste(align_file, ".bam", sep = ""))
}

#' Align microbiome reads to set of indexed Bowtie2 libraries
#'
#' This is the main MetaScope target library mapping function, using Rbowtie2
#' and multiple libraries. Aligns to each library separately, filters unmapped
#' reads from each file, and then merges and sorts the .bam files from each
#' library into one output file. If desired, output can be passed to
#' `filter_host_bowtie()` to remove reads that also map to filter library
#' genomes.
#'
#' The default parameters are the same that PathoScope 2.0 uses.
#' "--very-sensitive-local -k 100 --score-min L,20,1.0"
#'
#' If you experience any issues with reading the input files, make sure that
#' the file(s) are not located in a read-only folder. This can be circumvented
#' by copying files to a new location before running the function.
#'
#' @param read1 Path to the .fastq file to align.
#' @param read2 Optional: Location of the mate pair .fastq file to align.
#' @param lib_dir Path to the directory that contains the Bowtie2 indexes.
#' @param libs The basename of the Bowtie2 indexes to align against (without
#'   trailing .bt2 or .bt2l extensions).
#' @param align_dir Path to the directory where the output alignment file should
#'   be created.
#' @param align_file The basename of the output alignment file (without trailing
#'   .bam extension).
#' @param  bowtie2_options Optional: Additional parameters that can be passed to
#'   the align_target_bowtie() function. To see all the available parameters use
#'   Rbowtie2::bowtie2_usage(). See Details for default parameters. NOTE: Users
#'   should pass all their parameters as one string and if optional parameters
#'   are given then the user is responsible for entering all the parameters to
#'   be used by Bowtie2. The only parameter that should NOT be specified here is
#'   the number of threads.
#' @param threads The number of threads that can be utilized by the function.
#'   Default is 1 thread.
#' @param overwrite Whether existing files should be overwritten. Default is
#'   FALSE.
#' @param quiet Turns off most messages. Default is \code{TRUE}.
#'
#' @return Returns the path to where the output alignment file is stored.
#'
#' @export
#'
#' @examples
#' #### Align example reads to an example reference library using Rbowtie2
#' \donttest{
#' ## Create temporary directory to store file
#' target_ref_temp <- tempfile()
#' dir.create(target_ref_temp)
#'
#' namesText<-c(
#'   "3052345\t|\tMorbillivirus hominis\t|\t|\tscientific name\t|",
#'   "11229\t|\tMorbillivirus\t|\t|\tscientific name\t|",
#'   "2560076\t|\tOrthoparamyxovirinae\t|\t|\tscientific name\t|",
#'   "11158\t|\tParamyxoviridae\t|\t|\tscientific name\t|",
#'   "11157\t|\tMononegavirales\t|\t|\tscientific name\t|",
#'   "2497574\t|\tMonjiviricetes\t|\t|\tscientific name\t|",
#'   "2497570\t|\tHaploviricotina\t|\t|\tscientific name\t|",
#'   "2497569\t|\tNegarnaviricota\t|\t|\tscientific name\t|",
#'   "2732396\t|\tOrthornavirae\t|\t|\tscientific name\t|",
#'   "2559587\t|\tRiboviria\t|\t|\tscientific name\t|",
#'   "10239\t|\tViruses\t|\t|\tscientific name\t|"
#' )
#'
#'
#' nodesText<-c(
#'   "1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|",
#'   "11234\t|\t3052345\t|\tno rank",
#'   "3052345\t|\t11229\t|\tspecies",
#'   "11229\t|\t2560076\t|\tgenus",
#'   "2560076\t|\t11158\t|\tsubfamily",
#'   "11158\t|\t11157\t|\tfamily",
#'   "11157\t|\t2497574\t|\torder",
#'   "2497574\t|\t2497570\t|\tclass",
#'   "2497570\t|\t2497569\t|\tsubphylum",
#'   "2497569\t|\t2732396\t|\tphylum",
#'   "2732396\t|\t2559587\t|\tkingdom",
#'   "2559587\t|\t10239\t|\tclade",
#'   "10239\t|\t1\t|\tsuperkingdom"
#' )
#'
#' tmp_accession<-tempfile()
#' taxonomizr::read.names.sql(textConnection(namesText),tmp_accession, overwrite = TRUE)
#' taxonomizr::read.nodes.sql(textConnection(nodesText),tmp_accession, overwrite = TRUE)
#'
#' ## Dowload reference genome
#' MetaScope::download_refseq("Morbillivirus hominis",
#'                            reference = FALSE,
#'                            representative = FALSE,
#'                            compress = TRUE,
#'                            out_dir = target_ref_temp,
#'                            caching = TRUE,
#'                            accession_path = tmp_accession
#' )
#'
#' ## Create temporary directory to store the indices
#' index_temp <- tempfile()
#' dir.create(index_temp)
#'
#' ## Create bowtie2 index
#' MetaScope::mk_bowtie_index(
#'   ref_dir = target_ref_temp,
#'   lib_dir = index_temp,
#'   lib_name = "target",
#'   overwrite = TRUE
#' )
#'
#' ## Create temporary directory for final file
#' output_temp <- tempfile()
#' dir.create(output_temp)
#'
#' ## Get path to example reads
#' readPath <- system.file("extdata", "virus_example.fastq",
#'                         package = "MetaScope")
#'
#' ## Align to target genomes
#' target_map <-
#'   MetaScope::align_target_bowtie(
#'     read1 = readPath,
#'     lib_dir = index_temp,
#'     libs = "target",
#'     align_dir = output_temp,
#'     align_file = "bowtie_target",
#'     overwrite = TRUE,
#'     bowtie2_options = "--very-sensitive-local"
#'   )
#'
#' ## Remove extra folders
#' unlink(target_ref_temp, recursive = TRUE)
#' unlink(index_temp, recursive = TRUE)
#' unlink(output_temp, recursive = TRUE)
#' unlink(tmp_accession, recursive = TRUE)
#' }

align_target_bowtie <- function(read1, read2 = NULL, lib_dir, libs,
                                align_dir, align_file, bowtie2_options = NULL,
                                threads = 1, overwrite = FALSE, quiet = TRUE) {
    # Convert user specified paths to absolute paths for debugging purposes
    lib_dir <- tools::file_path_as_absolute(lib_dir)
    align_dir <- tools::file_path_as_absolute(align_dir)
    # If user does not specify parameters, specify for them
    if (is.null(bowtie2_options)) {
      bowtie2_options <- paste(
        "--local -R 2 -N 0 -L 25 -i S,1,0.75 -k 5 --score-min L,0,1.7",
        "--threads", threads)
    } else bowtie2_options <- paste(bowtie2_options, "--threads", threads)
    bam_files <- numeric(length(libs))
    for (i in seq_along(libs)) {
      # Don't attach .bam extension - Rbowtie2 does this already
      bam_files[i] <- file.path(
        align_dir, paste(basename(tools::file_path_sans_ext(read1)),
                         ".", libs[i], sep = ""))
      if (!quiet) message("Attempting to perform Bowtie2 alignment on ",
                          libs[i], " index")
      Rbowtie2::bowtie2_samtools(
        bt2Index = file.path(lib_dir, libs[i]),
        output = bam_files[i],
        outputType = "bam",
        seq1 = read1,
        seq2 = read2,
        overwrite = overwrite,
        ... = bowtie2_options
      )
      unlink(".bowtie2.cerr.txt")
      # Attach .bam extension to bam files to call function
      if (!quiet) message("Filtering unmapped reads")
      filter_unmapped_reads(paste0(bam_files[i], ".bam"))
    }
    if (!quiet) message("Library alignment complete")
    outputFile <- file.path(align_dir, paste0(align_file, ".bam"))
    bam_files <- paste0(bam_files, ".bam")
    # If more than one library was aligned, then combine bam files
    if (length(bam_files) > 1) {
      if (!quiet) message("Merging the bam files into ", align_file, ".bam")
      merge_bam_files(bam_files, tools::file_path_sans_ext(outputFile), quiet=quiet)
    } else file.rename(bam_files, outputFile)
    if (!quiet) message("DONE! Alignments written to ", outputFile)
    return(outputFile)
  }
