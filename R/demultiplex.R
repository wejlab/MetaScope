#' Helper function for demultiplexing
#'
#' Helper function for demultiplexing sequencing reads, designed in a way to
#' allow for parallelization across barcodes (parallel extraction of reads by
#' barcode). This function takes a specific barcode (numeric index) from lists
#' of sample names/barcodes, a \code{Biostrings::DNAStringSet} of barcodes by
#' sequence header, and a \code{Biostrings::QualityScaledXStringSet} of reads
#' corresponding to the barcodes. Based on the barcode index given, it extracts
#' all reads for the indexed barcode and writes all the reads from that barcode
#' to a separate .fastq file.
#'
#' @param barcodeIndex Which barcode (integer number or index) in the barcodes
#'   or sample name to use for read extraction.
#' @param barcodes A list of all barcodes in the sequencing dataset. Correlates
#'   and in same order as \code{sampleNames}.
#' @param sampleNames A list of sample names or identifiers associated with each
#'   barcode in the barcodes list.
#' @param index A \code{Biostrings::DNAStringSet} that contains the read headers
#'   and barcode sequence for each header in the sequence slot.
#' @param reads A \code{Biostrings::QualityScaledXStringSet} that has the same
#'   headers and order as the index file, but contains the read sequences and
#'   their quality scores.
#' @param location A directory location to store the demultiplexed read files.
#'   Defaults to generate a new subdirectory at './demultiplex_fastq'
#' @param rcBarcodes Should the barcode indices in the barcodes list be reverse
#'   complemented to match the sequences in the index DNAStringSet? Defaults to
#'   \code{TRUE}.
#' @param hDist Uses a Hamming Distance or number of base differences to allow
#'   for inexact matches for the barcodes/indexes. Defaults to 0. Warning: if
#'   the Hamming Distance is >=1 and this leads to inexact index matches to more
#'   than one barcode, that read will be written to more than one demultiplexed
#'   read files.
#' @param quiet Turns off most messages. Default is \code{TRUE}.
#'
#' @return Writes a single .fastq file that contains all reads whose index
#'   matches the barcode specified. This file will be written to the location
#'   directory, and will be named based on the specified sampleName and barcode,
#'   e.g. './demultiplex_fastq/SampleName1_GGAATTATCGGT.fastq.gz'
#' @export
#' @examples
#'
#' ## Create temporary directory
#' ref_temp <- tempfile()
#' dir.create(ref_temp)
#'
#' ## Load example barcode, index, and read data into R session
#' barcodePath <- system.file("extdata", "barcodes.txt", package = "MetaScope")
#' bcFile <- read.table(barcodePath, sep = "\t", header = TRUE)
#'
#' indexPath <- system.file("extdata", "virus_example_index.fastq",
#' package = "MetaScope")
#' inds <- Biostrings::readDNAStringSet(indexPath, format = "fastq")
#'
#' readPath <- system.file("extdata", "virus_example.fastq",
#'                         package = "MetaScope")
#' reads <- Biostrings::readQualityScaledDNAStringSet(readPath)
#'
#' ## Extract reads from the first barcode
#' results <- extract_reads(1, bcFile[, 2], bcFile[, 1], inds, reads,
#'                         rcBarcodes = FALSE, location = ref_temp)
#'
#' ## Extract reads from multiple barcodes
#' more_results <- lapply(1:6, extract_reads, bcFile[, 2], bcFile[, 1], inds,
#'                        reads, rcBarcodes = FALSE, location = ref_temp)
#'
#' ## Remove temporary directory
#' unlink(ref_temp, recursive = TRUE)
#'

extract_reads <- function(barcodeIndex, barcodes, sampleNames, index, reads,
                          location = "./demultiplex_fastq", rcBarcodes = TRUE,
                          hDist = 0, quiet = TRUE) {
  barcode <- barcodes[barcodeIndex]
  sampleName <- sampleNames[barcodeIndex]
  if (!quiet) message("Finding reads for barcode: ", barcode)
  if (rcBarcodes) {
    rci <- as.character(Biostrings::reverseComplement(
      Biostrings::DNAString(barcode)))
  } else {
    rci <- barcode
  }
  ind_match <- utils::adist(as.character(index), rci) <= hDist
  numReads <- sum(ind_match)
  outFileName <- paste(location, "/", sampleName, "_", barcode,
                       ".fastq.gz", sep = "")
  if (numReads == 0 && !quiet) {
    message("\tFound 0 reads for this barcode; no file will be written.")
  } else {
    if (!quiet) message("\tFound ", sum(ind_match),
                        " reads, writing reads to: ", outFileName)
    Biostrings::writeQualityScaledXStringSet(
      reads[c(ind_match)], outFileName, compress = TRUE)
  }
  return(list(output_file = outFileName, numberOfReads = numReads,
              matchedIndexes = ind_match))
}

errmessages <- function(barcodes, samNames, numReads) {
  warning("Did not find any reads for the following barcodes: ",
          paste(barcodes[numReads == 0], collapse = " "))
  warning("Did not find any reads for the following samples: ",
          paste(samNames[numReads == 0], collapse = " "))
  write(paste("Did not find any reads for the following barcodes:",
              paste(barcodes[numReads == 0], collapse = " "), "\n",
              "Did not find any reads for the following samples: ",
              paste(samNames[numReads == 0], collapse = " ")),
        file = "demultiplex_fastq/unmapped_barcodes_samples.txt")
}

#' Demultiplexing sequencing reads
#'
#' Function for demultiplexing sequencing reads arranged in a common format
#' provided by sequencers (such as Illumina) generally for 16S data. This
#' function takes a matrix of sample names/barcodes, a .fastq file of barcodes
#' by sequence header, and a .fastq file of reads corresponding to the barcodes.
#' Based on the barcodes given, the function extracts all reads for the indexed
#' barcode and writes all the reads from that barcode to separate .fastq files.
#' @param barcodeFile Path to a file containing a .tsv matrix with a
#'   header row, and then sample names (column 1) and barcodes (column 2).
#' @param indexFile Path to a .fastq file that contains the barcodes for
#'   each read. The headers should be the same (and in the same order) as
#'   \code{readFile}, and the sequence in the \code{indexFile} should be the
#'   corresponding barcode for each read. Quality scores are not considered.
#' @param readFile Path to the sequencing read .fastq file that corresponds
#'   to the \code{indexFile}.
#' @param rcBarcodes Should the barcode indexes in the barcodeFile be reverse
#'   complemented to match the sequences in the \code{indexFile}? Defaults to
#'   \code{TRUE}.
#' @param location A directory location to store the demultiplexed read files.
#'   Defaults to generate a new temporary directory.
#' @param threads The number of threads to use for parallelization
#'   (BiocParallel). This function will parallelize over the barcodes and
#'   extract reads for each barcode separately and write them to separate
#'   demultiplexed files.
#' @param hammingDist Uses a Hamming Distance or number of base differences to
#'   allow for inexact matches for the barcodes/indexes. Defaults to \code{0}.
#'   Warning: if the Hamming Distance is \code{>=1} and this leads to inexact
#'   index matches to more than one barcode, that read will be written to more
#'   than one demultiplexed read files.
#' @param quiet Turns off most messages. Default is \code{TRUE}.
#'
#' @return Returns multiple .fastq files that contain all reads whose index
#'   matches the barcodes given. These files will be written to the location
#'   directory, and will be named based on the given sampleNames and barcodes,
#'   e.g. './demultiplex_fastq/SampleName1_GGAATTATCGGT.fastq.gz'
#'
#' @export
#'
#' @examples
#'
#' ## Get barcode, index, and read data locations
#' barcodePath <- system.file("extdata", "barcodes.txt", package = "MetaScope")
#' indexPath <- system.file("extdata", "virus_example_index.fastq",
#'                          package = "MetaScope")
#' readPath <- system.file("extdata", "virus_example.fastq",
#'                          package = "MetaScope")
#'
#' ## Demultiplex
#' demult <- meta_demultiplex(barcodePath, indexPath, readPath, rcBarcodes = FALSE,
#'                       hammingDist = 2)
#' demult
#'

meta_demultiplex <- function(barcodeFile, indexFile, readFile, rcBarcodes = TRUE,
                             location = NULL, threads = 1,
                             hammingDist = 0, quiet = TRUE) {
  if (is.null(location)) {
    location <- tempfile()
    dir.create(location)
  }
  if (!quiet) message("Reading Sample Names and Barcodes from: ",
                      barcodeFile)
  bcFile <- utils::read.table(barcodeFile, sep = "\t", header = TRUE)
  barcodes <- bcFile[, 2]
  samNames <- bcFile[, 1]
  if (!quiet) message("\tFound information for ", length(barcodes),
                      " samples/barcodes")
  if (!quiet) message("Reading Index File: ", indexFile)
  inds <- Biostrings::readDNAStringSet(indexFile, format = "fastq")
  if (!quiet) message("\tFound indexes for ", length(inds), " reads")
  if (!quiet) message("Reading Sequence File: ", readFile)
  reads <- Biostrings::readQualityScaledDNAStringSet(readFile)
  if (!quiet) message("\tFound ", length(reads), " reads")
  ## make output directory if necessary
  if (!dir.exists(location)) dir.create(location)
  # Loop over barcodes
  numReads <- NULL
  ind_no_match <- numeric(length(reads))
  for (i in seq_along(barcodes)) {
    extracted <- extract_reads(
      i, barcodes, samNames, inds, reads, rcBarcodes = rcBarcodes,
      location = location, hDist = hammingDist, quiet = quiet)
    numReads <- c(numReads, extracted$numberOfReads)
    ind_no_match <- ind_no_match + extracted$matchedIndexes
  }
  if (!quiet) message(sum(ind_no_match > 1))
  ind_no_match <- (ind_no_match == 0)
  # number of reads for each barcode
  if (any(numReads == 0)) {
    errmessages(barcodes, samNames, numReads)
  }
  # Track reads without matches, and write them to an 'orphan' file
  if (!quiet) message("Found ", sum(ind_no_match),
                      " reads without a matching barcode (",
                      100 * round(mean(ind_no_match), 4), "%), writing reads to: ",
                      location, "/orphans.fastq.gz")
  Biostrings::writeQualityScaledXStringSet(
    reads[c(ind_no_match)], paste(location, "/orphans.fastq.gz",
                                  sep = ""), compress = TRUE)
  summaryMat <- cbind(bcFile[seq_along(barcodes), ],
                      NumberOfReads = numReads)
  utils::write.table(summaryMat,
                     file = paste(location, "/summary.txt", sep = ""),
                     col.names = FALSE, row.names = TRUE, quote = FALSE)
  return(summaryMat)
}
