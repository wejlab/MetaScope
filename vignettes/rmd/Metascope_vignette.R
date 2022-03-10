## -----------------------------------------------------------------------------
readPath <- system.file("extdata", "SRR606366.fastq", package = "MetaScope")

readPath

## ----ref lib, eval = FALSE----------------------------------------------------
#  ## Code block not run
#  
#  download_refseq('viral', compress = FALSE)

## ----ref lib bacterial, eval = FALSE------------------------------------------
#  ## Code block not run
#  
#  # Representative MUST be set to TRUE for the archaea library
#  download_refseq('archaea', compress = FALSE, representative = TRUE)

## ----demultiplex, eval = FALSE------------------------------------------------
#  ## Code block not run
#  
#  # Get barcode, index, and read data locations
#  barcodePath <- system.file("extdata", "barcodes.txt", package = "MetaScope")
#  indexPath <- system.file("extdata", "virus_example_index.fastq",
#                           package = "MetaScope")
#  readPath <- system.file("extdata", "virus_example.fastq", package = "MetaScope")
#  
#  # Get barcode, index, and read data locations
#  demult <- demultiplex(barcodePath, indexPath, readPath, rcBarcodes = FALSE,
#                        hammingDist = 2)

## ----make indexes, eval = FALSE-----------------------------------------------
#  ## Code block not run
#  
#  mk_subread_index('bacteria.fasta')
#  
#  
#  mk_subread_index('archaea.fasta')

## ----alignment align, eval = FALSE--------------------------------------------
#  ## Code block not run
#  
#  bacteria_map <- align_target(readPath, libs = "bacteria",
#                            project_name = "bacteria_example")

## ----alignment filter, eval = FALSE-------------------------------------------
#  ## Code block not run
#  
#  readPath2 <- system.file("extdata", "bacteria_example.bam", package = "MetaScope")
#  final_map <- filter_host(readPath2, libs = "archaea")
#  # Produces file bacteria_example.filtered.bam

## ----identification, eval = FALSE---------------------------------------------
#  MetaScope_id(final_map)

## ----session info-------------------------------------------------------------
sessionInfo()


