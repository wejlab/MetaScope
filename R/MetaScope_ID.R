#' MetaScope ID
#'
#' This function will read in a .bam file, annotate the taxonomy and genome
#' names, reduce the mapping ambiguity using a mixture model, and output a
#' .csv file with the results. Currently, it assumes that the genome
#' library/.bam files use NCBI accession names for reference names (rnames in
#' .bam file). 
#'
#' @param bam_file The .bam file that needs to be summarized, annotated, and
#' needs removal of ambiguity.
#' @param out_file The name of the .csv output file. Defaults to the bam_file
#' basename plus ".metascope_id.csv".
#' @param EMconv The convergence parameter of the EM algorithm. Default set at
#' \code{1/10000}.
#' @param EMmaxIts The maximum number of EM iterations, regardless of whether
#' the EMconv is below the threshhold. Default set at \code{50}.
#' If set at \code{0}, the algorithm skips the EM step and summarizes the .bam
#' file 'as is'
#' 
#' @return
#' This function returns a .csv file with annotated read counts to genomes with
#' mapped reads. The function iself returns the output .csv file name.
#'
#' @examples
#' ## Get a reference genome library
#' download_refseq('viral', compress = FALSE)
#'
#' ## Make and align to a single a reference genome library
#' mk_subread_index('viral.fasta')
#' readPath <- system.file("extdata", "virus_example.fastq",
#' package = "MetaScope")
#' viral_map <- align_target(readPath, "viral", "virus_example")
#'
#' #### Apply MetaScope ID:
#' metascope_id(viral_map)
#'
#' @export

metascope_id <- function(bam_file, 
                         out_file = paste(tools::file_path_sans_ext(bam_file),
                                          ".metascope_id.csv", sep = ""),
                         EMconv = 1/10000, EMmaxIts = 25) {
  set.seed(99)
  now <- Sys.time()
  message("Reading .bam file: ", bam_file)
  reads <- Rsamtools::scanBam(bam_file, 
                              param = Rsamtools::ScanBamParam(what = c("qname",
                                                                       "rname",
                                                                       "cigar")))
  unmapped <- is.na(reads[[1]]$rname)
  mapped_qname <- reads[[1]]$qname[!unmapped]
  mapped_rname <- reads[[1]]$rname[!unmapped]
  mapped_cigar <- reads[[1]]$cigar[!unmapped]
  read_names <- unique(mapped_qname)
  accessions <- unique(mapped_rname)
  message("\tFound ", length(read_names), " reads aligned to ",
          length(accessions), " NCBI accessions")
  
  ## convert accessions to taxids and get genome names
  message("Obtaining taxonomy and genome names")
  suppressMessages(tax_id_all <- taxize::genbank2uid(id = accessions))
  taxids <- sapply(tax_id_all, function(x) x[1])
  unique_taxids <- unique(taxids)
  taxid_inds <- match(taxids, unique_taxids)
  genome_names <- sapply(tax_id_all, function(x) attr(x, "name"))
  unique_genome_names <- genome_names[!duplicated(taxid_inds)]
  message("\tFound ", length(unique_taxids), " unique NCBI taxonomy IDs")
  
  ## make an aligment matrix (rows: reads, cols: unique taxids)
  message("Setting up the EM algorithm")
  qname_inds <- match(mapped_qname, read_names)
  rname_inds <- match(mapped_rname, accessions)
  rname_tax_inds <- taxid_inds[rname_inds]
  
  #order based on read names
  rname_tax_inds <- rname_tax_inds[order(qname_inds)]
  qname_inds <- sort(qname_inds)
  
  
  input_distinct <- distinct(tibble(cbind(qname_inds, rname_tax_inds)))
  qname_inds_2 <- input_distinct$`cbind(qname_inds, rname_tax_inds)`[, 1]
  rname_tax_inds_2 <- input_distinct$`cbind(qname_inds, rname_tax_inds)`[, 2]
  gammas <- Matrix::sparseMatrix(qname_inds_2, rname_tax_inds_2,
                                 x = 1) #align_scores_cigar)
  
  pi_old <- rep(1 / nrow(gammas), ncol(gammas))
  pi_new <-  Matrix::colMeans(gammas)
  conv <- max(abs(pi_new - pi_old) / pi_old)
  it <- 0
  
  message("Starting EM iterations")
  while (conv > EMconv & it < EMmaxIts) {
    # Expectation Step: Estimate expected value for each read to each genome
    pi_mat <- Matrix::sparseMatrix(qname_inds_2, rname_tax_inds_2,
                                   x = pi_new[rname_tax_inds_2])
    weighted_gamma <- gammas * pi_mat
    weighted_gamma_sums <- Matrix::rowSums(weighted_gamma)
    gammas_new <- weighted_gamma/weighted_gamma_sums
    
    # Maximization step: proportion of reads to each genome
    pi_new <- Matrix::colMeans(gammas_new)
    
    # Check convergence
    it <- it + 1
    conv <- max(abs(pi_new - pi_old) / pi_old, na.rm = TRUE)
    pi_old <- pi_new
    print(c(it, conv))
  }
  message("\tDONE! Converged in ", it, " interations.")
  
  ## Collect results
  hit_which <- qlcMatrix::rowMax(gammas_new, which = TRUE)$which
  best_hit <- Matrix::colSums(hit_which)
  names(best_hit) <- seq_along(best_hit)
  best_hit <- best_hit[best_hit != 0]
  
  hits_ind <- as.numeric(names(best_hit))
  
  final_taxids <- unique_taxids[hits_ind]
  final_genomes <- unique_genome_names[hits_ind]
  
  proportion <- best_hit / sum(best_hit)
  gammasums <- Matrix::colSums(gammas_new)
  EMreads <- round(gammasums[hits_ind], 1)
  EMprop <- gammasums[hits_ind] / sum(gammas_new)
  
  aubs <- results <- cbind(TaxonomyID = final_taxids, Genome = final_genomes,
                           read_count = best_hit, Proportion = proportion,
                           EMreads = EMreads,
                           EMProportion = EMprop)
  results <- results[order(best_hit, decreasing = TRUE), ]
  message("Found reads for ", length(best_hit), " genomes")
  
  ## Write to file
  write.csv(results, file = out_file, row.names = FALSE)
  message("Results written to ", out_file)

  return(list(out_file, results))
}
