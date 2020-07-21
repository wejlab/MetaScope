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
#' basename plus ".MetaScopeID.csv".
#' @param EMconv The convergence parameter of the EM algorithm. Default set at
#' \code{0.001}.
#' @param EMmaxIts The maximum number of EM iterations, regardless of whether
#' the EMconv is below the threshhold. Default set at \code{50}.
#' If set at \code{0}, the algorithm skips the EM step and summarizes the .bam file 'as is'
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
                                          ".MetaScopeID.csv", sep = ""),
                         EMconv = 1/10000, EMmaxIts = 25) {
  ## read in .bam file
  message("Reading .bam file: ", bam_file)
  reads <- Rsamtools::scanBam(bam_file, 
                              param = Rsamtools::ScanBamParam(what = c("qname",
                                                                       "rname")))
  unmapped <- is.na(reads[[1]]$rname)
  mapped_qname <- reads[[1]]$qname[!unmapped]
  mapped_rname <- reads[[1]]$rname[!unmapped]
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
  tax_query <- data.frame( rname_tax_inds , qname_inds )
  
  #keep only one alignment to each taxid
  tax_query <- tax_query %>% distinct(rname_tax_inds, qname_inds, .keep_all = TRUE)
  
  #split unique and multimapped reads
  tax_query <- tax_query %>% group_by(qname_inds) %>% mutate(duplicate.flag = n() > 1)
  tax_query_unique <- tax_query %>% filter(!duplicate.flag)
  tax_query_mulitmap <- tax_query %>% filter(duplicate.flag)
  
  #summarize results for uniquely mapped reads
  tax_count <- sapply(1:max(tax_query$rname_tax_inds),function(i){sum(tax_query_unique$rname_tax_inds==i)})
  total_unique <- sum(tax_count)
  
  # Functions for updating gamma and pi
  pi_fun <- function(x, gammas, rname , tax_count, total_unique){ ( tax_count[x]+sum(gammas[rname == x])) /( total_unique + sum(gammas) ) }
  gamma_fun <- function(x, qname_inds, pis){  pis[qname_inds == x]/sum(pis[qname_inds == x]) }

  ## EM algorithm for reducing abiguity in the alignments
  gammas <- rep(1, nrow(tax_query_mulitmap))
  pi_old <- 1 / max(tax_query_mulitmap$rname_tax_inds)
  pi_new <- sapply(1:max(tax_query$rname_tax_inds), pi_fun, gammas, tax_query_mulitmap$rname_tax_inds, tax_count, total_unique)
  conv <- max(abs(pi_new - pi_old) )
  it <- 0

  message("Starting EM iterations")
  while (conv > EMconv & it < EMmaxIts) {
    ## Expectation Step: Estimate the expected value for each multi-mapped read to each genome
    gammas <- unlist(sapply(unique(tax_query_mulitmap$qname_inds), gamma_fun, tax_query_mulitmap$qname_inds, pi_new[tax_query_mulitmap$rname_tax_inds]))

    ## Maimization step: proportion of reads to each genome 
    pi_new <- sapply(1:max(tax_query$rname_tax_inds), pi_fun, gammas, tax_query_mulitmap$rname_tax_inds, tax_count, total_unique)
    
    ## Check convergence
    it <- it + 1
    conv <- max(abs(pi_new - pi_old), na.rm = TRUE)
    pi_old <- pi_new
  }
  message("\tDONE! Converged in ", it, " interations.")

  ## Collect results
  gamma_count <- rep(0,length(tax_count)); for(i in unique(tax_query_mulitmap$rname_tax_inds)){gamma_count[i]=sum(gammas[tax_query_mulitmap$rname_tax_inds==i])}
  
  max_fun <- function(i, gam, qname){max_gam <- max(gam[qname ==i]); 1*(gam[qname ==i]==max_gam)/sum(gam[qname ==i]==max_gam)}
  gamma_max <- unlist(sapply(unique(tax_query_mulitmap$qname_inds), max_fun, gammas, tax_query_mulitmap$qname_inds))
  gamma_max_count <- rep(0,length(tax_count)); for(i in unique(tax_query_mulitmap$rname_tax_inds)){gamma_max_count[i]=sum(gamma_max[tax_query_mulitmap$rname_tax_inds==i])}
    
  best_hits <- round(tax_count + gamma_max_count, 1)
  proportion <- best_hits / sum(best_hits)
  
  EM_hits <- round(tax_count + gamma_count, 1)
  EMprop <- EM_hits / sum(EM_hits)
  
  results <- cbind(TaxonomyID = unique_taxids, Genome = unique_genome_names,
                   read_count = best_hits, Proportion = proportion,
                   EMreads = EM_hits,
                   EMProportion = EMprop)
  results <- results[order(best_hits, decreasing = TRUE), ]
  
   message("Found reads for ", length(best_hits), " genomes")

  ## Write to file
  write.csv(results, file = out_file, row.names = F)
  message("Results written to ", out_file)

  return(list(out_file, results))
}
