#' Count the number of base lengths in a CIGAR string for a given operation
#' 
#' The ‘CIGAR’ (Compact Idiosyncratic Gapped Alignment Report) string is how
#' the SAM/BAM format represents spliced alignments. This function will accept
#' a CIGAR string for a single read and a single character indicating the
#' operation to be parsed in the string. An operation is a type of column that
#' appears in the alignment, e.g. a match or gap. The integer following the
#' operator specifies a number of consecutive operations. The `count_matches()`
#' will identify all occurrences of the operator in the string input, add them,
#' and return an integer number representing the total number of operations
#' for the read that was summarized by the input CIGAR string.
#' 
#' This function is best used on a vector of CIGAR strings using an apply
#' function (see examples).
#' 
#' @param x Character. A CIGAR string for a read to be parsed. Examples of
#' possible operators include "M", "D", "I", "S", "H", "=", "P", and "X".
#' @param char A single letter representing the operation to total for the
#' given string.
#' 
#' @return an integer number representing the total number of alignment
#' operations for the read that was summarized by the input CIGAR string.
#' 
#' @export
#' 
#' @examples 
#' # A single cigar string: 3M + 3M + 5M
#' cigar1 <- "3M1I3M1D5M"
#' count_matches(cigar1, char = "M")
#' 
#' # Parse with operator "P": 2P
#' cigar2 <- "4M1I2P9M"
#' count_matches(cigar2, char = "P")
#' 
#' # Apply to multiple strings: 1I + 1I + 5I
#' cigar3 <- c("3M1I3M1D5M", "4M1I1P9M", "76M13M5I")
#' sapply(cigar3, count_matches, char = "I")
#'

count_matches <- function(x, char = "M") {
    if (length(char) != 1) {
        stop("Please provide a single character operator with which to parse.")
    } else if (length(x) != 1) {
        stop("Please provide a single CIGAR string to be parsed.")
    }
    pattern <- paste("\\d+", char , sep = "")
    ind <- gregexpr(pattern, x)[[1]]
    start <- as.numeric(ind)
    end <- start + attr(ind, "match.length") - 2
    out <- sum(as.numeric(apply(cbind(start, end), 1,
                                function(y) substr(x, start = y[1], stop = y[2]))))
    return(data.table::fifelse(is.na(out[1]), yes = 0, no = out[1]))
}

#' Determines if a read is multi-mapping 
#' 
#' A read should be multi-mapping if it appears more than once in a bam file. 
#' The function accepts the read group size which will either be 1 if the read
#' is mapped uniquely or >1 if the read is multi-mapping. The function will 
#' return 0 if the read is unique and 1 if the read is multi-mapping. 
#' The function is intended to be used to create the uniqueness indicator 
#' (y_ind_2) which is used for calculating theta. Not meant to be used outside
#' of function.
#' 
#' @param x Integer. An integer representing the group size of a specific read.
#' @return either 0 (unique) or 1 (multi-mapping)

unique_identifier <- function(x)
{
    if(x == 1){
        return(0)
    }
    else{
        return(1)
    }
}

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
#' @param aligner The aligner which was used to create the bam file. Default is 
#' "bowtie" but can also be set to "subread" or "other"
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
#' mapped reads. The function itself returns the output .csv file name.
#'
#' @export
#'
#' @examples
#' 
#' ## Create object with path to example reference fasta file
#' refPath <- system.file("extdata","Mononegavirales.fasta", package = "MetaScope")
#' 
#' ## Copy the example reference fasta file to the current directory
#' file.copy(from = refPath, to = file.path(".", "Mononegavirales.fasta"))
#'
#' ## Make and align to a single a reference genome library
#' mk_subread_index('Mononegavirales.fasta')
#' readPath <- system.file("extdata", "virus_example.fastq", package = "MetaScope")
#' viral_map <- align_target(readPath, "Mononegavirales", project_name = "virus_example")
#'
#' #### Apply MetaScope ID:
#' metascope_id(viral_map, aligner="subread")
#' 


metascope_id <- function(bam_file, aligner = "bowtie", out_file = paste(tools::file_path_sans_ext(bam_file),".metascope_id.csv", sep = ""), EMconv = 1/10000, EMmaxIts = 25) {
    
    # Check to make sure valid aligner is specified
    if (aligner != "bowtie" && aligner != "subread" && aligner != "other")
        stop("Please make sure aligner is set to either 'bowtie', 'subread', or 'other'")
    
    message("Reading .bam file: ", bam_file)
    
    # If the aligner is set to bowtie then also extract the AS (alignment score) tag
    if (aligner == "bowtie")
        params <- Rsamtools::ScanBamParam(what = c("qname", "rname", "cigar", "qwidth"), tag = c("AS"))
    # If the aligner is set to subread then also extract the NM (edit score) tag
    else if (aligner == "subread")
        params <- Rsamtools::ScanBamParam(what = c("qname", "rname", "cigar", "qwidth"), tag = c("NM"))
    # If the aligner is set to other then we do not extract any additional tags
    else if (aligner == "other")
        params <- Rsamtools::ScanBamParam(what = c("qname", "rname", "cigar", "qwidth"))
    
    reads <- Rsamtools::scanBam(bam_file, param = params)
    unmapped <- is.na(reads[[1]]$rname)
    mapped_qname <- reads[[1]]$qname[!unmapped]
    mapped_rname <- reads[[1]]$rname[!unmapped]
    
    # If aligner is set to bowtie then we have to do additional processing to get accession numbers
    if (aligner == "bowtie")
        mapped_rname <- gsub(".*accession\\|","",mapped_rname)
    
    mapped_cigar <- reads[[1]]$cigar[!unmapped]
    mapped_qwidth <- reads[[1]]$qwidth[!unmapped]
    
    
    # If aligner is set to bowtie then we want to get the mapped alignment scores
    if (aligner == "bowtie")
        mapped_alignment <- reads[[1]][["tag"]][["AS"]][!unmapped]
    
    # If aligner is set to subread then we want to get the mapped edit scores
    else if (aligner == "subread")
        mapped_edit <- reads[[1]][["tag"]][["NM"]][!unmapped]
    
    
    read_names <- unique(mapped_qname)
    accessions <- unique(mapped_rname)
    
    message("\tFound ", length(read_names), " reads aligned to ",
            length(accessions), " NCBI accessions")
    
    # Convert accessions to taxids and get genome names
    message("Obtaining taxonomy and genome names")
    suppressMessages(tax_id_all <- taxize::genbank2uid(id = accessions))
    taxids <- vapply(tax_id_all, function(x) x[1], character(1))
    unique_taxids <- unique(taxids)
    taxid_inds <- match(taxids, unique_taxids)
    genome_names <- vapply(tax_id_all, function(x) attr(x, "name"), character(1))
    unique_genome_names <- genome_names[!duplicated(taxid_inds)]
    message("\tFound ", length(unique_taxids), " unique NCBI taxonomy IDs")
    
    # Make an aligment matrix (rows: reads, cols: unique taxids)
    message("Setting up the EM algorithm")
    qname_inds <- match(mapped_qname, read_names)
    rname_inds <- match(mapped_rname, accessions)
    rname_tax_inds <- taxid_inds[rname_inds]
    # Don't uncomment below line 
    #cigar_strings <- mapped_cigar[rname_inds]
    qwidths <- mapped_qwidth[rname_inds]
    
    # Order based on read names
    rname_tax_inds <- rname_tax_inds[order(qname_inds)]
    # Modified this from original
    cigar_strings <- mapped_cigar[order(qname_inds)]
    
    # If aligner is set to bowtie we want to order the alignment scores
    if (aligner == "bowtie")
        scores <- mapped_alignment[order(qname_inds)]
    # If aligner is set to subread we want to order the edit scores
    else if (aligner == "subread")
        scores <- mapped_edit[order(qname_inds)]
    # If aligner is set to other then we make no assumptions about scores
    else if (aligner == "other")
        scores <- 1
    
    qwidths <- qwidths[order(qname_inds)]
    qname_inds <- sort(qname_inds)
    
    
    
    # If aligner is subread then to get an alignment score we subtract the edit score
    # from the number of nucleotide matches which is taken from the cigar string.
    if (aligner == "subread"){
        num_match <- unlist(vapply(cigar_strings, count_matches, USE.NAMES = FALSE, double(1)))
        alignment_score <- num_match - scores
        relative_alignment_score <- alignment_score - min(alignment_score)
        exp_alignment_score <- 2^relative_alignment_score
    }
    
    # If aligner is bowtie then we already have an alignment score 
    else if (aligner == "bowtie"){
        relative_alignment_score <- scores - min(scores)
        exp_alignment_score <- 2^relative_alignment_score
    }
    
    # If aligner is other then we make no assumptions about the alignment score
    else if (aligner == "other"){
        exp_alignment_score <- 1
    }
    
    
    combined <- dplyr::bind_cols("qname" = qname_inds, "rname" = rname_tax_inds, "scores" = exp_alignment_score)
    
    input_distinct <- dplyr::distinct(combined, qname, rname, .keep_all = TRUE)
    qname_inds_2 <- input_distinct$qname
    rname_tax_inds_2 <- input_distinct$rname
    
    # Normalize the exponential alignment scores so that the scores 
    # represent a probability that the read originates from that genome
    by_read <- dplyr::group_by(input_distinct, qname)
    scores_2 <- dplyr::summarize(by_read, scores_2 = scores/(sum(scores)))$scores_2
    
    # Uniqueness indicator vector (1 = multimapping, 0 = unique or non multimapping)
    y_ind_2 <- dplyr::summarize(by_read, multimapping_2 = unique_identifier(dplyr::n()))$multimapping_2
    
    
    gammas <- Matrix::sparseMatrix(qname_inds_2, rname_tax_inds_2, x = scores_2)
    
    pi_old <- rep(1 / nrow(gammas), ncol(gammas))
    pi_new <-  Matrix::colMeans(gammas)
    
    # Added this
    theta_new <- Matrix::colMeans(gammas)
    
    
    conv <- max(abs(pi_new - pi_old) / pi_old)
    it <- 0
    
    message("Starting EM iterations")
    while (conv > EMconv & it < EMmaxIts) {
        # Expectation Step: Estimate expected value for each read to each genome
        pi_mat <- Matrix::sparseMatrix(qname_inds_2, rname_tax_inds_2, x = pi_new[rname_tax_inds_2])
        
        # Added this
        theta_mat <- Matrix::sparseMatrix(qname_inds_2, rname_tax_inds_2, x = theta_new[rname_tax_inds_2])
        
        # Modified this
        weighted_gamma <- gammas * pi_mat * theta_mat
        
        weighted_gamma_sums <- Matrix::rowSums(weighted_gamma)
        gammas_new <- weighted_gamma/weighted_gamma_sums
        
        # Maximization step: proportion of reads to each genome
        pi_new <- Matrix::colMeans(gammas_new)
        
        #Added this
        theta_new <- (Matrix::colSums(y_ind_2*gammas_new)+1) / (nrow(gammas_new)+1)
        
        
        # Check convergence
        it <- it + 1
        conv <- max(abs(pi_new - pi_old) / pi_old, na.rm = TRUE)
        pi_old <- pi_new
        print(c(it, conv))
    }
    message("\tDONE! Converged in ", it, " interations.")
    
    # Collect results
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
    
    aubs <- results <- cbind(TaxonomyID = final_taxids, Genome = final_genomes, read_count = best_hit, Proportion = proportion, EMreads = EMreads, EMProportion = EMprop)
    results <- results[order(best_hit, decreasing = TRUE), ]
    message("Found reads for ", length(best_hit), " genomes")
    
    # Write to file
    utils::write.csv(results, file = out_file, row.names = FALSE)
    message("Results written to ", out_file)
    return(results)
}
