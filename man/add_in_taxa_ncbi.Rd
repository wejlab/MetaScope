% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metascope_blast.R
\name{add_in_taxa_ncbi}
\alias{add_in_taxa_ncbi}
\title{Adds in taxa if NCBI database
Returns MetaScope Table with NCBI taxa in separate columns}
\usage{
add_in_taxa_ncbi(combined_pre, NCBI_key)
}
\arguments{
\item{combined_pre}{MetaScope ID file with NCBI taxa qnames}

\item{NCBI_key}{(character) NCBI Entrez API key. optional.
See taxize::use_entrez(). Due to the high number of requests made to NCBI,
this function will be less prone to errors if you obtain an NCBI key.}
}
\description{
Adds in taxa if NCBI database
Returns MetaScope Table with NCBI taxa in separate columns
}
