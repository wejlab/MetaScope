The example reads provided in this package are simulated reads obtained through the wgsim function of SAMtools (https://doi.org/10.1093/bioinformatics/btp352). Using wgsim, 1000 simulated reads of Staphylococcus aureus strain RF122 were obtained. In addition, 500 simulated reads of Staphylococcus epidermidis strain RP62A were also obtained to act as a filter target.

Below are the steps to recreating the reads.

## Inside R    

library(MetaScope)
download_refseq(taxon = "Staphylococcus aureus RF122", representative = FALSE, reference = FALSE, compress = FALSE)
download_refseq(taxon = "Staphylococcus epidermidis RP62A", representative = FALSE, reference = FALSE, compress = FALSE)

## Outside R

mv Staphylococcus\ aureus\ RF122.fasta Staphylococcus_aureus_RF122.fasta
mv Staphylococcus\ epidermidis\ RP62A.fasta Staphylococcus_epidermidis_RP62A.fasta
wgsim -N 1000 -1 250 -2 250 -S 32 Staphylococcus_aureus_RF122.fasta read1.fq read2.fq 
wgsim -N 500 -1 250 -2 250 -S 32 Bacillus_cereus_AH1237.fasta read3.fq read4.fq  

cat read1.fq read3.fq > reads.fastq

Files created: reads.fastq

The target fasta file is created by merging the following strains into one fasta file: Staphylococcus aureus strain RF122, Staphylococcus aureus strain ST398, Staphylococcus aureus strain N315, Staphylococcus aureus strain Mu50, Staphylococcus aureus strain Mu3, and Staphylococcus aureus strain Newman

Below are the steps to recreate the target fasta file.

## Inside R

library(MetaScope)
download_refseq(taxon = "Staphylococcus aureus RF122", reference = FALSE, representative = FALSE, compress = FALSE)
download_refseq(taxon = "Staphylococcus aureus subsp. aureus ST398", reference = FALSE, representative = FALSE, compress = FALSE)
download_refseq(taxon = "Staphylococcus aureus subsp. aureus N315", reference = FALSE, representative = FALSE, compress = FALSE)
download_refseq(taxon = "Staphylococcus aureus subsp. aureus Mu50", reference = FALSE, representative = FALSE, compress = FALSE)
download_refseq(taxon = "Staphylococcus aureus subsp. aureus Mu3", reference = FALSE, representative = FALSE, compress = FALSE)
download_refseq(taxon = "Staphylococcus aureus subsp. aureus str. Newman", reference = FALSE, representative = FALSE, compress = FALSE)

## Outside R (assuming directory contains only the downloaded fasta files)

cat *.fasta > target.fasta

Files created: target.fasta

The filter fasta file is the fasta file of Bacillus cereus strain AH1273. 

Below are the steps to recreate the filter fasta file.

## Inside R

download_refseq(taxon = "Staphylococcus epidermidis RP62A", representative = FALSE, reference = FALSE, compress = FALSE)

## Outside R

mv Staphylococcus\ epidermidis\ RP62A.fasta filter.fasta

Files created: filter.fasta










