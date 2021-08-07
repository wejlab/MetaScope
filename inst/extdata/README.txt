#### Data Generation Steps

## Creating demultiplexed_virus_example.fastq
 
barcodePath <- system.file("extdata", "barcodes.txt", package = "MetaScope")

indexPath <- system.file("extdata", "virus_example_index.fastq",
                         package = "MetaScope")

readPath <- system.file("extdata", "virus_example.fastq",
                         package = "MetaScope")
                
demult <- demultiplex(barcodePath, indexPath, readPath, rcBarcodes = FALSE,
                      hammingDist = 2)

cd demultiplex_fastq

Uncompress the *.fastq.gz files

cat *.fastq > demultiplexed_virus_example.fastq

## Creating Mononegavirales.fasta file

download_refseq(taxon = "Mononegavirales", representative = FALSE, reference = FALSE)

## Creating the Morbillivirus.fasta file

download_refseq(taxon = "Morbillivirus", representative = FALSE, reference = FALSE)







