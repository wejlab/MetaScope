# MetaScope (dev version)

## Bug Fixes
* Identified k as being doubled when specified for bowtie filter or align steps. Fixed the 16S params object accordingly.

# MetaScope (Bioc 3.18)

## Bug Fixes
* Fixed taxonomy table function to output correctly formatted table
* Fixed examples for various functions that were calling genomes with download_refseq but genomes were not able to be found.
* Fixed plot generation for `metascope_id()`
* Fixed premature stopping of `download_refseq` for strains labeled as "no rank" in NCBI.
* Fixed identification of reads as unknown genomes (due to outdated reference databases that identify genomes now removed from NCBI). The unknown genomes will now be distinctly identified based on NCBI accession ID, which will separate them in the final results. IDs can also be looked up manually (on NCBI website) to see what the reads were aligning to.
* Incorporated unknown taxa (removed from NCBI databases) into output of `convert_animalcules()`

## Major changes
* Added ability to identify reads from databases other than NCBI (known to work for Silva and Greengenes2)

# MetaScope 0.99.0 (4/19/2022)

* Pre-Release version of MetaScope

## Bug Fixes
* Fixed check error message about `data.table::fread` for reading .gz files by adding `R.utils` to imports.

## Major Changes
* Submitted to Bioconductor

## Minor Changes
