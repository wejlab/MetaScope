# MetaScope (dev version)

## Bug Fixes
* Fixed taxonomy table function to output correctly formatted table
* Fixed examples for various functions that were calling genomes with download_refseq but genomes were not able to be found.
* Fixed plot generation for `metascope_id()`
* Fixed premature stopping of `download_refseq` for strains labeled as "no rank" in NCBI.

# MetaScope 0.99.0 (4/19/2022)

* Pre-Release version of MetaScope

## Bug Fixes
* Fixed check error message about `data.table::fread` for reading .gz files by adding `R.utils` to imports.

## Major Changes
* Submitted to Bioconductor

## Minor Changes
