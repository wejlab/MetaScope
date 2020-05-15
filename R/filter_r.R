#' Align reads against one or more filter libraries and subsequently
#' remove mapped reads
#'
#' Description goes here
#'
#' @param
#'
#' @return
#' 
#' @examples

filter_r <- function() {
  # Input: name of merged BAM file aligned to target library
  # Locate that file in user's working directory

  # Initialize list of names

  # for loop (or helper function?); for each filter lib:
    # Align BAM to the lib
    # Generate new BAM file
    # save (into running list?) header/read names that succesfully hit that filter
    # throw away bam file
  # end loop
  
  # Now use Rsamtools; steal code from align_target??
  # Go into original BAM file and remove anything with same read name
    # To do this - check helper function in align_target

  # output that filtered BAM file

}

##### My notes

# Think: how to get list of headers with which to filter
# Goal: run in 5-10 minutes or less
# Note that output of align_R is already sorted by chromosome
  # Maybe then we sort the BAM files from filtering and do it that way?