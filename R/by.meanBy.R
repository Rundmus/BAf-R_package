# -----------------------------------------------------------------------------#
#' Mean within subset
#' 
#' Calculate an average of a target within the subset of samples divided by the
#' specified column in \code{@sinfo}
#' 
#' @param baf an object of the \code{\link{BAf-class}}
#' @param by the column in \code{@sinfo} by which samples will be stratified 
#' @param geometric if the geometic mean should be calculated 
#' 
#' @return 
#' A \code{"matrix"} of which the columns are levels in the columns of 'by' and
#' the rows are targets
#' 
#' @examples
#' data(sba)
#' meanBy(sba, "plate")		# average by "plate"
#' meanBy(sba, "dis", geometric = TRUE)		# geometric mean by "dis"
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @export
# -----------------------------------------------------------------------------#
# created  : 2011-07-13 by Mun-Gwan
# modified : 
#   2012-05-04 by Mun-Gwan : generalize from only for plate to any column in 
#     @sample
#   2013-07-03 by Mun-Gwan : adapt to "SBAe" class
# -----------------------------------------------------------------------------#

meanBy <- function(baf, by, geometric = FALSE) {
  stopifnot(inherits(baf, "BAf"))
  stopifnot(by %in% colnames(baf@sinfo))
  
  byCol <- baf@sinfo[, by] %>% unlist %>% factor		# the column of 'by'
  
  if(geometric) {
    # take care of log(0)
    sX(baf)[sX(baf) <= 0] <- (min(sX(baf)[sX(baf) > 0]) - 1) %>% max(1)
    
    sapply(levels(byCol), function(j) {
      sX(baf)[byCol == j, ] %>% log() %>% apply(2, mean, na.rm = TRUE) %>% exp()
    })
  } else {
    sapply(levels(byCol), function(j) {
      sX(baf)[byCol == j, ] %>% apply(2, mean, na.rm = TRUE)
    })
  }
}

