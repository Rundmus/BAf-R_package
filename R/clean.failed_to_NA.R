# -----------------------------------------------------------------------------#
#' Failed sample/binder to NA
#' 
#' replace the data points of the failed samples/binders with NA
#' 
#' 
#' @param baf an object of the \code{\link{BAf-class}}
#' @param wise sample-wise or binder-wise.
#' @param fail_flag the name of the data frame that contains fail flags in
#'   \code{@assy_s} or \code{@assy_b}
#' @param failMark what is the label in the flag for failed sample
#' @param show_count Whether the number of "failed" is shown
#' @param by if show_count is TRUE, then "by" indicates by which column the 
#'   counts are stratified, e.g. "plate", "assay". This can be given as a vector
#'   like the result of \code{batch(baf, "sinfo")}, which is the default.
#'        
#' @return an object of the \code{\link{BAf-class}}
#' @examples 
#' data(sba)
#' failed_to_NA(sba, wise= "sinfo", by= "plate")
#' failed_to_NA(sba)         # default 'by' is batch
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @export
# -----------------------------------------------------------------------------#
# created  : 2011-04-17 by Mun-Gwan
# modified : 
#   2017-08-20 by Mun-Gwan : to adapt to BAf-class
# -----------------------------------------------------------------------------#
failed_to_NA <- function(baf, 
                         wise= c("sinfo", "binder"),
                         fail_flag = "fail_flag",
                         failMark = "failed",
                         show_count = TRUE,
                         by) {
  wise <- s_b_switch(wise)    # allow "sample", too
  fF <- switch(
    wise,
    sinfo = {
      stopifnot(fail_flag %in% names(baf@assy_s))
      baf@assy_s[[fail_flag]]
    }, 
    binder = {
      stopifnot(fail_flag %in% names(baf@assy_b))
      baf@assy_b[[fail_flag]]
    }
  )
  
  if(show_count) {
		cat("The number of \"failed\"\n")
		print(apply(fF, 2, function(eC) sum(eC == failMark)))

		if(missing(by)) by <- batch(baf, wise)
		print(apply(fF, 2, function(eC) {
		  index_grouped_by_cat(slot(baf, wise), by) %>% 
		    sapply(., . %>% { sum(eC[.] == failMark) })
		}))
	}
	
	##  replacing all data of failed samples / binders
  switch(
    wise, 
    sinfo = {
      for(j in colnames(fF)) {
        sX(baf)[fF[, j] == failMark, batch(baf, "binder") == j] <- NA
      }
    }, 
    binder = {
      for(j in colnames(fF)) {
        sX(baf)[batch(baf, "sinfo") == j, fF[, j] == failMark] <- NA
      }
    }
  )
	return(invisible(baf))
}




