# -----------------------------------------------------------------------------#
#' summary of the BAf object
#' 
#' This shows the summary of \code{@sinfo}, \code{@binder}, and \code{@.Data}
#' 
#' @param object an object of the \code{\link{BAf-class}}
#' @param ... not used
#' 
#' @return a list of summarized elements
#' 
#' @examples
#' data(sba)
#' summary(sba)
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @rdname summary
#' @export
# -----------------------------------------------------------------------------#
# created  : 2011-07-12 by Mun-Gwan
# modified : 2013-06-27 by Mun-Gwan : adapt to "SBAe" class
#            2016-01-22 by Mun-Gwan : change from S3 to S4 for 'print'
#            2017-08-17 by Mun-Gwan : to adapt to BAf-class
# -----------------------------------------------------------------------------#

setGeneric("summary")

# -----------------------------------------------------------------------------#
#' @rdname summary
#' @export
# -----------------------------------------------------------------------------#
setMethod("summary", signature(object = "BAf"), function(object, ...) {
  validObject(object)
  
  list('sinfo' = summary(object@sinfo),          # 'table' chr class
       'binder'= summary(object@binder),
       'dim'   = dim(object@.Data),
       'part'  = summary(object@.Data[, 1:min(ncol(object@.Data), 5)])
  ) %>% 
    `class<-`("summary.BAf")
})
setClass("summary.BAf")


#' @export
print.summary.BAf <- function(x, ...) {
  cat0 <- function(...) cat(..., sep= "")
  cat0("Dimension of this BAf data is [", paste(x$dim, collapse= ", "), "].\n")
  
  show_max_2_lines <- function(x, var) {
    cat0("@", var, "\n")           # show what
    xv <- x[[var]]
    icol <- xv %>% 
      apply(2, . %>% na.omit %>% nchar %>% max) %>%      # width of each column
      cumsum %>% 
      #  1st line (if shorter than width)
      { . - .[max(which(. < getOption("width")))] } %>%
      { which(. < getOption("width")) }          # max 2 lines
    print(xv[, icol])
    if((o_n <- ncol(xv) - max(icol)) > 0) {
      oV <- colnames(xv)[-icol]
      if(length(oV) > 10) c(oV[1:10], "...")
      cat0("# Note: There are ", o_n, " more columns in @", var, 
           paste(oV, collapse= ","), "\n")
    }
    cat("\n")
  }
  
  ## show @sinfo
  show_max_2_lines(x, "sinfo")
  ## show @binder
  show_max_2_lines(x, "binder")
  ## show @.Data
  cat("@.Data\n")
  if(x$dim[2L] > 5) 
    cat0("The summary of the first 5 columns are shown below.\n")
  # show the first 5 columns only
  print(x$part)
}
