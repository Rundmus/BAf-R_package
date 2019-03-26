# -----------------------------------------------------------------------------#
#' divide by a column in @sinfo
#' 
#' It produces a list of BAf objects divided by the values in the given column
#' of \code{@sinfo}
#' 
#' @param baf an object of the \code{\link{BAf-class}}
#' @param by the name of the column by which sample will be divided
#' @param abbr_name if the name of the list element should be shortened
#' 
#' @return A \code{list} of \code{\link{BAf-class}} objects
#' @examples
#' data(sba)
#' # the data "sba" is divided based on the values in "plate" column.
#' sba2 <- divideBy(sba, "plate")
#' str(sba2, 2)		# show the structure of "B"
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @export
# -----------------------------------------------------------------------------#
# created  : 2012-04-13 by Mun-Gwan
# modified : 2012-04-25 by Mun-Gwan : add '
# -----------------------------------------------------------------------------#

divideBy <- function(baf, by, abbr_name = FALSE) {
  stopifnot(inherits(baf, "BAf"))
  stopifnot(by %in% colnames(baf@sinfo))
  
  by_col <- baf@sinfo[, by]
  out <- by_col %>% 
    unlist() %>% as.character() %>% unique() %>%
    na.omit() %>%
    sapply(function(eL) {
      baf[!is.na(by_col) & by_col == eL, ]
    }, simplify= F) %>% { 
      c(list("all"= baf), .) 
    }
  
  ## abbreviate names
  if(abbr_name) {
    names(out) <- names(out) %>% abbreviate(1)
  }
  return(out)
}

