# -----------------------------------------------------------------------------#
#' index vectors grouped by category
#' 
#' Get the index vectors grouped by category given by \code{grp}
#'
#' @param tbl a table to search for
#' @param grp one or more characters of the column names of \code{tbl} or a
#'   vector of the length of the number of columns of the \code{tbl}
#'
#' @return a list of index vectors
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @noRd
# -----------------------------------------------------------------------------#
# created  : 2017-08-20 by Mun-Gwan
# modified : 
# -----------------------------------------------------------------------------#
index_grouped_by_cat <- function(tbl, grp) {
  t_len <- nrow(tbl)
  stopifnot(is.data.frame(tbl))
  if(is.character(grp)) {
    # if one or more characters are given for 'grp', then find the columns
    # having the names and use them for stratification
    if(all(grp %in% names(tbl))) {
      return( by(seq_len(t_len), tbl[grp], c) )
    } else if(length(grp) < 6) {     # small number but no matched column
      stop("'", paste(grp, collapse= "', '"), 
           "' column cannot be found in the given 'tbl'.")
    }
  }
  
  if(length(grp) != t_len) 
    stop("length('grp') should be same as nrow('tbl').")
  
  by(seq_len(t_len), grp, c)
}


# -----------------------------------------------------------------------------#
#' get index in s_b by category
#' 
#' Get a vector of indices grouped by a category variable in \code{@sinfo} or
#' \code{@binder}
#'
#' @param baf an object of \code{\link{BAf-class}}
#' @param by_x by which batch IDs ('sinfo' or 'binder'). The argument must be
#'   given using 'by_s' or 'by_b' because the type of the batch ID is determined
#'   by reading the 4th letter.
#'
#' @return a list of index vectors
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @noRd
# -----------------------------------------------------------------------------#
# created  : 2017-08-19 by Mun-Gwan
# modified : 
# -----------------------------------------------------------------------------#
.index_grouped_by_s_or_by_b <- function(baf, by_x) {
  sb_in <- deparse(substitute(by_x))
  if(sb_in == "") stop("The argument 'by_x' should not be missed. ")
  s_b <- sb_in %>% 
    substr(4, 4) %>% 
    s_b_switch()
  i_len <- switch(s_b, sinfo = nrow(baf), binder = ncol(baf))
  if(is.null(by_x)) return(list(seq_len(i_len)))
  
  index_grouped_by_cat(slot(baf, s_b), by_x)
}



