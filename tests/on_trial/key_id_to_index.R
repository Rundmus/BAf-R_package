# -----------------------------------------------------------------------------#
#' Key ID to index
#' 
#' Change key ID to index (\code{numeric}). This allows to accept different form
#' of key IDs such as logical vector.
#' 
#' @param baf an object of the \code{\link{BAf-class}}
#' @param k_id one or more key IDs. It accepts a numerical index or a logical
#'   vector. If this is missing, return TRUE.
#' @param wise sample IDs or binder IDs.
#' @param exact This is valud only when \var{k_id} was given as characters. If
#'   it is FALSE (as default), the indices of all \var{k_id} ignoring the repeat
#'   indicator will be returned. Otherwise (\code{TRUE}), the indices of exactly
#'   matched key IDs will be returned.
#' 
#' @examples
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
# -----------------------------------------------------------------------------#
# created  : 2017-08-31 by Mun-Gwan
# modified : 
# -----------------------------------------------------------------------------#

.key_id_to_index <- function(baf,
                             k_id,
                             wise = c("sinfo", "binder"),
                             exact = FALSE) {
  
  stopifnot(inherits(baf, "BAf"),
            !missing(k_id))
  
  wise <- s_b_switch(wise)    # allow "sample", too
  ids <- switch(wise,
                sinfo  = sid(baf, exact= exact), 
                binder = bid(baf, exact= exact))
  
  if(is.character(k_id)) {
    # warn if exact == F but any ID has the repeat indicator 
    if(!exact && any(grepl("\\*[[:digit:]]+$", k_id))) {
      warning("The repeat indicator in the given key IDs will be ignored ",
              "when the 'exact' is FALSE.")
      k_id <- rm_repeat_indicator(k_id)
    }
    
    if(! all(k_id %in% ids))
      stop("All given key IDs should be in @", wise)
    
    match(k_id, ids)
    
  } else if(is.logical(k_id)) {
    if(length(k_id) != length(ids))
      stop("The given key IDs should have same length as IDs in @", wise)
    
    which(k_id)
    
  } else {
    seq_len(length(ids))[k_id]
  } 
}
  