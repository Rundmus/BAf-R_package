# -----------------------------------------------------------------------------#
#' 'sinfo' and 'binder' switch
#' 
#' By pattern matching, a shorter character is extended to a full name
#'
#' @param s_b a character that can be matched to 'sinfo' or 'binder'.
#'   Exceptionally, 'sample' is matched to 'sinfo' in order to allow more
#'   readable code.
#'
#' @return a character with full name. i.e. 'sinfo' or 'binder'
#' @examples
#' s_b_switch("s")
#' @noRd
# -----------------------------------------------------------------------------#
# created  : 2017-08-19 by Mun-Gwan
# modified : 
# -----------------------------------------------------------------------------#
s_b_switch <- function(s_b= c("sinfo", "binder")) {
  stopifnot(is.character(s_b))
  if(s_b[1] == "sample") return("sinfo")
  match.arg(s_b)
}


# -----------------------------------------------------------------------------#
#' Handle repeat indicator of the key IDs
#' 
#' The key IDs (\code{'id'} in \code{@sinfo} or \code{@binder}) should be unique
#' throughout the samples and binders, in order to distinguish among them. So,
#' any \code{id}s repeatedly used are amended by adding "*" and the digit
#' indicating how many times the same ids have been employed ahead as a suffix
#' (e.g. HPA0001, HPA0001*1)
#' \code{rm_repeat_indicator} - Remove the additional suffix (repeat indicator)
#' 
#' @name repeat_indicator
#' @param ids a \code{vector} of the key IDs (\code{'id'} in \code{@sinfo} or
#'   \code{@binder}). It is by default converted to \code{character}.
#' @return \code{rm_repeat_indicator} A \code{vector} of \code{character} of
#'   the 'id' after removal of the indicator
#'   
#' @examples
#' rm_repeat_indicator(c("HPA001", "HPA001*1"))
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @rdname repeat_indicator
#' @export
# -----------------------------------------------------------------------------#
# created  : 2014-11-04 by Mun-Gwan
# modified :
# -----------------------------------------------------------------------------#

rm_repeat_indicator <- function(ids) {
  ids %>% 
    as.character() %>% { 
      sub("\\*[[:digit:]]+$", "", .) 
    }
}

# -----------------------------------------------------------------------------#
#' @description 
#' \code{make.unique_key_ids} generate unique 'id's by adding the suffices. When
#' necessary, the indicators previously attached are removed first.
#' 
#' @return \code{make.unique_key_ids} A \code{vector} of \code{character} of 
#'   cleaned 'id's with the suffix
#' 
#' @examples
#' make.unique_key_ids(c("HPA001", "HPA003", "HPA001", "HPA001", "HPA003"))
#' 
#' @rdname repeat_indicator
#' @export
# -----------------------------------------------------------------------------#
# required fn
#   - rm_repeat_indicator
# created  : 2014-11-04 by Mun-Gwan
# modified :
# -----------------------------------------------------------------------------#

make.unique_key_ids <- function(ids) {
  ids %>%
    rm_repeat_indicator() %>%
    make.unique(sep="*")
}

# -----------------------------------------------------------------------------#
#' Reset key Ids
#' 
#' Reset both key Ids (\code{'id'} in \code{@sinfo} and \code{@binder}). The
#' default key IDs are those from the 'id' in @sinfo or in @binder rather than
#' those in @.Data.
#' 
#' @param baf an object of the \code{\link{BAf-class}}
#' @return an object of the class \link{BAf-class} after modifying the key IDs
#'   and checking the validity
#'   
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @noRd
# -----------------------------------------------------------------------------#
# required fn
#   - make.unique_key_ids
# created  : 2017-07-31 by Mun-Gwan
# modified : 2017-08-17 by Mun-Gwan : to adapt to BAf-class
# -----------------------------------------------------------------------------#

reset_key_ids <- function(baf, sid= baf@sinfo$id, bid= baf@binder$id) {
  stopifnot(inherits(baf, "BAf"))
  
  rownames(baf@.Data) <- baf@sinfo$id  <- make.unique_key_ids(sid)
  colnames(baf@.Data) <- baf@binder$id <- make.unique_key_ids(bid)
  
  for(ii in seq_along(baf@assy_s)) {
    rownames(baf@assy_s[[ii]]) <- rownames(baf@.Data) 
  }
  for(ii in seq_along(baf@assy_b)) {
    rownames(baf@assy_b[[ii]]) <- colnames(baf@.Data) 
  }
  validObject(baf)
  baf
}



# -----------------------------------------------------------------------------#
#' Access key ids 
#' 
#' @name key_id_access   
#'   
#' @param x an object of the \code{\link{BAf-class}}
#' @param value a \code{\link{vector}} of key ids
#' @param exact whether the exact ID including repeat indicator, i.e. \code{*1},
#'   should be returned or not
#' 
#' @return 
#' \tabular{ll}{
#'    \code{sid}\tab a \code{vector} of sample IDs\cr
#'    \code{bid}\tab a \code{vector} of binder IDs\cr
#' }
#' 
#' @examples
#' data(sba)
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @rdname key_id_access
#' 
# -----------------------------------------------------------------------------#


# -----------------------------------------------------------------------------#
#' @description 
#' \code{sid} returns a vector of sample IDs.\cr
#' \code{sid<-} replaces it with new \code{value}.
#' 
#' @examples
#' head(sid(sba))
#' 
#' @aliases sid
#' @rdname key_id_access
#' @export
# -----------------------------------------------------------------------------#
# created  : 2017-07-31 by Mun-Gwan
# modified :
# -----------------------------------------------------------------------------#
setGeneric("sid", function(x, exact = FALSE) standardGeneric("sid"));

# -----------------------------------------------------------------------------#
#' @rdname key_id_access
#' @export
# -----------------------------------------------------------------------------#
# required fn
#   - rm_repeat_indicator

setMethod("sid", signature(x = "BAf"), function(x, exact = FALSE) {
  rownames(x@.Data) %>% { 
    if(exact) . else rm_repeat_indicator( . )
  }
})

# -----------------------------------------------------------------------------#
#' @examples
#' sid(sba)[1] <- "sA"
#' 
#' @rdname key_id_access
#' @export
# -----------------------------------------------------------------------------#
# created  : 2012-10-24 by Mun-Gwan
# modified : 2013-06-27 by Mun-Gwan : adapt to "SBAe" class
# -----------------------------------------------------------------------------#
setGeneric("sid<-", function(x, value) standardGeneric("sid<-"));

# -----------------------------------------------------------------------------#
#' @rdname key_id_access
#' @export
# -----------------------------------------------------------------------------#
# required fn
#   - reset_key_ids

setMethod("sid<-", signature(x = "BAf"), function(x, value) {
  stopifnot(nrow(x) == length(value))
  if(!is.character(value)) value <- as.character(value)
  
  x %>% reset_key_ids(sid= value)
})


# -----------------------------------------------------------------------------#
#' @description 
#' \code{bid} returns a vector of binder IDs.\cr
#' \code{bid<-} replace it with new \code{value}.
#' 
#' @examples 
#' bid(sba)[1:10]
#' 
#' @rdname key_id_access
#' @export
# -----------------------------------------------------------------------------#
# created  : 2011-07-18 by Mun-Gwan
# modified : 2013-06-27 by Mun-Gwan : adapt to "SBAe" class
# -----------------------------------------------------------------------------#
setGeneric("bid", function(x, exact = FALSE) standardGeneric("bid"));

# -----------------------------------------------------------------------------#
#' @rdname key_id_access
#' @export
# -----------------------------------------------------------------------------#
# required fn
#   - rm_repeat_indicator

setMethod("bid", signature(x = "BAf"), function(x, exact = FALSE) {
  colnames(x@.Data) %>% { 
    if(exact) . else rm_repeat_indicator( . ) 
  }
})


# -----------------------------------------------------------------------------#
#' @rdname key_id_access
#' @export
# -----------------------------------------------------------------------------#
# created  : 2012-10-24 by Mun-Gwan
# modified : 2013-06-27 by Mun-Gwan : adapt to "SBAe" class
# -----------------------------------------------------------------------------#
setGeneric("bid<-", function(x, value) standardGeneric("bid<-"));

# -----------------------------------------------------------------------------#
#' @rdname key_id_access
#' @export
# -----------------------------------------------------------------------------#
# required fn
#   - reset_key_ids

setMethod("bid<-", signature(x = "BAf"), function(x, value) {
  stopifnot(ncol(x) == length(value))
  if(!is.character(value)) value <- as.character(value)
  
  x %>% reset_key_ids(bid= value)
})

