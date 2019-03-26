# -----------------------------------------------------------------------------#
#' Access data stored in a BAf object
#' 
#' Short commands to access the information stored in a \code{\link{BAf-class}}.
#' \code{sI}, \code{sX}, \code{sB}, \code{sA}, \code{codebook}, and \code{note}
#' invoke the stored information in \code{@sinfo}, \code{@.Data},
#' \code{@binder}, \code{@assy_s} or \code{@assy_b}, \code{@codebook}, and
#' \code{note}, respectively. The details of the contents can be found in
#' \link{BAf-class}. The corresponding assignment functions e.g. \code{sI<-},
#' \code{sX<-}, \code{sB<-}, are also available.
#' 
#' @name data_access   
#'   
#' @param x an object of the \code{\link{BAf-class}}
#' @param value a \code{\link{tbl_df}}, \code{\link{data.frame}} or
#'   \code{\link{matrix}} that fits to the accessing slot
#' @param ... for other functions with same name
#' 
#' @return 
#' \tabular{ll}{
#'    \code{sI}\tab a \code{\link{tbl_df}} of sample information of the given
#'    BAf\cr 
#'    \code{sX}\tab a \code{\link{matrix}} or \code{\link{vector}} of numerical
#'    measurement values obtained from binder arrays\cr 
#'    \code{sB}\tab a \code{\link{data.frame}} of the information about
#'    binders\cr 
#'    \code{sA}\tab a \code{list} of \code{\link{data.frame}} of batch (or 
#'    assay) dependent information for the samples or binders\cr
#'    \code{codebook}\tab a \code{list} of two \code{\link{data.frame}}s named 
#'    as \code{sinfo} and \code{binder} unless \code{wise} is given. Each
#'    contains the codebook for the \code{@sinfo} or \code{@binder} above.\cr
#'    \code{note}\tab any note in any form
#' }
#' 
#' @examples
#' data(sba)
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
# -----------------------------------------------------------------------------#
#' @examples
#' sI(sba)
#' 
#' @aliases sI
#' @rdname data_access
#' @export
# -----------------------------------------------------------------------------#
# created  : 2013-03-28 by Mun-Gwan
# modified : 2017-08-17 by Mun-Gwan : to adapt to BAf-class
# -----------------------------------------------------------------------------#
setGeneric("sI", function(x, ...) standardGeneric("sI"));

# -----------------------------------------------------------------------------#
#' @rdname data_access
#' @export
# -----------------------------------------------------------------------------#
setMethod("sI", signature(x = "BAf"), function(x) x@sinfo)


# -----------------------------------------------------------------------------#
#' @rdname data_access
#' @export
# -----------------------------------------------------------------------------#
setGeneric("sI<-", function(x, ..., value) standardGeneric("sI<-"));

# -----------------------------------------------------------------------------#
#' @rdname data_access
#' @export
# -----------------------------------------------------------------------------#
setMethod("sI<-", signature(x = "BAf"), function(x, value) {
  if(! inherits(value, "tbl_df")) value <- dplyr::as_tibble(value)
  x@sinfo <- value
  
  validObject(x)
  return(x)
})


# -----------------------------------------------------------------------------#
#' @param b_id one or more binder IDs. This allows to find the data of those
#'   binders ignoring repeat indicator, like \code{sX(x)[, match(b_id, bid(x,
#'   exact = F))]}. For an exact match, use \code{`[`} as \code{sX(x)[, b_id]},
#'   instead.
#' @param drop same as \code{\link{[}}
#' 
#' @examples
#' sX(sba)[1:5, 1:5]           # (entire matrix)[1:5, 1:5]
#' sX(sba, "T2")[1:10]         # only the values obtained by 'T2' binder [1:10]
#' dim(sX(sba, "T2", drop= FALSE))
#' 
#' @rdname data_access
#' @export
# -----------------------------------------------------------------------------#
setGeneric("sX", function(x, ...) standardGeneric("sX"));

# -----------------------------------------------------------------------------#
#' @rdname data_access
#' @export
# -----------------------------------------------------------------------------#
setMethod("sX", signature(x = "BAf"), function(x,
                                               b_id,
                                               drop = TRUE) {
  
  if(missing(b_id)) return(x@.Data)
  
  stopifnot(is.character(b_id))
  
  if(any(grepl("\\*[[:digit:]]+$", b_id))) {
    warning("The repeat indicator in the given 'b_id' will be ignored.")
    b_id <- rm_repeat_indicator(b_id)
  }
  
  stopifnot(all(b_id %in% bid(x, exact = FALSE)))
  
  x@.Data[, match(b_id, bid(x, exact = FALSE)), drop= drop]
})


# -----------------------------------------------------------------------------#
#' @rdname data_access
#' @export
# -----------------------------------------------------------------------------#
setGeneric("sX<-", function(x, ..., value) standardGeneric("sX<-"));

# -----------------------------------------------------------------------------#
#' @rdname data_access
#' @export
# -----------------------------------------------------------------------------#
setMethod("sX<-", signature(x = "BAf"), function(x, value) {
  stopifnot(is.numeric(value))
  x@.Data <- value
  
  validObject(x)
  return(x)
})


# -----------------------------------------------------------------------------#
#' @examples
#' sB(sba)
#' 
#' @rdname data_access
#' @export
# -----------------------------------------------------------------------------#
setGeneric("sB", function(x, ...) standardGeneric("sB"));

# -----------------------------------------------------------------------------#
#' @rdname data_access
#' @export
# -----------------------------------------------------------------------------#
setMethod("sB", signature(x = "BAf"), function(x) x@binder)

# -----------------------------------------------------------------------------#
#' @rdname data_access
#' @export
# -----------------------------------------------------------------------------#
setGeneric("sB<-", function(x, ..., value) standardGeneric("sB<-"));

# -----------------------------------------------------------------------------#
#' @rdname data_access
#' @export
# -----------------------------------------------------------------------------#
setMethod("sB<-", signature(x = "BAf"), function(x, value) {
  x@binder <- value
  validObject(x)
  return(x)
})

# -----------------------------------------------------------------------------#
#' @param wise which batch, sample-wise or binder-wise.
#' 
#' @examples
#' head(sA(sba, "s")$fail_flag)
#' sA(sba, "sinfo")$fail_flag[1:5, ]
#' sA(sba, "sample")$fail_flag[1, ] <- FALSE      # 'sample' instead of 'sinfo'
#' sA(sba, "sinfo")$fail_flag[1:5, ]
#' head(sA(sba, "binder")$fail_flag)
#' 
#' @rdname data_access
#' @export
# -----------------------------------------------------------------------------#
setGeneric("sA", function(x, ...) standardGeneric("sA"));

# -----------------------------------------------------------------------------#
#' @rdname data_access
#' @export
# -----------------------------------------------------------------------------#
setMethod(
  "sA", signature(x = "BAf"),
  function(x, wise= c("sinfo", "binder")) {
    s_b_switch(wise) %>%     # allow "sample", too
      switch(., 
             sinfo = x@assy_s, 
             binder = x@assy_b)
  }
)

# -----------------------------------------------------------------------------#
#' @rdname data_access
#' @export
# -----------------------------------------------------------------------------#
setGeneric("sA<-", function(x, ..., value) standardGeneric("sA<-"));

# -----------------------------------------------------------------------------#
#' @rdname data_access
#' @export
# -----------------------------------------------------------------------------#
setMethod(
  "sA<-", signature(x = "BAf"), 
  function(x, wise= c("sinfo", "binder"), value) {
    wise <- s_b_switch(wise)    # allow "sample", too

    if(wise == "sinfo") x@assy_s <- value
    else x@assy_b <- value
    
    validObject(x)
    return(x)
  }
)


# -----------------------------------------------------------------------------#
#' @examples
#' head(codebook(sba))              # both codebooks for samples & binders
#' codebook(sba, "bin")         # short name
#' codebook(sba, "sample")      # 'sample' instead of 'sinfo'
#' 
#' @rdname data_access
#' @export
# -----------------------------------------------------------------------------#
setGeneric("codebook", function(x, ...) standardGeneric("codebook"));

# -----------------------------------------------------------------------------#
#' @rdname data_access
#' @export
# -----------------------------------------------------------------------------#
setMethod(
  "codebook", signature(x = "BAf"),
  function(x, wise= c("sinfo", "binder")) {
    if(all(wise == c("sinfo", "binder"))) return(x@codebook)

    s_b_switch(wise) %>%     # allow "sample", too
      x@codebook[[ . ]]
  }
)

# -----------------------------------------------------------------------------#
#' @rdname data_access
#' @export
# -----------------------------------------------------------------------------#
setGeneric("codebook<-", function(x, ..., value) standardGeneric("codebook<-"));

# -----------------------------------------------------------------------------#
#' @rdname data_access
#' @export
# -----------------------------------------------------------------------------#
setMethod(
  "codebook<-", signature(x = "BAf"), 
  function(x, wise= c("sinfo", "binder"), value) {
    if(all(wise == c("sinfo", "binder"))) {
      x@codebook <- value
    } else {
      wise <- s_b_switch(wise)    # allow "sample", too
      x@codebook[[ wise ]] <- value
    }
    validObject(x)
    return(x)
  }
)

# -----------------------------------------------------------------------------#
#' @rdname data_access
#' @export
# -----------------------------------------------------------------------------#
setGeneric("note", function(x) standardGeneric("note"));

# -----------------------------------------------------------------------------#
#' @rdname data_access
#' @export
# -----------------------------------------------------------------------------#
setMethod("note", signature(x = "BAf"), function(x) x@note)

# -----------------------------------------------------------------------------#
#' @rdname data_access
#' @export
# -----------------------------------------------------------------------------#
setGeneric("note<-", function(x, value) standardGeneric("note<-"));

# -----------------------------------------------------------------------------#
#' @rdname data_access
#' @export
# -----------------------------------------------------------------------------#
setMethod("note<-", signature(x = "BAf"), function(x, value) {
  x@note <- value
  validObject(x)
  return(x)
})


# -----------------------------------------------------------------------------#
#' as.matrix.BAf
#' 
#' This is defined because \code{as.matrix} returns a \code{BAf} object when it
#' was given as a input of the function.
#' 
#' @param x an object of the \code{\link{BAf-class}}
#' @param ... not used
#' 
#' @return a \code{\link{matrix}} that contains all measurement values
#' @examples
#' data(sba)
#' as.matrix(sba)
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @export
# -----------------------------------------------------------------------------#
# created  : 2017-08-18 by Mun-Gwan
# modified : 
# -----------------------------------------------------------------------------#

as.matrix.BAf <- function(x, ...) x@.Data
