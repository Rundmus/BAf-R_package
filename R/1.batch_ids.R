# -----------------------------------------------------------------------------#
#' Access batch IDs 
#' 
#' \code{batch} gets and \code{batch<-} assigns a vector of batch IDs for
#' binders/samples. The ID is used to distinguish separate assay runs.
#' 
#' @param x an object of the \code{\link{BAf-class}}
#' @param value a \code{\link{vector}} of new batch IDs. It should be either a
#'   single value or in exactly same dimension as the previous values.
#' @param wise which batch, sample-wise or binder-wise.
#' @param ... for other functions with same name
#' 
#' @return 
#' \code{batch}: a \code{vector} of batch IDs for binders/samples
# -----------------------------------------------------------------------------#
#' @examples
#' data(sba)
#' batch(sba, "binder")               # get batch IDs
#' batch(sba, "binder") <- rep(c("AY10", "AY11"), each= ncol(sba) / 2)
#' batch(sba, "b")               # short name is also allowed.
#' 
#' batch(sba, "sinfo")
#' batch(sba, "sinfo") <- rep(LETTERS[1:4], each= 96)
#' batch(sba, "sinfo")[1:96] <- "A"
#' batch(sba, "sample")     # 'sample' instead of 'sinfo' is also allowed.
#' 
#' @aliases batch
#' @rdname batch
#' @export
# -----------------------------------------------------------------------------#
# created  : 2017-08-17 by Mun-Gwan
# modified : 
# -----------------------------------------------------------------------------#
setGeneric("batch", function(x, ...) standardGeneric("batch"));

# -----------------------------------------------------------------------------#
#' @rdname batch
#' @export
# -----------------------------------------------------------------------------#
setMethod(
  f= "batch", signature(x = "BAf"), 
  function(x, wise= c("sinfo", "binder")) {
    wise <- s_b_switch(wise)

    slot(x, wise)[[ x@batch_c[[wise]] ]] %>% as.character()
  }
)


# -----------------------------------------------------------------------------#
#' @rdname batch
#' @export
# -----------------------------------------------------------------------------#
# created  : 2017-08-17 by Mun-Gwan
# modified : 
# -----------------------------------------------------------------------------#
setGeneric("batch<-", function(x, ..., value) standardGeneric("batch<-"));

# -----------------------------------------------------------------------------#
#' @rdname batch
#' @export
# -----------------------------------------------------------------------------#
setMethod(
  f= "batch<-", signature(x = "BAf"), 
  function(x, wise= c("sinfo", "binder"), value) {
    wise <- s_b_switch(wise)
    
    switch(wise,
           sinfo = stopifnot(length(value) == 1 | nrow(x) == length(value)),
           binder = stopifnot(length(value) == 1 | ncol(x) == length(value)))

    if(!is.factor(value)) value <- as.factor(value)
    s_as <- switch(wise, sinfo = "assy_b", binder = "assy_s")
    prev <- batch(x, wise)      # character, not factor
    
    ##  Replace batch ID in @assy_?
    if(length(value) == 1L) {
      ##  if 'value' is single
      
      #  check a single batch ID before replacement
      if(length(unique(prev)) != 1L)
        stop("A single 'value' is allowed only when the 'x' had",
             " one batch ID previously.")
      
      for(j in seq(slot(x, s_as))) {      # per element of @assy_?
        colnames(slot(x, s_as)[[j]]) <- as.character(value)
      }
    } else {
      ##  if 'value' has a vector of multiple values
      
      #  check if one-to-one match is possible
      if(length(unique(prev)) != length(levels(value))) 
        stop("Only one-to-one replacement is allowed.")
      
      bid_link <- tibble(new = levels(value)) %>% 
        mutate(old = character(n()))
      ##  For multiple batch IDs
      i_ob <- by(seq(prev), prev, c)      # divide by old batch ID
      
      for(j in seq(i_ob)) {      # per old assay ID
        new_bid <- value[ i_ob[[j]] ] %>% as.character
        if(length(new_bid %>% unique) != 1) # only one-to-one match
          stop("All samples/binders of a batch should be replaced ",
               "with a single batch ID.")
        bid_link$old[ bid_link$new == new_bid[1] ] <- names(i_ob)[j]
      }
      for(j in seq(slot(x, s_as))) {      # per element of @assy_?
        colnames(slot(x, s_as)[[j]]) <- 
          bid_link$new[match(colnames(slot(x, s_as)[[j]]), bid_link$old)]
      }
    }
    
    #  batch IDs in @binder
    slot(x, wise)[[ x@batch_c[[wise]] ]] <- value
    
    validObject(x)
    return(x)
  }
)


# -----------------------------------------------------------------------------#
#' \code{batch_colname} returns the column name having sample/binder batch IDs
#' in \code{@sinfo} / \code{@batch}
#' 
#' @return 
#'  \code{batch_colname}: a character of the name of the column that contains
#'  sample/binder batch IDs\cr
#'  
#' @rdname batch
#' @export
# -----------------------------------------------------------------------------#
# created  : 2017-08-20 by Mun-Gwan
# modified : 
# -----------------------------------------------------------------------------#
setGeneric("batch_colname", function(x, ...) standardGeneric("batch_colname"))

# -----------------------------------------------------------------------------#
#' @rdname batch
#' @export
# -----------------------------------------------------------------------------#
setMethod(
  "batch_colname", signature(x = "BAf"), 
  function(x, wise= c("sinfo", "binder")) {
    s_b_switch(wise) %>%     # allow "sample", too
      x@batch_c[[ . ]]
  }
)



# -----------------------------------------------------------------------------#
#' \code{n_batch} returns the number of batches
#' 
#' @return 
#'  \code{n_batch}: a list of the numbers of sample/binder batches
#'
#' @examples
#' n_batch(sba)
#'  
#' @rdname batch
#' @export
# -----------------------------------------------------------------------------#
# created  : 2017-09-19 by Mun-Gwan
# modified : 
# -----------------------------------------------------------------------------#
setGeneric("n_batch", function(x, ...) standardGeneric("n_batch"))

# -----------------------------------------------------------------------------#
#' @rdname batch
#' @export
# -----------------------------------------------------------------------------#
setMethod(
  "n_batch", signature(x = "BAf"), 
  function(x) {
    sapply(c("sinfo", "binder"), function(wise) {
      batch(x, wise) %>% 
        unique() %>% 
        length()
    }, simplify= FALSE)
  }
)

