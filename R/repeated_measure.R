# -----------------------------------------------------------------------------#
#' Transform to landscape form 
#' 
#' When a \code{\link{BAf-class}} object contains the data of multiple repeats
#' in one assay, this function transform the data into landscape form, in which
#' the repeats are extended into columns.
#' 
#' @param baf an object of \code{\link{BAf-class}}
#' 
#' @return \code{\link{BAf-class}} object in landscape form
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @importFrom Useful2me combine_multi_elements
#' @export
# -----------------------------------------------------------------------------#
# created  : 2016-10-18 by Mun-Gwan
# modified : 
#   2017-08-22 by Mun-Gwan : to adapt to BAf-class
#   2018-04-20 by Mun-Gwan : 
#     1) use 'sid()' instead of 'sid_i' argument
#     2) no 'set_i'
# -----------------------------------------------------------------------------#

transform_to_landscape_for_repeats <-
  function(baf) {
    
    ##  Check given arguments
    stopifnot(inherits(baf, "BAf"))

    ##  * FIND OUT REPEATED SAMPLES *  ##
    
    #  n_rep = number of repeats, (names = ID of repeated samples only)
    n_rep <- table(sid(baf)) %>% 
      .[. > 1]
    if(length(n_rep) == 0) stop("No repeated sample ID")
    if(!all(n_rep == n_rep[1])) 
      stop("The number of repeats should be constant for all samples.")
    
    isRepeated <- sid(baf) %in% names(n_rep)
    n_rep <- n_rep[1]
    
    baf_rep <- baf[isRepeated, ] %>%            # only repeated samples
      .[order(sid(.)), ]
    
    ##  * Divide data by repeat number *  ##
    baf_rep <- lapply(1:n_rep, function(ii) {
      baf_rep[(seq_len(nrow(baf_rep) / n_rep) - 1) * n_rep + ii, ]
    })
    
    ##  Confirm all sample IDs are identical  ##
    for(ii in 2:n_rep) {
      stopifnot(identical(sid(baf_rep[[1]]), sid(baf_rep[[ii]])))
    }
    
    ###   Preparation for 'cbind'   ###
    
    ##  Combine sample information for same samples
    tmp <- colnames(sI(baf_rep[[1]]))
    sinfo_n <- lapply(tmp, function(eC) {
      # eCb <- do.call("cbind", lapply(baf_rep, . %>% {.@sinfo[eC]}))
      do.call("cbind", lapply(baf_rep, . %>% {.@sinfo[eC]})) %>% 
        apply(1, combine_multi_elements) 
    }) %>% 
      setNames(tmp) %>% 
      as_tibble() %>% 
      mutate_at(tmp[sapply(sI(baf_rep[[1]]), is.factor)], as.factor)

    for(ii in seq_along(baf_rep)) { 
      sI(baf_rep[[ii]]) <- sinfo_n %>% 
        mutate(id= sid(baf_rep[[ii]], exact= T))
      batch(baf_rep[[ii]], "binder") <- batch(baf_rep[[ii]], "binder") %>% 
        paste0(".", ii)
    }
    
    #  ***  'cbind' to make it landscape form  ***  #
    out <- do.call("cbind.BAf", baf_rep)
    
    ##  Add 'rep_binder' in @binder to be used in other functions 
    ##  (e.g. hist_correlation)
    tmp <- make.unique(c(colnames(sB(out)), "rep_binder"))      # to avoid same name
    sB(out)[tmp[length(tmp)]] <- bid(baf_rep[[1]]) %>%
      make.unique_key_ids() %>%           # refresh
      rep(., length(baf_rep)) %>%
      factor()
    
    return(out)
  }


