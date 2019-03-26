assy_new_rm_duplicates <- function(xs, s_as) {
  xs1_as <- slot(xs[[1]], s_as)
  if(length(xs1_as) == 0) {
    return(list())
  } else {
    #  batches of all BAf objects stored in column names of @assy_? 
    batid <- xs %>% 
      sapply(. %>% 
               {slot(., s_as)[[1]]} %>% 
               names()
      ) %>% unlist()
    iDupl_bat <- duplicated(batid)
    
    lapply(seq_along(xs1_as), function(j) { # per property
      #  bind the jth property of every BAf in column wise
      do.call("cbind", 
              lapply(xs, . %>% 
                       {slot(., s_as)[[j]]}
              )
      ) %>%
        `names<-`(batid) %>% {
          x <- .
          if(any(iDupl_bat)) {		# if any batch ID in common
            ## check identity of the contents for the same batch ID
            for(iD in which(iDupl_bat)) {
              if( any(x[, match(batid[iD], batid)] != x[, iD]) )
                stop("All data in @", s_as, " for same batch must be same.")
            }
            x[, !iDupl_bat, drop= F]
          } else x
        }
    }) %>% 
      `names<-`(names(xs1_as))
  }
}





# -----------------------------------------------------------------------------#
#' Combine BAf objects by Rows or Columns
#' 
#' Both functions works in almost same way as generic forms of them.\cr 
#' \itemize{ 
#' \item \code{cbind.BAf} binds BAf objects in column-wise. It checks if it is 
#' plausible to bind two or more BAf objects first, such as identical sample
#' information stored in \code{@sinfo} and sample batch-wide data in in
#' \code{@assy_s}.\cr
#' \item \code{rbind.BAf} binds in row-wise. To confirm it is possible to
#' combine multiple data, this checks the consistency in variable names in 
#' \code{@sinfo} and in \code{@assy_s}, the contents in \code{@assy_b} and in
#' \code{@binder}\cr
#' }
#' 
#' @name cbind_rbind
#' 
#' @param ... one or more objects of the \code{\link{BAf-class}}
#' @param deparse.level not used. a variable included in the generic function
#' 
#' @return an object of the \code{\link{BAf-class}} after binding
#' @examples
#' data(sba)
#' summary(sba@binder)
#' 
#' sba2 <- cbind(sba[, 1:3], sba[, 10:20])
#' summary(sba2)
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @export
#' @rdname cbind_rbind
# -----------------------------------------------------------------------------#
# created  : 2011-07-12 by Mun-Gwan
# modified : 
#   2013-06-27 by Mun-Gwan : adapt to "SBAe" class
#   2016-01-22 by Mun-Gwan : change from S3 to S4 for 'print'
#   2017-08-17 by Mun-Gwan : to adapt to BAf-class
# -----------------------------------------------------------------------------#

cbind.BAf <- function(..., deparse.level = 1) {
  xs <- list(...)
  if(! length(xs)) return(NULL)
  stop_if_any_is_not_BAf(xs)
  sapply(xs, validObject)		# check validity of the objects
  
  if(length(xs) == 1) return(xs[[1]])
  
  # vN = the variable names in dots
  vN <- as.character(match.call())[-1]
  
  xs1 <- xs[[1]]
  
  #>>> Check consistency >>>--------------
  for(j in 2:length(xs)) {
    x_j <- xs[[j]]
    # check @sinfo of all entries are same
    if(! isTRUE( all.equal(xs1@sinfo, x_j@sinfo, check.attributes= FALSE) ))
      stop(vN[j], "@sinfo is not identical to ", vN[1], "@sinfo!")
    # check @batch_c$sinfo / $binder
    for(k in c("sinfo", "binder")) {
      if(xs1@batch_c[[k]] != x_j@batch_c[[k]])
        stop("All @batch_c$", k, " must be identical.")
    }
    # check consistency of colnames of @binder
    if(! identical(colnames(xs1@binder), colnames(x_j@binder))) 
      stop(vN[j], "@binder contains different variables from ", vN[1])
    # check consistency of names in @assy_s
    if(! identical(names(xs1@assy_s), names(x_j@assy_s)))
      stop("names(", vN[j], "@assy_s) is not identical to ",
           "names(", vN[1], "@assy_s) !")
    if(! identical(colnames(xs1@assy_b), colnames(x_j@assy_b)))
      stop("The colnames of the elements in @assy_b of all BAfs should be",
           "identical!")
  }
  #<<< check <<<--------------------------
  
  ##  new @binder
  binder_new <- do.call("rbind", lapply(xs, slot, "binder"))
  ##  new @assy_b
  assy_b_new <- lapply(seq_along(xs1@assy_b), function(ii) {
    do.call("rbind", lapply(xs, . %>% {.@assy_b[[ii]]}))
  }) %>% 
    `names<-`(names(xs1@assy_b))
  
  ##  new @assy_s
  assy_s_new <- assy_new_rm_duplicates(xs, "assy_s")

  ##  new @.Data
  new.Data= do.call("cbind", lapply(xs, function(ea) ea@.Data))
  
  BAf(new.Data, 
      sinfo = xs1@sinfo,		# note NO added row or column
      binder = binder_new,
      sinfo_batch_i  = xs1@batch_c$sinfo, 
      binder_batch_i = xs1@batch_c$binder,
      assay_sinfo  = assy_s_new,
      assay_binder = assy_b_new,
      codebook_sinfo  = codebook_bind(xs, "sinfo"), 
      codebook_binder = codebook_bind(xs, "binder"), 
      note = do.call("c", lapply(xs, slot, "note"))
  )
}


# -----------------------------------------------------------------------------#
#' @examples
#' sba3 <- rbind(sba[1:5, ], sba[10:20, ])
#' summary(sba3)
#' 
#' @export
#' @importFrom Useful4me combine_multi_elements
#' @rdname cbind_rbind
# -----------------------------------------------------------------------------#
# modified : 
#   2015-05-18 by Mun-Gwan : allow different assay IDs when same BAf
#   2016-10-20 by Mun-Gwan : allow some difference in @binder
#   2017-08-17 by Mun-Gwan : to adapt to BAf-class
#   2018-04-11 by Mun-Gwan : 
#     1) Remove 'varying' 
#     2) The contents of any varying columns in 'binder' are combined.
# -----------------------------------------------------------------------------#
rbind.BAf <- function(..., deparse.level = 1) {
  xs <- list(...)
  if(! length(xs)) return(NULL)
  stop_if_any_is_not_BAf(xs)
  sapply(xs, validObject)		# check validity of the objects
  
  if(length(xs) == 1) return(xs[[1]])
  
  xs1 <- xs[[1]]
  
  #  has any columns have varying contents
  has_varying_col_binder <- FALSE


  ##  Check consistency of individual elements -------------------------------

  for(j in 2:length(xs)) {
    x_j <- xs[[j]]
    
    if(! identical(names(xs1@sinfo), names(x_j@sinfo))) 
      stop('[', j, "]@sinfo contains different variables from [1]@sinfo")
    
    if(! identical(names(xs1@binder), names(x_j@binder))) 
      stop('[', j, "]@binder contains different variables from [1]@binder")
    
    if(! identical(bid(xs1, exact= T), bid(x_j, exact= T)))
      stop("The binder IDs of all BAfs should be identical!")
    
    #  Whether @binder of all entries are same
    if(! has_varying_col_binder &&
       ! isTRUE( all.equal(xs1@binder, x_j@binder, check.attributes= FALSE) )) {
      message(
        'Some columns in @binder contain varying contents.', "\n",
        'The values in those columns are combined.'
      )
      has_varying_col_binder <- TRUE
    }

    #  check @batch_c$binder / $sinfo
    for(k in c("binder", "sinfo")) {
      if(xs1@batch_c[[k]] != x_j@batch_c[[k]])
        stop("All @batch_c$", k, " must be identical.")
    }
    #  check consistency of names in @assy_b
    if(! identical(names(xs1@assy_b), names(x_j@assy_b)))
      stop("names([", j, "]@assy_b) is not identical to ",
           "names([1]@assy_b) !")
    #  check consistency of names in @assy_b
    if(! identical(names(xs1@assy_s), names(x_j@assy_s)))
      stop("The names of @assy_s of all BAfs should be identical!")
  }
  
  
  ###  New @binder -----------------------------------------------------------

  binder_new <- xs1@binder
  
  
  ##  When @binder includes varying columns ----------------
  if(has_varying_col_binder) {
    
    ##  Find varying variables
    i_vary <- names(xs1@binder) %>% 
      map_lgl(function(ii) {
        map_lgl(2:length(xs), . %>% {
          ! identical(xs1@binder %>% select_at(ii),
                      xs[[.]]@binder %>% select_at(ii))
        }) %>% 
          any()       # any difference
      }) %>% 
      { names(xs1@binder)[.] }
    
    ##  Combine binder information for same binders
    binder_new[i_vary] <- sapply(i_vary, function(ii) {
        #  extract one varying column and bind_cols
        out <- do.call("bind_cols", map(xs, . %>% sB() %>% select_at(ii))) %>% 
          #  combine row wise
          apply(1, . %>% 
                  unlist() %>% 
                  combine_multi_elements()
          ) 
        if(is.factor(xs1@binder[[ii]])) factor(out) else out
      }, simplify= F) %>% 
      as_tibble
    
    
    ##  Fix binder batch IDs if this is one of the varying variables
    if(batch_colname(xs1, "b") %in% i_vary) {
      for(ii in seq_along(xs)) {
        batch(xs[[ii]], "binder") <- binder_new[[batch_colname(xs1, "b")]]
      }
    }
  }
  
  

  ###  New other slots --------------------------------------------------------

  ##  new @assy_s
  assy_s_new <- lapply(seq_along(xs1@assy_s), function(ii) {
    do.call("rbind", lapply(xs, . %>% {.@assy_s[[ii]]}))
  }) %>% 
    `names<-`(names(xs1@assy_s))
  
  ##  new @assy_b
  assy_b_new <- assy_new_rm_duplicates(xs, "assy_b")

  names(xs) <- NULL		# to avoid 1.r1,... 2.r20 for "sample" by rbind.data.frame
  
  new.Data <- do.call("rbind", lapply(xs, function(ea) ea@.Data))
  
  BAf(new.Data, 
      sinfo = do.call("rbind", lapply(xs, . %>% {.@sinfo})),
      binder = binder_new,
      sinfo_batch_i  = xs1@batch_c$sinfo, 
      binder_batch_i = xs1@batch_c$binder,
      assay_sinfo  = assy_s_new,
      assay_binder = assy_b_new,
      codebook_sinfo  = codebook_bind(xs, "sinfo"), 
      codebook_binder = codebook_bind(xs, "binder"), 
      note = do.call("c", lapply(xs, slot, "note"))
  )
}
