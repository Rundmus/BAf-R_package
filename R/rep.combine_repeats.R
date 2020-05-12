# -----------------------------------------------------------------------------#
#' Function to combine repeated measures
#' 
#' Merge repeated measures by taking one representative value (e.g. average) of 
#' them. It assumes the measurements were conducted in separate binder batches
#' (or assays), so that the data are stored in different columns of \code{baf}.
#' If repeated measurements were conducted in a same binder batch, transform the
#' data using \code{\link{transform_to_landscape_for_repeats}}.
#'
#' @param baf an object of the \code{\link{BAf-class}}
#' @param rep_binder a \code{factor} (or a \code{vector}) that indicates which
#'   binders were repeated, so their values should be combined. 
#' @param rm_grp_rep_binder if the \code{'grp'} and \code{'rep_binder'} columns
#'   in \code{@binder} that were added by 
#'   \code{\link{transform_to_landscape_for_repeats}}, they should be removed.
#' @param center,scale if the difference between assays should be adjusted using
#'   \code{\link{scale}}
#' @param FUN the function to compute a reference value such as mean and median 
#' @param ... not used
#'
#' @return
#' A \code{\link{BAf-class}} in which the two observed values were replaced with
#' the representative value computed by the \code{FUN}.
#' 
#' @seealso 
#' \code{\link{transform_to_landscape_for_repeats}}
#'
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @export
#' @importFrom Useful2me combine_multi_elements
# -----------------------------------------------------------------------------#
# created  : 2016-10-19 by Mun-Gwan
# modified : 2017-08-22 by Mun-Gwan : to adapt to BAf-class
# -----------------------------------------------------------------------------#

combine_repeats <-
  function(baf,
           rep_binder,
           rm_grp_rep_binder = T,
           center = T,
           scale = F,
           FUN = function(x) mean(x, na.rm = T),
           ...) {
    
    ##  Check given arguments
    stopifnot(inherits(baf, "BAf"),
              is.vector(rep_binder) | is.factor(rep_binder),
              dim(baf)[2L] == length(rep_binder))
    
    ##  Make sure [rep_binder] is a factor, in which levels are in the same
    ##  order of appearances
    rep_binder <- as.character(rep_binder) %>%
      factor(., levels= unique(.))
    
    tmp <- table(rep_binder)
    if(!all(tmp == 2)) 
      stop("Currently, only the data sets with duplicated measures are allowed.")
    
    #  index for each repeated binder
    i_by_b <- by(1:ncol(baf), rep_binder, c)
    
    ##  adjustement for the difference between assays using scale
    sX_adj <- scale(sX(baf), center= center, scale= scale)
    
    scale_attr <- function(x, var, def) {
      scaled_center <- attr(sX_adj, paste0("scaled:", var)) %>% 
      {    
        if(is.null(.)) rep(def, length(i_by_b)) 
        else sapply(i_by_b, function(ii) FUN(.[ii]))
      }
    }
    scaled_center <- scale_attr(sX_adj, "center", 0)
    scaled_scale  <- scale_attr(sX_adj, "scale",  1)
    
    
    ##  ***   Compute   ***  ##
    X <- sapply(i_by_b, function(ii) {
      apply(sX_adj[, ii], 1, FUN) %>% {
        if(!is.vector(.)) stop("The FUN produced a non-vector object.")
        else .
      } %>% 
        `[<-`(is.na(.), NA)    # NaN
    }) %>% 
      t() %>% 
      {. * scaled_scale + scaled_center} %>%
      t()
    
    ##    Binder    ##
    binder <- lapply(i_by_b, function(ii) {
      runs <- sB(baf)[ii, ] %>% 
        mutate(id = rm_repeat_indicator(id))
      if(rm_grp_rep_binder) {
        iGrpCol <- grep("^grp\\.?[[:digit:]]*$", colnames(sB(baf)))
        iRepCol <- grep("^rep_binder\\.?[[:digit:]]*$", colnames(sB(baf)))
        iRmCol <- c(max(iGrpCol, 0), max(iRepCol, 0))
        iRmCol <- iRmCol[iRmCol != 0]
        runs <- if(length(iRmCol) == 0) runs else runs[, -c(iRmCol)]
      }
      runs %>% summarise_all(Useful2me::combine_multi_elements)
    }) %>% 
      do.call("rbind", .)
    
    ##  To let factor variables have same levels  
    for(ii in 1:ncol(binder)) {
      if(is.factor(binder[, ii])) {
        o_sB <- sB(baf)[unclass(rep_binder)[!duplicated(rep_binder)], ii]
        if(all.equal(as.character(binder[, ii]), as.character(o_sB)) == T) {
          binder[, ii] <- o_sB
        }
      }
    }
    
    batch_b <- sapply(i_by_b, function(ii) batch(baf, "binder")[ii])
    ##  binder batch IDs merging
    batch_b_merge <- apply(batch_b, 2, Useful4me::combine_multi_elements) %>%
      factor(., levels= unique(.))
    
    ##  @assy_s
    assy_s_new <- lapply(sA(baf, "sample"), function(eA) {
      lapply(levels(batch_b_merge), FUN= function(ii) {
        tmp <- batch_b[, batch_b_merge == ii][, 1]
        sapply(1:nrow(eA), function(jj) {
          Useful2me::combine_multi_elements(unlist(eA[jj, unique(tmp)]))
        })
      }) %>%
        as.data.frame() %>%
        `colnames<-`(levels(batch_b_merge)) %>%
        `rownames<-`(rownames(eA))
    })
    
    ##  @assy_b
    assy_b_new <- lapply(sA(baf, "binder"), function(eA) {
      lapply(i_by_b, function(ii) {
        summarise_all(eA[ii, ], Useful2me::combine_multi_elements)
      }) %>% 
        do.call("rbind", .)
    })

    BAf(X, 
        sinfo = baf@sinfo,		# note NO added row or column
        binder = binder,
        sinfo_batch_i  = baf@batch_c$sinfo, 
        binder_batch_i = baf@batch_c$binder,
        assay_sinfo  = assy_s_new,
        assay_binder = assy_b_new,
        codebook_sinfo  = codebook(baf, "sinfo"), 
        codebook_binder = codebook(baf, "binder"), 
        note = note(baf)
    )
  }
