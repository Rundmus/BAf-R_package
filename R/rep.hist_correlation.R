# -----------------------------------------------------------------------------#
#' Histogram of correlation
#' 
#' Histogram of correlation values between repeated measures, together with a
#' density plot of those between random pairs of measures
#' 
#' @param baf an object of the \code{\link{BAf-class}} or a matrix
#' @param rep_binder a \code{factor} (or a \code{vector}) that indicates which
#'   binders were repeated, so their values should be combined. 
#' @param method,use same as those of \code{\link{cor}}, but different defaults
#' @param ... parameters that passed over to \code{\link{hist}} 
#' 
#' @return a \code{list} that has those correlation values
#' 
#' @seealso 
#' \code{\link{cor}}
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @importFrom Useful2me overwrite_par toupper_1st_char shuffle
#' @export
# -----------------------------------------------------------------------------#
# created  : 2016-10-18 by Mun-Gwan
# modified : 2017-08-22 by Mun-Gwan : to adapt to BAf-class
# -----------------------------------------------------------------------------#

hist_correlation <- function(baf,
                             rep_binder,
                             method = c("spearman", "pearson", "kendall"),
                             use= "complete.obs",
                             ...) {
  
  ##  Check given arguments
  stopifnot(inherits(baf, "BAf") || is.matrix(baf),
            is.vector(rep_binder) | is.factor(rep_binder),
            ncol(baf) == length(rep_binder))
  if(!is.factor(rep_binder)) rep_binder <- as.factor(rep_binder)
  method <- match.arg(method)
  
  n_rep <- table(rep_binder)
  if(!all(n_rep <= 2)) 
    stop("Currently, only the data sets with duplicated measures are allowed.")
  
  # list of indicies of the duplicated binders
  i_rep_binder <- levels(rep_binder) %>% 
    .[. %in% names(n_rep)[n_rep > 1]] %>%         # no single binder
    sapply(function(x) which(rep_binder == x), simplify= F)    # index
  
  if(inherits(baf, "BAf")) {
    #  batch
    ach <- batch(baf, "binder")
    baf <- sX(baf)
  } else {
    #  batch
    ach <- rep("B", ncol(baf))
  }

  ##  check if NA in a batch are given sample-wise
  if(use != "pairwise.complete.obs") {
    has_sample_wise_NA <- ach[unlist(i_rep_binder)] %>% 
      unique() %>%
      sapply(function(eb) {           # per binder batch
        tmp <- is.na(baf[, ach == eb])
        all(apply(tmp, 2, identical, tmp[, 1])) 
      })
    if(!all(has_sample_wise_NA))
      warning(
        "The sample with NA for a binder was not same as", 
        " for another binder in a batch [", 
        paste(names(has_sample_wise_NA)[has_sample_wise_NA], collpase= ","),
        "]."
      )
  }
  
  all_corr <- baf %>% 
    stats::cor(., method= method, use= use) %>% 
    `diag<-`(NA)
  
  ##  **  Correlation between repeated pairs  **  ##
  corr <- sapply(i_rep_binder, function(ii) {
    all_corr[ii, ii] %>% 
    {.[upper.tri(.)]}
  })
  
  ##  **  Correlation with a random pair in same assay as partner  **  ##
  rnd_corr <- sapply(i_rep_binder, function(ii) {
    #  Same batch as partner but different binder
    which(
      ach == ach[ii][2] & !unwhich(ii, length(ach))
    ) %>% {
      # randomly selected 10 or all others
      if(sum(.) > 10) shuffle(., 10) else .
    } %>% 
      all_corr[ii[1], .]
  })
  
  #  Default parameters for 'hist'
  hist_par <- list(
    x = corr,
    main = "Duplicated assays", xlim= c(0, 1), 
    xlab= paste0(toupper_1st_char(method), "'s correlation")
  )
  dots <- if(missing(...)) list() else list(...)
  hist_par <- overwrite_par(prev_par= hist_par, new_par= dots, excl= "x")
  
  #  ***    Histogram    ***  #
  # suppress warnings for the case 'plot= F'
  out_hist <- suppressWarnings(do.call("hist", hist_par)) 
  
  if(dev.cur() > 1) {
    ## *   Density plot of random correlation   * ##
    density(as.vector(unlist(rnd_corr))) %>% {
      .$y <- .$y * (par()$usr[4] * 0.95) / max(.$y)       # adjust y range
      .
    } %>% 
      lines(col= "cadetblue")
  }
  
  return(invisible(list(corr= corr, rnd_corr= rnd_corr, hist= out_hist)))
}
