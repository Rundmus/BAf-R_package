# -----------------------------------------------------------------------------#
#' Signal distribution across binders
#' 
#' Depict a plot that displays signal distribution across samples
#' 
#' @param baf an object of the \code{\link{BAf-class}}
#' @param color_tbl_b the colors of the binders employed for positive or 
#'   negative control on plot. The information should be provided as a 
#'   \code{\link{data.frame}}, in which columns are \code{name}, \code{reg}
#'   (regular expression to find the binder), \code{col} (color), and
#'   \code{type} (\code{positive} , \code{negative}, or \code{measure}). If this
#'   is given as NULL, it results in no coloring.
#' @param guide_line if some guide lines of controls are displayed.
#' @param main,xlab,ylab,... same as \code{\link{boxplot}}
#' 
#' @return a \code{list} of \code{color_tbl_b}
#' 
#' @examples
#' data(sba)
#' plot_QC_binder_signal_boxplot(sba)
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' 
#' @include plot_QC_0_underlying.R
#' @export
# -----------------------------------------------------------------------------#
# created  : 2014-11-24 by Mun-Gwan
# modified : 
#   2017-08-20 by Mun-Gwan : to adapt to BAf-class
# -----------------------------------------------------------------------------#

plot_QC_binder_signal_boxplot <- 
  function(baf,
           color_tbl_b = NULL,
           main = deparse(substitute(baf)),
           xlab = "Binders",
           ylab = "Signal Intensity",
           guide_line = TRUE,
           ...) {
    
    stopifnot(inherits(baf, "BAf"))
    force(main)       # avoid lazy evaluation after change of 'baf'
    
    baf <- replace_0(baf, value= 1, show_count= F)
    
    c_l <- .boxplot_signal_color_tbl(
      x= sX(baf), ids= bid(baf), color_tbl= color_tbl_b, 
      main = main, xlab = xlab, ylab = ylab, ...
    ) 
    
    ##  Guide lines to show the levels of control binders
    if(!is.null(c_l$used_col) && nrow(c_l$used_col) > 0 && guide_line) {
      cTbl <- c_l$used_col %>%
        mutate(light= Useful2me::lighten_color(col, 100))
      
      for(ii in 1:nrow(cTbl)) {
        sX(baf)[, c_l$fg_col == cTbl$col[ii], drop= F] %>% 
          apply(., 2, median, na.rm= T) %>% 
          abline(h= ., col= cTbl$light[ii])
      }
    }
    invisible(list("color_tbl_b"= c_l$used_col))
  }

# -----------------------------------------------------------------------------#
#' @describeIn plot_QC_binder_signal_boxplot when \var{baf} has the data from
#'   suspension bead arrays
#' @export
# -----------------------------------------------------------------------------#
# created  : 2017-08-29 by Mun-Gwan
# modified : 
# -----------------------------------------------------------------------------#

plot_QC_binder_signal_boxplot.SBA <- 
  function(baf,
           color_tbl_b = default_color_SBA_bead(T),
           main = paste(unique(batch(baf, "binder")), collapse= ", "),
           xlab = "Beads",
           ylab = "Raw MFI",
           guide_line = TRUE,
           ...) {
    plot_QC_binder_signal_boxplot(
      baf = baf, color_tbl_b = color_tbl_b, 
      main = main, xlab = xlab, ylab = ylab, 
      guide_line = guide_line, ...)
  }


# -----------------------------------------------------------------------------#
#' Correlation with a binder
#' 
#' This depicts the distribution of correlation between a binder and all others.
#' 
#' @param baf an object of the \code{\link{BAf-class}}
#' @param b_id one binder IDs in numeric or in character with repeat indicator
#'   if necessary.
#' @param incl whether each column is included or not. This can be used to
#'   exclude other control binders
#' @param method,use refer to \code{\link{cor}}
#' @param main,xlab,... pass over to \code{\link{hist}}
#' @param show_count whether the count is displayed or not.
#' 
#' @return a list having two elements, computed correlation values and the
#'   output of \code{\link{hist}}
#' 
#' @examples
#' data(sba)
#' plot_QC_binder_corr_hist(sba, "rabbit IgG")
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' 
#' @include plot_QC_0_underlying.R
#' @export
# -----------------------------------------------------------------------------#
# created  : 2017-08-28 by Mun-Gwan
# modified :
# -----------------------------------------------------------------------------#

plot_QC_binder_corr_hist <- 
  function(baf, 
           b_id, 
           incl = rep(T, ncol(baf)),
           method = c("spearman", "pearson", "kendall"),
           use = "complete.obs",
           xlab = paste("Corr. with", b_id),
           main = paste("Histogram of", xlab),
           ...,
           show_count = TRUE) {
    
    stopifnot(inherits(baf, "BAf"),
              is.numeric(b_id) || is.character(b_id),
              length(b_id) == 1,
              is.logical(incl),
              length(incl) == ncol(baf))
    
    method <- match.arg(method)
    
    force(xlab)
    
    if(is.character(b_id)) {
      stopifnot(b_id %in% bid(baf, exact= TRUE))
      b_id <- which(bid(baf, exact= T) == b_id)   # to index
    }
    
    incl[b_id] <- F
    
    mbaf <- as.matrix(baf) %>% {
      cbind(.[, b_id, drop = F], .[, incl, drop= F])
    }

    corr <- stats::cor(mbaf, method= method, use= use)[-1, 1]
    
    for(ii in which(!incl)) {
      corr <- append(corr, NA, ii - 1)
      names(corr)[ii] <- bid(baf, exact= T)[ii]
    }

    out <- hist(corr, 
                breaks= seq(-1, 1, by= 0.1),
                xlab= xlab,
                main= main, 
                ...)
    out$xname <- xlab
    
    n <- out$counts
    text(out$mids[n > 0], 
         y= 0, 
         labels= n[n > 0], 
         adj= c(0.5, -1),
         col= "coral2",
         cex= 0.8)
    
    return(invisible(list("corr"= corr, "hist"= out)))
  }

