# -----------------------------------------------------------------------------#
#' Signal distribution across samples
#' 
#' Draw a plot that displays signal distribution across samples
#' 
#' @inheritParams .reorder_BAf_and_replace_0
#' @param color_tbl_s colors of samples on the plot. The information should be 
#'   provided as a \code{\link{data.frame}} that has \code{name}, \code{reg} 
#'   (regular expression to find the sample), \code{col} (color), \code{type} 
#'   (\code{negative} or \code{replicated}).\cr If this is NULL, then no
#'   coloring is applied.
#' @param main,xlab,ylab,... same as \code{\link{boxplot}}
#' 
#' @return a \code{list} of the "\var{color_tbl_s}"
#' 
#' @examples
#' data(sba)
#' plot_QC_sample_signal_boxplot(sba)
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @include plot_QC_0_underlying.R
#' @export
# -----------------------------------------------------------------------------#
# created  : 2014-11-24 by Mun-Gwan
# modified : 
#   2017-08-20 by Mun-Gwan : to adapt to BAf-class
# -----------------------------------------------------------------------------#

plot_QC_sample_signal_boxplot <- 
  function(baf,
           sample_order = order(batch(baf, "sample")),
           color_tbl_s = NULL,
           main = deparse(substitute(baf)),
           xlab = "Sample (sorted)",
           ylab = "Signal Intensity",
           ...) {
    stopifnot(inherits(baf, "BAf"))
    force(main)       # avoid lazy evaluation after change of 'baf'
    
    baf <- .reorder_BAf_and_replace_0(baf, sample_order, F)
    
    .boxplot_signal_color_tbl(
      x= t(sX(baf)), ids= sid(baf), color_tbl= color_tbl_s, 
      main = main, xlab = xlab, ylab = ylab, ...
    ) %>% {
      invisible(list("color_tbl_s"= .$used_col))
    }
  }

# -----------------------------------------------------------------------------#
#' @describeIn plot_QC_sample_signal_boxplot when \var{baf} has the data from
#'   suspension bead arrays
#' @export
# -----------------------------------------------------------------------------#
# created  : 2017-08-29 by Mun-Gwan
# modified : 
# -----------------------------------------------------------------------------#

plot_QC_sample_signal_boxplot.SBA <- 
  function(baf,
           sample_order = order(batch(baf, "sample")),
           color_tbl_s = default_color_SBA_ctrl_samples(sid(baf)),
           main = paste(unique(batch(baf, "binder")), collapse= ", "),
           xlab = "Sample (sorted)",
           ylab = "Raw MFI",
           ...) {
    plot_QC_sample_signal_boxplot(baf = baf, sample_order = sample_order,
      color_tbl_s = color_tbl_s, main = main, xlab = xlab, ylab = ylab, ...)
  }
