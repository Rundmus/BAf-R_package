# -----------------------------------------------------------------------------#
#' Color by category
#' 
#' This returns color codes grouped by categories
#' 
#' @param x a \code{vector} of a variable
#' @param cat the categories in the \code{x} by which colors should be given
#' @param color the color for each \code{cat}egory
#' @param def_color the color for no-match
#' @param grep if T, then this function takes \code{cat} as regular expressions.
#'   The catergories are identified after changing \code{x} to lower cases and
#'   removing non-alphanumeric characters.
#' @return a \code{\link{vector}} of colors
#' @examples
#' data(sba)
#' summary(farg(sba@sinfo$plate))
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @noRd
# -----------------------------------------------------------------------------#

farg <- function(x,
                 cat = NULL,
                 color = NULL,
                 def_color = "black",
                 grep = FALSE) {
  if(length(cat) == 0) return(rep(def_color, length(x)))
  
  stopifnot(length(cat) <= length(color),
            inherits(x, "character") | inherits(x, "factor"))
  x <- as.factor(x)
  
  if(grep) {
    x <- .fit_to_std_format(x)    # prepare for regular expression
    out <- vapply(cat, grepl, rep(T, length(x)), x)
  } else {
    out <- vapply(cat, function(ec) {
      x == ec
    }, rep(T, length(x)))
  }
  
  
  out <- ( if(grep) {
    x <- .fit_to_std_format(x)    # prepare for regular expression
    vapply(cat, grepl, rep(T, length(x)), x)
  } else {
    vapply(cat, function(ec) {
      x == ec
    }, rep(T, length(x)))
  } ) %>%
    apply(1, function(eR) {
      if(any(eR)) color[eR][1] else def_color
    })
  
  return(out)
}



# -----------------------------------------------------------------------------#
#' Get color table
#' 
#' Get a table to color control binders / samples. If \code{x} is missing then
#' \code{def} is used.
#' 
#' @param x a data frame having a color table which have at least "name", "reg",
#'   and "col" columns
#' 
#' @return a data frame having a color table
#' @noRd
# -----------------------------------------------------------------------------#

.chq_color_tbl <- function(x) {
  if(is.null(x)) {
    return(NULL)
  } else {
    stopifnot(inherits(x, "data.frame"),
              all(c("name", "reg", "col") %in% colnames(x)))
    return(x %>% mutate_if(is.factor, as.character))
  }
}


# -----------------------------------------------------------------------------#
#' reorder sample and replace 0
#' 
#' Reorder the samples of an BAf according to given \code{sample_order} and
#' replace 0 with 1 to avoid any error in log-scale plots
#' 
#' @param baf an object of the \code{\link{BAf-class}}
#' @param sample_order the order of sample in plots. The default is the order by
#'   sample batch in order to easily check batch differences.
#' @param show_count_of_0 whether the number of 0 should be displayed
#' @return an object of the \code{\link{BAf-class}}
# -----------------------------------------------------------------------------#

.reorder_BAf_and_replace_0 <-
  function(baf,
           sample_order = order(batch(baf, "sample")),
           show_count_of_0 = FALSE) {
    baf[sample_order,] %>% 
      replace_0(., value= 1, show_count = show_count_of_0, 
                by_s= batch_colname(baf, "sample"))
  }



# -----------------------------------------------------------------------------#
#' boxplot using color table
#' 
#' Find colors from a given color table and draw a box-and-whisker plot
#' 
#' @param x a matrix 
#' @param ids a vector of IDs
#' @param color_tbl a \code{tbl_df} having coloring variables. All
#'   \code{'reg'}, \code{col}, and \code{name} are required.
#' @param grep same as \code{\link{farg}}
#' @param main,xlab,ylab,... same as \code{\link{boxplot}}
#' @return a list of \code{fg_col} and \code{used_col} (a subset of
#'   \var{color_tbl} of the colors used here)
# -----------------------------------------------------------------------------#

.boxplot_color_tbl <- function(x, ids, color_tbl, grep= FALSE, 
                             main, xlab, ylab, ...) {
  # foreground color
  fg_col <- farg(ids, color_tbl$reg, color_tbl$col, grep = grep)
  
  col0 <- Useful2me::lighten_color(fg_col, 70)
  
  boxplot(x, log= "y", cex= 0.5, 
          outcol= fg_col, boxcol= fg_col, 
          outpch= ifelse(fg_col == "black", 1, 16), 
          col= col0, outbg= col0,
          main= main, xlab= "", ylab= ylab,        # xlab, its own 'mgp'
          xaxt= "n", ...)
  
  title(xlab= xlab, mgp= c(1, 1, 0))
  
  used_col <- color_tbl[color_tbl$col %in% fg_col, ]
  
  if(!is.null(used_col) && nrow(used_col) > 0)      # if any special color was used,
    legend("bottomright", legend= used_col$name, #  show a legend
           text.col= used_col$col, bty= "n", cex= 0.5)
  return(list(fg_col= fg_col, used_col= used_col))
}


# -----------------------------------------------------------------------------#
#' boxplot of signal using color table
#' 
#' Draw a box-and-whisker plot that displays signal distribution. The varisous
#' colors of box and outlier points were used to distinguish special groups
#' 
#' @inherit .boxplot_color_tbl params return
#' 
#' @include plot_QC_0_underlying.R
#' @noRd
# -----------------------------------------------------------------------------#
.boxplot_signal_color_tbl <- function(x, ids, color_tbl, 
                                      main, xlab, ylab, ...) {
  # color_tbl
  color_tbl <- .chq_color_tbl(color_tbl) 
  
  opar <- par(mar= c(3, 5, 4, 2) + 0.1, ask= dev.interactive())
  on.exit(par(opar))
  
  .boxplot_color_tbl(x= x, ids= ids, color_tbl= color_tbl, grep= TRUE,
                     main= main, xlab= xlab, ylab= ylab, ...) 
}


