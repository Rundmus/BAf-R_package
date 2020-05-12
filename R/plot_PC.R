# -----------------------------------------------------------------------------#
#' Plot of principal components
#' 
#' The plot of a pair of principal components is displayed for a matrix or a
#' single \code{\link{BAf-class}} object. Colors or symbols of the points can be
#' diffrentiated by the variable given in \var{col_by} or \var{pch_by} to assist
#' to find hidden substructure of sample. The principal components were obtained
#' by (\code{\link{prcomp}}) with scaling.
#' 
#' @param x an object of the \code{\link{prcomp}}, \code{\link{matrix}} or
#'   \code{\link{BAf-class}}
#' @param pc a vector that consists of two integer elements for two principal
#'   components
#' @param main,xlab,ylab same as \code{\link{plot.default}}
#' @param ... pass to \code{plot}
#' 
#' @return \code{plot_PC} \code{\link{prcomp}} output object
#' 
#' @author Mun-Gwan Hong \email{mun-gwan.hong@scilifelab.se}
#' @seealso \code{\link{prcomp}}
#' @examples 
#' data(sba)
#' plot_PC(sba)
#' plot_PC(sba, col_by= "plate")    # no legend
#' plot_PC(sba, col_by= "plate", pch_by= "sex", legend_arg= list())
#' 
#' @aliases plot_PC
#' @export
# -----------------------------------------------------------------------------#
# created  : 2012-05-02 by Mun-Gwan
# modified : 
#   2013-07-03 by Mun-Gwan : adapt to "SBAe" class
#   2015-05-29 by Mun-Gwan : add percentage of variances
#   2017-08-22 by Mun-Gwan : fail if NA, split by input class
# -----------------------------------------------------------------------------#

setGeneric("plot_PC", function(x, pc = c(1, 2), main = "", xlab, ylab, ...) {
  standardGeneric("plot_PC")
})

# -----------------------------------------------------------------------------#
#' @describeIn plot_PC The input is the output of 'prcomp'
# Commented out due to warning message: 
#   class "prcomp" is not exported by 'namespace:stats' 
# #' @importClassesFrom stats prcomp   
#' @export
# -----------------------------------------------------------------------------#
setMethod(
  "plot_PC", signature(x = "prcomp"), 
  function(x, pc = c(1, 2), main = "", xlab, ylab, ...) {
    
    stopifnot(max(pc) < ncol(x$x))
    pca <- x
    
    paste_percent <- function(pc_i, var1) {
      paste0("PC", pc_i, " (", round(var1[pc_i] / sum(var1) * 100, 1), "%)")
    }
    varvec <- pca$sdev^2
    if(missing(xlab)) xlab <- paste_percent(pc[1], varvec)
    if(missing(ylab)) ylab <- paste_percent(pc[2], varvec)
    
    plot(pca$x[, pc[1:2]], main = main, xlab= xlab, ylab= ylab, ...)
    
    return(invisible(pca))
  }
)


# -----------------------------------------------------------------------------#
#' @describeIn plot_PC The input is a 'matrix'
#' @export
# -----------------------------------------------------------------------------#
setMethod(
  "plot_PC", signature(x = "matrix"), 
  function(x, pc = c(1, 2), main = "", xlab, ylab, ...) {
    stopifnot(max(pc) < ncol(x))
    #*  PCA with scaling  *#
    pca <- na.fail(x) %>% prcomp(scale. = TRUE)
    
    plot_PC(pca, main = main, xlab= xlab, ylab= ylab, ...)
    
    return(invisible(pca))
  }
)

# -----------------------------------------------------------------------------#
#' @describeIn plot_PC The input is a 'BAf' object
#' 
#' @param col_by,pch_by a name of a column in \code{sI(x)}. Different color /
#'   symbols the points will have.
#' @param legend_arg a list of parameters to be used in adding a legend. If
#'   NULL, skip to add legend
#' 
#' @importFrom Useful2me overwrite_par
#' @export
# -----------------------------------------------------------------------------#
setMethod(
  "plot_PC", signature(x = "BAf"), 
  function(x,
           pc = c(1, 2),
           main = "",
           xlab,
           ylab,
           ...,
           col_by = NULL,
           pch_by = NULL,
           legend_arg = list(x= "topright")) {
    stopifnot(max(pc) < ncol(x))
    
    dots <- list(...)
    
    by2 <- list(col= col_by, pch= pch_by)
    has_by <- lapply(by2, . %>% {!is.null(.)})

    for(ii in names(by2)) {
      if(has_by[[ii]]) {
        e_by <- by2[[ii]]
        if(!e_by %in% colnames(x@sinfo))
          stop("'", ii, "_by' cannot be found in '@sinfo'.")
        if(!is.null(dots[[ii]])) warning("The '", ii, "' is ignored.")
        
        # by2 is possibly changed here.
        by2[[ii]] <- x@sinfo[[e_by]] %>% as.factor()
        dots[[ii]] <- by2[[ii]] %>% as.integer()
      } else {
        if(is.null(dots[[ii]])) dots[[ii]] <- rep(par()[[ii]], nrow(x))
      }
    }

    #*  plot_PC  *#
    plot_skip_col_pch <- function(..., col, pch, col_, pch_) {
      plot_PC(..., col= col_, pch= pch_)
    }
    pca <- plot_skip_col_pch(
      x= sX(x), main = main, xlab= xlab, ylab= ylab, ..., 
      col_= dots$col, pch_= dots$pch
    )

    ##  add legend
    has_by <- unlist(has_by)
    if(!is.null(legend_arg) && any(has_by)) {

      lg <- c(by2[has_by], dots[names(by2)]) %>% 
        do.call("data.frame", .) %>% 
        na.omit() %>% 
        mutate(lgd = apply(.[, c(1:2)[has_by], drop= F], 1, 
                           paste, collapse= ",")) %>% 
        distinct()
      
      # the legend parameters dependent to the plot above
      list(
        x = "topright",
        legend = lg$lgd,
        col    = lg$col.1,
        text.col = lg$col.1,
        pch    = lg$pch.1,
        bty    = "n"
      ) %>% 
        Useful2me::overwrite_par(new_par= legend_arg) %>% 
        do.call("legend", .)
    }
    
    return(invisible(pca))
  }
)
