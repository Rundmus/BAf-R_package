# -----------------------------------------------------------------------------#
#' Plot to compare a pair of vectors
#' 
#' The plot function to compare a pair of vectors. One of the main usages is to
#' check how homogeneous is the data generated on multiple plates. The points
#' are colored by group.
#' 
#' @param x a vector of numerical values
#' @param y a vector of values that are very likely correlated with values in 'x'
#' @param grp_id a vector that has the group IDs of which each will be plotted
#'   in different color\cr - It should have same length as \code{x} and \code{y}
#' @param square whether the ranges of the x- and y-axes in plot should same or
#'   not
#' @param log which coordinates are in log scale. Same as \code{'log'} in 
#'   \code{\link{plot.default}}.
#' @param main,xlab,ylab same as the ones for \code{'plot'} generic function 
#'   (refer \code{\link{par}})
#' @param col_ea_grp color of each group given in \code{'grp_id'}
#' @param ... any arguments that pass to \code{'plot'} and \code{'points'}
#' @param legend_arg a list of parameters to be used in adding a legend. If 
#'   \code{grp_id} is missing then no legend will be displayed.
#'   
#' @return A plot is generated and a \code{"list"} of x- and y- coordinate
#'   limits is returned.
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @export
# -----------------------------------------------------------------------------#
# created  : 2011-07-06 by Mun-Gwan
# modified : 
#   2011-09-16 by Mun-Gwan : 
#     1) replace 'xPosRel', 'yPosRel' with 'legend_arg' to let the function more
#        flexible
#     2) replace many arguments to '...'
#   2016-10-17 by Mun-Gwan : 
#     1) change function name. 
#     2) Remove 'xlim', 'ylim'
#     3) add 'col_ea_grp'
#   2017-08-23 by Mun-Gwan : to adapt to BAf-class
# -----------------------------------------------------------------------------#

plot_a_pair_vectors <-
  function(x,
           y,
           grp_id,
           square = TRUE,
           log = "xy",
           main = paste("Scatter plot of", xlab, "&", ylab),
           xlab = deparse(substitute(x)),
           ylab = deparse(substitute(y)),
           col_ea_grp,
           ...,
           legend_arg = list(x = "bottomright")) {
    
    # check if the necessary data vectors are provided
    stopifnot(!missing(x) & !missing(y),
              is.numeric(x), 
              is.numeric(y),
              length(x) == length(y))			# check the length of two vectors
    
    lth <- length(x)		# the length of the vector 'x'
    
    ## if 'grp_id' is missing, no legend
    if(missing(grp_id)) {
      nGroupID <- 1
    } else {
      stopifnot(length(x) == length(grp_id))		
      # NO sorting of levels
      grp_id <- as.factor(grp_id)
      nGroupID <- nlevels(grp_id)	# how many varieties of group IDs
    }

    ###   limits of x and y axes in plot
    # change the names of arguments to full names in 'plot.default'
    arg_plot <- match.call(plot.default)
    xlim <- eval(arg_plot$xlim, parent.frame(n= 1))
    ylim <- eval(arg_plot$ylim, parent.frame(n= 1))
    if(square && !is.null(xlim) && !is.null(ylim)) {	## check conflict
      stopifnot(square == TRUE & identical(xlim, ylim))
    }
    if(square) {		# if square plot
      rg_y <- rg_x <- if(is.null(xlim)) { 
        if(is.null(ylim)) range(c(x, y)) else ylim 
      } else xlim
    } else {
      rg_x <- if(is.null(xlim)) range(x, na.rm= T) else xlim
      rg_y <- if(is.null(ylim)) range(y, na.rm= T) else ylim
    }
    
    ###   change PALETTE to have gradual color changes ###
    if(missing(col_ea_grp)) {
      if(nGroupID > 7) {
        opal <- palette(rainbow(nGroupID, start = 0.4, end = 0.1))
        on.exit(palette(opal))
      }
    } else {
      opal <- palette(if(is.numeric(col_ea_grp)) palette()[col_ea_grp] 
                      else col_ea_grp)
      on.exit(palette(opal))
    }
    
    
    #------------------------#
    # <<<<<    PLOT    >>>>> #
    #------------------------#
    
    # initialize the plot with an empty area in order to draw 'abline' first
    plot(rg_x, rg_y, log = log, main=main, xlab=xlab, ylab=ylab, type="n", ...)
    
    # add line of identity or basis line
    if(square) abline(0, 1, col="cadetblue1") else abline(h=0, col="cadetblue1")
    
    
    ###  Scatter points  ###
    dots <- list(...)
    if(!is.null(dots$col)) warning("The \'col\' is ignored.")
    # arguments for [points], ignore all parameters for [plot.default] and 'col'
    points_dots <- function(..., col, log, axes, frame.plot, 
                            panel.first, panel.last) {
      points(..., col= if(nGroupID == 1) 1 else unclass(grp_id))
    }
    # scatter plot
    points_dots(x, y, ...)
    
    
    ###   add LEGEND   ###
    
    if(! missing(grp_id) ) {
      stopifnot(class(legend_arg) == "list")
      
      # the legend parameters dependent to the plot above
      list(
        legend = levels(grp_id),
        col    = 1:nGroupID,
        text.col = 1:nGroupID,
        pch    = c(as.list(arg_plot), par())$pch,
        pt.cex = c(as.list(arg_plot), par())$cex
      ) %>% 
        ## over-write the default parameters over 'legend_arg'
        c(., legend_arg, list(bty= "n")) %>% 
        { .[!duplicated(names(.))] } %>% 
        do.call("legend", .)
    }
    
    # return the limit of x and y coordinates
    return(invisible(list(x = rg_x, y = rg_y)))
  }




# -----------------------------------------------------------------------------#
#' MA-plot with coloring by group
#' 
#' @param x1,x2 the first and second vectors to show 
#' @param grp_id a vector that has the group IDs of which each will be plotted
#'   in different color\cr - It should have same length as \var{x1} and \var{x2}
#' @inheritParams plot_a_pair_vectors
#' @param x1name,x2name the names of \var{x1} and \var{x2}. These will be used
#'   in the labels of x-, y-axis
#' @param log whether the \code{MA} should be computed in log-scale or not. If
#'   this is TRUE, 'x' coordinate is displayed in log scale
#' @param ... any arguments that pass to \code{plot_a_pair_vectors}
#' 
#' @return A plot is generated and a \code{"list"} of x- and y- coordinate
#'   limits is returned.
#'
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @export
# -----------------------------------------------------------------------------#
# required : plot_a_pair_vectors()
#
# created  : 2011-07-06 by Mun-Gwan
# modified : 
#   2017-08-23 by Mun-Gwan : to adapt to BAf-class
# -----------------------------------------------------------------------------#


plot_a_pair_vectors_in_MA <- 
  function(x1, x2, grp_id, 
           x1name = deparse(substitute(x1)), 
           x2name = deparse(substitute(x2)), 
           log = TRUE, main, xlab, ylab, ..., legend_arg) {
    
    stopifnot(!missing(x1) & !missing(x2))
    
    if(log) {
      x1 <- log2(x1)
      x2 <- log2(x2)
    }
    
    ## calculate M and A
    M2 <- (x1 + x2) /2
    A <- x2 - x1
    
    if(log) {
      M2 <- 2^ M2
    }

    ## default title and x-, y-axis labels
    if(missing(main)) main <- paste("MA-plot of", x1name, "&", x2name)
    if(missing(xlab)) {
      xlab <- paste0(if(log) "geometic" else "",
                     "mean of (", x1name, ") & (", x2name, ")")
    }
    if(missing(ylab)) {
      ylab <- if(log) {
        substitute(paste(log[2], "(", n2, ") - ", log[2], "(", n1, ")"),
                   list(n1=x1name, n2=x2name))
      } else {
        paste0("(", x2name, ") - (", x1name, ")")
      }
    }

    ### plot ###
    out <- plot_a_pair_vectors(
      x= M2, y= A, grp_id= grp_id, square = FALSE, 
      log= (if(log) "x" else ""), 
      main= main, xlab= xlab, ylab= ylab, 
      ..., legend_arg= legend_arg)

    return(out)
  }


# -----------------------------------------------------------------------------#
#' Two plots to check a pair of vectors
#' 
#' Two plots to check how homogeneous are the data in the pair of vectors. This 
#' function shows a scatter plot and MA-plot colored by group in a single
#' window.
#' 
#' @inheritParams plot_a_pair_vectors_in_MA
#' @param xlim same as that for \code{\link{plot}} function
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @noRd
# -----------------------------------------------------------------------------#
# required : plot_a_pair_vectors(), plot_a_pair_vectors_in_MA()
# created  : 2011-07-06 by Mun-Gwan
# modified :
#   2017-08-23 by Mun-Gwan : to adapt to BAf-class
# -----------------------------------------------------------------------------#

plot_a_pair_vectors_in_2_plots <- 
  function(x1, x2, grp_id, 
           x1name = deparse(substitute(x1)), 
           x2name = deparse(substitute(x2)), 
           log = TRUE, main, xlab, ylab, ..., legend_arg) {
    
    # layout to have two plots in a window, 
    # unless the layout has been modified by other function invoking this.
    if(par()$mfcol == c(1,1) && par()$mfrow == c(1,1)) { 
      opar <- par(mfcol = c(2,1)); on.exit(par(opar)) 
    }
    
    ## default names of variables and titles of the plots
    if(missing(main)) {
      main1 <- paste("Scatter plot of", x1name, "&", x2name)
      main2 <- paste("MA-plot of", x1name, "&", x2name)
    } else {
      # if the 'main' is provided, combind it to the type of plot for title
      main1 <- paste("Scatter plot of", main)
      main2 <- paste("MA-plot of", main)
    }
    
    arg_plot <- match.call(plot.default)
    xlim <- eval(arg_plot$xlim, parent.frame(n= 1))
    ylim <- eval(arg_plot$ylim, parent.frame(n= 1))
    if(!is.null(ylim)) warning("The given 'ylim' is ignored.")
    # let the two plots have same range of x coordinates
    if(is.null(xlim)) xlim <- range(x1, x2, na.rm= T)
    
    ## The both plot
    plot_a_pair_vectors(x1, x2, grp_id, square = TRUE, 
                        log= (if(log) "xy" else ""), 
                        main1, x1name, x2name, 
                        xlim= xlim, ylim= xlim, ..., legend_arg)
    rg <- plot_a_pair_vectors_in_MA(x1, x2, grp_id, 
                                    x1name, x2name, 
                                    log, main2, 
                                    xlim= xlim, ..., legend_arg)
    
    # return the limit of x and y coordinates
    return(invisible(list(x = xlim, y = rg$y)))
  }
