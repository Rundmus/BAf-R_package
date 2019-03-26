# -----------------------------------------------------------------------------#
# A collection of the functions that make a plot in different kinds 
# -----------------------------------------------------------------------------#





# -----------------------------------------------------------------------------#
#' To compare means in every plates
#' 
#' plot to observe correlation between the estimates (e.g. mean) that represent intensities for each
#' antibody in a single window coloring each assay separately.
#' 
#' @name Plot_of_mean_comparison
#' 
#' @param X a \code{\link{BAf}} object or a matrix in which each ROW represents each TARGET and each
#'   COLUMN indicates each COLLECTION of samples to be paired. The column name should indicate the
#'   name of the collections.
#' @param groupID a vector of the group ids. Each group will be plotted in different color. It
#'   should have same number of rows of \code{X}
#' @param geometric if the geometric mean is computed instead of arithmetic one.
#' @param by by which column the samples should be devided during calculation of the mean
#' @param log logically indicates if the coordinates are in log scale
#' @param MAplot whether an MA-plot is wanted
#' @param main same as the ones for \code{'plot'} generic function (refer \code{\link{par}})
#' @param legend_par a list of parameter to be used in adding a legend. If \code{groupID} is missing
#'   then no legend will be added.
#' @param prefix_plate the prefix added to the plate ID shown in X- and Y-axes
#' @param ... any arguments that pass to \code{\link{plot}} and \code{\link{points}}
#' 
#' @return 
#' A plot is generated and a \code{"list"} of x- and y- coordinate limits is returned.
#' @examples
#' data(sba)
#' \dontrun{plot_compareMeanOfEach_inPairs.BAf(sba)}
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @export
#' @aliases plot_compareMeanOfEach_inPairs
#' @rdname plot_compareMeanOfEach_inPairs
# -----------------------------------------------------------------------------#
# required : plot_a_pair_vectors()
# return : if(MAplot) 'list(x=xlim, y=ylim) of current plot' else 'xlim (=ylim)' 
#
# created  : 2011-05-23 by Mun-Gwan
# modified : 2011-07-06 by Mun-Gwan : use plot_a_pair_vectors instead of plotByAssay
#            2012-12-10 by Mun-Gwan : add title on outer margin
#            2013-07-03 by Mun-Gwan : adapt to "SBAe" class
# -----------------------------------------------------------------------------#

setGeneric("plot_compareMeanOfEach_inPairs", function(X, ...) 
    standardGeneric("plot_compareMeanOfEach_inPairs"));

#' @export
#' @rdname plot_compareMeanOfEach_inPairs
setMethod(
    "plot_compareMeanOfEach_inPairs", signature(X = "matrix"), 
    function(X, groupID, log = TRUE, MAplot = FALSE, main, legend_par= list(x= "bottomright"), ...) {
        
        if(!missing(groupID)) stopifnot(dim(X)[1L] == length(groupID))
        
        nCln <- dim(X)[2L]		# the number of collections
        
        ## multiple plots in a window
        opar <- par(mfrow = c(ceiling(nCln/2), nCln - 1), 
                    oma= c(0, 0, 5, 0), mar= c(5, 5, 0, 1)+0.1)
        on.exit(par(opar))
        
        ## labels
        lab <- c()
        for (i in 1:nCln) {
            lab[i] <- if(is.null(colnames(X)[i])) i else {	# the default is sequencial number
                tmp <- colnames(X)[i]
                paste0(toupper(substr(tmp, 1, 1)), substring(tmp, 2))
            }
        }
        
        pairs <- combn(nCln, 2) %>% as.data.frame %>% as.list
        if(MAplot) {
            ## plot limit
            rg_x <- range(sapply(pairs, function(j) { (X[,j[1]] * X[,j[2]])^(1/2) }))
            rg_y <- range(sapply(pairs, function(j) { log2(X[,j[1]]) - log2(X[,j[2]]) }))
            for(j in pairs) {
                plot_a_pair_vectors_in_MA(X[,j[2]], X[,j[1]], groupID, lab[j[2]], lab[j[1]], 
                                               log, "", , , rg_x, rg_y, legend_par, ...)
            }
            mtext(text= main, outer= TRUE, padj= -1, cex= 1.5, font= 2)
            # return the limit of x and y coordinates
            return(list(x=rg_x, y=rg_y))
        } else {
            # range of the values in X for x and y axis limits. This will let the range of all the plots be same
            rg <- range(X)
            for(j in pairs) {
                plot_a_pair_vectors(X[,j[2]], X[,j[1]], groupID, TRUE, 
                                                 (if(log) "xy" else ""), 
                                                 main= "", xlab= lab[j[2]], ylab= lab[j[1]], 
                                                 xlim= rg, ylim= rg, legend_par= legend_par, ...)
            }
            mtext(text= main, outer= TRUE, padj= -1, cex= par()$cex.main, font= 2)
            # return the limit of x and y coordinates
            return(invisible(rg))
        }
    }
)


# -----------------------------------------------------------------------------#
#' @include by.meanBy.R
#' @export
#' @rdname plot_compareMeanOfEach_inPairs
# -----------------------------------------------------------------------------#
# created  : 2011-07-14 by Mun-Gwan
# modified : 2013-07-03 by Mun-Gwan : adapt to "SBAe" class
# -----------------------------------------------------------------------------#

setMethod(
    "plot_compareMeanOfEach_inPairs", signature(X = "BAf"), 
    function(X,
             geometric = FALSE,
             by = "plate",
             log = TRUE,
             MAplot = FALSE,
             main = paste0("Comparing", if (geometric) " geometric" else "", " means"),
             legend_par = list(x = "bottomright"),
             prefix_plate = by,
             ...) {
        M <- meanBy(X, by, geometric)
        if(!is.null(prefix_plate) || prefix_plate == "") 
            colnames(M) <- paste(prefix_plate, colnames(M))
        groupID <- b_batch(X)
        plot_compareMeanOfEach_inPairs(M, groupID, log, MAplot, main, legend_par, ...)
    }
)

