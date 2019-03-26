# -----------------------------------------------------------------------------#
#' Trim outlier out by Robust PCA
#' 
#' Trim outlier out based on the orthogonal and score distances computed by 
#' robust principal components analysis (PCA). After log-transformation, like
#' ordinary PCA, the values are scaled, but using robust statistics such as
#' median instead of mean and MAD instead SD. Please note this function is
#' applied to each binder batch separately and omits the data of the samples
#' having any NA first.
#' 
#' @param X a \code{\link{matrix}} or a \code{\link{BAf-class}} object
#' @param ... in order to encourage to use all parameter names
#' @param alpha,kmax the parameter for \code{rrcov::\link{PcaHubert}}
#' @param cutoff.od.P,cutoff.sd.P the probability threshold for the orthogonal
#'   and score distances
#' @param coord the distance coordinates to be used in outlier classification.
#'   "\code{o&s}" indicates the points beyond the cutoffs on both coordinates
#'   are labelled as outliers. In other options, the points over any cutoff on
#'   any coordinates are marked as outliers.
#' @param plotit if plots were wanted
#' @param main title of plot
#' 
#' @return The BAf object after outlier removal
#' 
#' @references 
#' Hubert, M., Rousseeuw, P. J., Branden, K. V., (2005) ROBPCA: A New Approach 
#' to Robust Principal Component Analysis. Technometrics 47, 64-79 
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @seealso \code{\link{apply_per_group}}
#' @examples 
#' data(sba)
#' B <- trim_outlier(sba, applyBy= "plate", plotit = FALSE)
#' summary(B)
#' @keywords BAf Suspension Bead Array
#' @import rrcov
#' @rdname trim_outlier
#' @export
# -----------------------------------------------------------------------------#
# created  : 2012-04-23 by Mun-Gwan
# modified : 
#   2012-10-15 by Mun-Gwan : fix the bug with the data including NAs disapearing
#     row names
#   2012-10-30 by Mun-Gwan : 
#     1) fix the error that appears when 'applyBy' is missing.
#     2) fix the case 'nrow(tg) < 10'
#   2013-05-13 by Mun-Gwan : 
#     1) handling "coord" variable
#     2) minor changes in diagnostic plot
#   2013-05-21 by Mun-Gwan : change 'kmax' for 'PcaHubert' to default 'kmax= 10' 
#     in rrcov package
#   2013-07-03 by Mun-Gwan : adapt to "SBAe" class
#   2013-07-25 by Mun-Gwan : coord = "o&s" limits the outliers in the top right
#     corner
#   2013-12-06 by Mun-Gwan : take 'kmax' as an argument
#   2014-09-17 by Mun-Gwan : select appropriate "plot" more specifically
#   2015-11-04 by Mun-Gwan : When MAD is 0, then skip scaling
#   2017-07-18 by Mun-Gwan : fix the problem of showing wrong sample ID when NA
#     was already included
#   2017-11-09 by Mun-Gwan : add a function to handle matrix
# -----------------------------------------------------------------------------#

setGeneric("trim_outlier", function(X, ...) standardGeneric("trim_outlier"));

# -----------------------------------------------------------------------------#
#' @rdname trim_outlier
#' @export
# -----------------------------------------------------------------------------#
setMethod(
  "trim_outlier", 
  signature(X = "matrix"), 
  function(X, 
           ...,
           alpha = 0.9,
           cutoff.od.P = 0.025,
           cutoff.sd.P = 0.025,
           coord = c("o&s", "o", "s", "os"), 
           plotit = FALSE,
           main = "",
           kmax= 10) {
    
    # Stop when any infinite value included
    stopifnot(all(is.finite(X[!is.na(X)])))
             
    coord <- match.arg(coord)
    
    ##  
    if(plotit) {
      # two rows in plot and ask = TRUE if it continues to the next plate
      opar <- par(
        mfrow = c(2, 1),
        oma = c(0, 0, 3, 0),
        mar = c(4, 4, 2, 2) + 0.1,
        mgp = c(2.7, 1, 0),
        ask = if (dev.interactive(orNone = TRUE)) TRUE else FALSE
      )
      on.exit(par(opar))
    }
    
    # omit NA of failed samples
    tg <- na.omit(X)
    i_omit <- unclass(attr(tg, "na.action"))
    
    if(nrow(tg) < 10) 
      stop("Too small number (<10) of samples were left after na.omit.")
    if(length(i_omit) > (nrow(X)/2)) 
      stop("More than a half of samples have NAs.")
    
    # log-transform + centering and scaling like ordinary PCA but using median
    # and MAD instead.
    tg <- log(tg) %>% 
      scale(center = apply(., 2, median, na.rm = TRUE), 
            scale  = apply(., 2, mad, na.rm = TRUE) %>% 
              # no scaling when MAD == 0
              ( function(x) { x[x == 0] <- 1; x } )
      )
    
    if(any(tmp <- apply(tg, 2, mad, na.rm= T) == 0))
      warning("The columns (", paste(which(tmp), collapse= ","), 
              ") have constant values for all samples")
    
    ## robust PCA
    pca <- rrcov::PcaHubert(tg, kmax = eval(kmax), alpha = alpha)
    pca@cutoff.od <- if(grepl("o", coord)) {
      rrcov:::.crit.od(pca@od, crit= (1- cutoff.od.P)) 
    } else Inf	#2013-05-13
    
    pca@cutoff.sd <- if(grepl("s", coord)) {
      sqrt(qchisq(cutoff.sd.P, pca@k, lower.tail = FALSE)) 
    } else Inf	#2013-05-13
    
    ## find outliers
    isOut.od <- pca@od > pca@cutoff.od
    isOut.sd <- pca@sd > pca@cutoff.sd
    # 2013-05-13 / 2013-07-25
    isOut <- if(grepl("&", coord)) {
      isOut.od & isOut.sd 
    } else { 
      isOut.od | isOut.sd
    }
    
    ##  * Plot it *  -----------------------------------------------------------
    
    if(plotit) {
      ##  Box plot - show distribution of signals of each sample
      boxplot(t(tg),
              col = ifelse(isOut, 2, 0),
              outcol = ifelse(isOut, 2, par()$fg),
              main = "Intensity distribution",
              cex = 0.5,
              cex.main = 1,
              col.main = "gray",
              xlab = "",
              ylab = "Scaled intensity",
              cex.lab = 1,
              xaxt = "n"
      )
      
      ##   X-axis label
      title(xlab= "Sample", mgp= c(2, 0, 0), col.lab= "gray")
      axis(side= 1,
           at = 1:nrow(tg),
           labels = rownames(tg),
           las = 3,
           cex.axis = 0.3,
           lwd.ticks = 0.5,
           mgp = c(3, 0.7, 0)
      )
      # thicker at every 5 ticks
      axis(side= 1, at= seq(5, nrow(tg), 5), labels= FALSE)
      
      abline(h=0, col= "gray", lty= "dotted")
      if(nrow(tg) > 10) # vertical gray lines per 10 samples
        abline(v=seq(10, nrow(tg), 10), col= "gray", lty= "dotted")
      
      # >> diagonistic plot << #
      selectMethod("plot", c("Pca", "missing"))(
        x= pca,
        id.n.sd = sum(isOut.sd) + 1, 
        id.n.od = sum(isOut.od) + 1,
        main = "Outlier map", 
        xlim = c(0, max(pca@sd, if (grepl("s", coord)) pca@cutoff.sd) * 1.1),
        ylim = c(0, max(pca@od, if (grepl("o", coord)) pca@cutoff.od)),
        cex.main = 1, 
        col.main = "gray", 
        cex.lab = 1, 
        off = 0.03, 
        cex = 0.8,
        col = ifelse(isOut, 2, 1), 
        pch = ifelse(isOut, 19, 1)
      )
      
      title(main= main, outer= T)
    }
    
    # squeeze the omitted row into the 'isOut'
    if(!is.null(i_omit)) {
      for(j in i_omit) isOut <- append(isOut, FALSE, after= (j - 1))
      # return to original row names in the rows of NAs (2012-10-15)
      names(isOut)[i_omit] <- names(i_omit)
    }
    
    # >> replace all values with NA of the outlier samples << #
    X[isOut, ] <- NA
    
    attr(X, "is_out") <- isOut
    
    return(X)
  }
)


# -----------------------------------------------------------------------------#
#' @param by_s Robust PCA per sample set divided by this. If it is a character,
#'   then the column named as it in \code{@sinfo} wil be used for
#'   stratification. When a factor vector is given, then it is used as it is in
#'   dividing into groups. If it is NULL as the default, there will be no
#'   stratification. 
#'   
#' @rdname trim_outlier
#' @export
# -----------------------------------------------------------------------------#
# created  : 2012-04-23 by Mun-Gwan
# modified : 
#   2017-11-09 by Mun-Gwan : add a function to handle matrix
# -----------------------------------------------------------------------------#

setMethod(
  "trim_outlier", 
  signature(X = "BAf"), 
  function(X, 
           ...,
           by_s = NULL,
           alpha = 0.9,
           cutoff.od.P = 0.025,
           cutoff.sd.P = 0.025,
           coord = c("o&s", "o", "s", "os"), 
           plotit = FALSE,
           kmax= 10) {
    
    apply_per_group(
      baf= X, 
      by_s = by_s, 
      by_b = batch_colname(X, "binder"),    # per binder batch
      passBAf = T, 
      FUN= function(TG) {
        main <- if(is.null(by_s)) {
          batch(TG, "binder")[1]
        } else if(is.character(by_s) && length(by_s) == 1) {
          paste(batch(TG, "binder")[1], "-", by_s, attr(TG, "by_names")$by_s)
        } else {
          paste(batch(TG, "binder")[1], "-", attr(TG, "by_names")$by_s)
        }
        
        
        isOut <- trim_outlier(
          X= sX(TG), 
          alpha = alpha,
          cutoff.od.P = cutoff.od.P,
          cutoff.sd.P = cutoff.sd.P,
          coord = coord, 
          plotit = plotit,
          main = main,
          kmax = kmax
        ) %>% 
          attr("is_out")
        
        # >> replace all values with NA of the outlier samples << #
        sX(TG)[isOut, ] <- NA
        sA(TG, "sample")$outlier <- data.frame(isOut) %>% 
          `names<-`(batch(TG, "binder")[1]) 
        
        return(TG)
      }
    )
  }
)
  