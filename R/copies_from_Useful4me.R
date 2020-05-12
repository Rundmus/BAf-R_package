# -----------------------------------------------------------------------------#
# Old functions from 'Useful4me' package
# -----------------------------------------------------------------------------#



# -----------------------------------------------------------------------------#
#' Box plot and number
#' 
#' Make a box-and-whisker plot with showing the number of elements in each group
#' If Y-axis is in log-scale, outliers are determined after log-transformation, 
#' unlike \code{boxplot}.
#' 
#' @param formula,data,subset,na.action same as those for \code{\link{boxplot}}
#' @param ... arguments to be passed to \code{\link{boxplot}}
#' @param test the test to check the difference between groups. If this is 
#'   \code{NULL} or "", no test will be performed.
#' @param test_txt_par a list of parameters to control on displaying the test 
#'   result. refer to \code{\link{text}}. In here, \code{x} can be 
#'   \code{"center"}, \code{"right"}, or numerical vectors as \code{x} of 
#'   \code{\link{text}}. The default is \code{x = "center", y= par("usr")[4],
#'   labels = 'pval', adj = c(0.5, 1.5), col = "navy", cex = 0.8)}.
#' @param nEle_txt_par a parameter list for displaying the number of elements of
#'   each group. The default is \code{list(x = seq_along(nEle), y =
#'   par("usr")[3], labels = nEle, col = "brown", cex = 0.8, pos = 3)}.
#' 
#' @examples 
#' data <- data.frame("y"= rnorm(20), "x"= gl(2, 10))
#' boxplot_nEle(y ~ x, data, test= "t.test", test_txt_par= list("col"= "blue"))
#' 
#' @author Mun-Gwan Hong, \email{mun-gwan.hong@scilifelab.se}
#' @seealso \code{\link{boxplot}}
#' @include overwrite_par.R
#' @importFrom Useful2me overwrite_par
# -----------------------------------------------------------------------------#
# created  : 2011-09-13 by Mun-Gwan
# modified : 
#   2013-10-04 by Mun-Gwan : add more parameters, "subset", "na.action"
#   2016-11-22 by Mun-Gwan 
# -----------------------------------------------------------------------------#

boxplot_nEle <- function(formula,
                         data,
                         ...,
                         subset,
                         na.action,
                         test = c("kruskal.test", "aov", "t.test"),
                         test_txt_par = list(),
                         nEle_txt_par = list()) {
  
  stopifnot(inherits(formula, "formula"))
  
  ##   Show ‘boxplot' ----------------------------------------------------------
  
  mcl <- match.call()
  tmp <- match(c("nEle_txt_par", "test", "test_txt_par"), names(mcl), 0L)
  # exclude the parameters for nEle and test
  mcl <- if(all(tmp == 0)) mcl else mcl[-tmp]
  # replace 'boxplot_nEle' with 'boxplot'
  mcl[[1L]] <- quote(graphics::boxplot)
  
  ##   when (log has "y")   ####
  ##   If Y-axis is set to show in log-scale,
  ylog <- FALSE
  if(!is.null(mcl$log)) ylog <- grepl("y", mcl$log)        # log has "y"
  
  if(ylog) {
    mcl$formula <- update(formula, log10(.) ~ .)
    mcl$log <- sub("y", "", mcl$log)
    if(!is.null(mcl$ylim)) mcl$ylim <- log10(mcl$ylim)
    if(is.null(mcl$yaxt)) mcl$yaxt <- "n"
    
    # **  make a box-and-whisker plot hiding labels in Y-axis  **
    fn_out <- eval(mcl, parent.frame())
    
    aT <- axisTicks(par("usr")[3:4], log= T)
    # show labels in Y-axis
    axis(2, at= log10(aT), labels= aT)        
  } else {
    # **  make a box-and-whisker plot  **
    fn_out <- eval(mcl, parent.frame())
  }
  
  
  ##   Show the number of elements ---------------------------------------------
  
  ##   find how many elements are in each group from the given 'formula'
  mf <- match.call(model.frame.default, expand.dots = F)
  mf[[1L]] <- quote(stats::model.frame.default)
  mf <- mf[names(mf) != "..."]
  nEle <- table(eval(mf, parent.frame())[, 2])
  
  # the default parameters for 'nEle_txt_par' 
  e_par <- list(x = seq_along(nEle),
                y = par("usr")[3],
                labels = nEle,
                col = "brown",
                cex = 0.8,
                pos = 3
  )
  e_par <-	overwrite_par(e_par, nEle_txt_par)
  if(ylog) e_par$y <- log10(e_par$y)
  do.call("text", e_par)       # <<<  Show 'nEle'
  
  
  ##   Show the test result ----------------------------------------------------
  
  ##    Hypothesis test - equal means across groups
  # if "" or NULL, then skip
  if(test[1] == "" || is.null(test)) return(invisible(fn_out))
  
  test <- match.arg(test)
  #  only take these 4 arguments for test
  tmp <- c("formula", "data", "subset", "na.action")
  hypoT <- do.call(test, as.list(mcl)[match(tmp, names(mcl), 0L)], 
                   envir= parent.frame())
  
  out <- switch(test, 
                kruskal.test= list(pval= hypoT$p.value, test_pre= "KW test"), 
                aov= list(pval= summary(hypoT)[[1]][1, 5], test_pre= "ANOVA"),
                t.test= list(pval= hypoT$p.value, test_pre= "T-test"))
  
  ## show P-value of the test with H0, same mean across groups
  t_par <- list()
  t_par$labels <- if(out$pval < 0.05) {
    pvaltxt <- format(out$pval, digits = 3, scientific = TRUE)  #1.23e-02 format
    pvaltxt <- strsplit(pvaltxt, "e")[[1]]          #split into "1.23" and "-02"
    eval(substitute( 
      expression(paste(x1, " P =", pval %*% 10^{power})), 
      list(x1= out$test_pre, 
           pval= pvaltxt[1], 
           power= as.numeric(pvaltxt[2])
      )
    ))
  } else paste(out$test_pre, "P =", signif(out$pval, 2))
  
  if(is.null(test_txt_par$x)) test_txt_par$x <- "center"   # default is "center"
  t_par <- c(
    t_par, 
    if(is.numeric(test_txt_par$x)) {
      list(x = test_txt_par$x)
    } else {
      switch(test_txt_par$x, 
             "center" = list(x = mean(seq_along(nEle)), adj = c(0.5, 1.5)),
             "right"  = list(x = par("usr")[2L], adj = c(1.05, 1.5)),
             stop("'x' of 'test_txt_par' should be numeric", 
                  " or 'center' or 'right'."))
    },
    list(
      y= par("usr")[4], 
      col = "navy", 
      cex = 0.8
    )
  )
  t_par <- overwrite_par(t_par, test_txt_par, c("x"))
  if(ylog) t_par$y <- log10(t_par$y)
  do.call("text", t_par)
  
  fn_out$test <- hypoT
  fn_out$p.value <- out$pval
  return(invisible(fn_out))
}

