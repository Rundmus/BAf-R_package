# -----------------------------------------------------------------------------#
#' Replicates variation plot
#' 
#' Draw the plot that shows the distribution of variation statistics (e.g.
#' coefficient of variation(CV)) of replicates
#' 
#' @param baf an object of the \code{\link{BAf-class}}
#' @param i_repl an integer or logical index vector of the rows that have
#'   replicated samples
#' @param by_s by which the samples are grouped. The \var{varFUN} is applied to 
#'   each group and the results are combined. This should be given by a 
#'   \code{\link{character}} string of a column name in \code{@sinfo}.
#' @param varFUN a function to compute a variation statistics. The default is
#'   the function to compute the coefficienct of variation.
#' @param col_lines a named vector of line colors in character strings
#' @param main,xlab,ylab,xlim,ylim same as those for \code{\link{plot}}
#' @param lty,lwd the first value is for the line to show the variation of
#'   replicates, while the second is for the other samples. Individual values
#'   are same those for \code{\link{plot}}
#' @param ... the arguments passed over to \code{plot} and \code{lines}
#' 
#' @return a list that contains computed varaition statistics of the replicates
#'   and all the other samples
#' 
#' @examples
#' data(sba)
#' sba2 <- sba[sba@sinfo$cohort != "EMPTY", ]
#' plot_QC_repl_var(sba2, sba2@sinfo$cohort == "MIX_1")
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' 
#' @include plot_QC_0_underlying.R
#' @importFrom Useful2me toupper_1st_char
#' @export
# -----------------------------------------------------------------------------#
# created  : 2014-11-24 by Mun-Gwan
# modified : 2017-08-20 by Mun-Gwan : to adapt to BAf-class
# -----------------------------------------------------------------------------#

plot_QC_repl_var <- function(
  baf,
  i_repl = 1:nrow(baf),
  by_s = batch_colname(baf, "sample"),
  varFUN = function(x) { sd(x, na.rm= T) / mean(x, na.rm= T) }, 
  col_lines = palette(),
  main = deparse(substitute(baf)),
  xlab = "Coef. of var.",
  ylab = "Density",
  xlim = c(0, 0.5),
  ylim = NULL,
  lty = c("solid", "dashed"),
  lwd = c(2, 1),
  ...
) {
  stopifnot(inherits(baf, "BAf"),
            is.integer(i_repl) || is.logical(i_repl),
            is.character(by_s),
            is.numeric(varFUN(sX(baf)[, 1])),
            length(varFUN(sX(baf)[, 1])) == 1)
  if(is.integer(i_repl)) stopifnot(length(i_repl) > 0)
  if(is.logical(i_repl)) {
    stopifnot(length(i_repl) == nrow(baf),
              any(i_repl))
  }
  
  bafl <- list()
  bafl$repl   <- baf[i_repl, ]
  bafl$others <- if(is.logical(i_repl)) baf[!i_repl, ] else baf[-i_repl, ]
  
  # extract replicates and replace 0s.
  bafl <- lapply(bafl, replace_0, value= 1, show_count= F)
  
  opar <- par(ask= dev.interactive())
  on.exit(par(opar))
  
  # *  variation statistics computation
  var_batch <- lapply(bafl, apply_per_group, apply, 2, varFUN, by_s= by_s)
  var_all   <- lapply(bafl, apply_per_group, apply, 2, varFUN, by_s= NULL)
  
  nm_by_s <- names(.index_grouped_by_s_or_by_b(baf, by_s))
  n_batch <- length(nm_by_s)
  #  Check the output was for binders in every by_s
  stopifnot(identical(
    rownames(var_batch[[1]]), rep(bid(baf, exact= T), n_batch)
  ))
  
  var_per_batch <- sapply(names(var_batch), function(ii) {
    ea <- var_batch[[ii]]
    if(is.null(ea)) return(NULL)
    if(n_batch == 1) return(data.frame('All'= var_all[[ii]]))
    
    # transform the longtitudinal table to landscape
    lapply(seq_len(n_batch), function(j) {
      ea[ncol(baf) * (j - 1) + seq_len(ncol(baf)), , drop= F]
    }) %>% 
      do.call("cbind", .) %>% 
      `colnames<-`(nm_by_s) %>% 
      as.data.frame() %>% 
      cbind('All'= var_all[[ii]], .)     # include 'All'
  }, simplify= F)
  density_per_batch <- var_per_batch %>% 
    lapply(lapply, . %>% as.vector %>% {
      if(all(is.na(.))) NULL else na.omit(.) %>% density
    })
  

  #  compute 'ylim' to cover all columns
  yliml <- lapply(density_per_batch, . %>% 
                    sapply(., . %>% .$y) %>% 
                    unlist() %>% 
                    # suppress for the case of NULL
                    {suppressWarnings(range(., na.rm= T))} 
                  )
  if(any(!is.finite(yliml$repl))) return(invisible(var_per_batch))
  if(is.null(ylim)) ylim <- yliml$repl
  
  
  #  scale the 'others' density values to fit to the plot of replicates
  yod <- diff(yliml$others)
  yd <- diff(ylim)
  density_per_batch$others <- 
    lapply(density_per_batch$others, function(ea) {
      ea$y <- ea$y * yd / yod * min(yod / yd, 1.2)
      ea
    })
  
  
  #  if 'col_lines' was shorter than the available groups
  col_lines <- length(density_per_batch[[1]]) %>% {
    rep(col_lines, ceiling(. / length(col_lines)))[seq_len(.)]
  }
  
  #  default names of 'col_lines'
  if(is.null(names(col_lines))) 
    names(col_lines) <- names(density_per_batch[[1]])
  
  if(missing(main)) {
    main <- paste(unique(batch(baf, "binder")), collapse= ", ")
  }
  
  plot(density_per_batch$repl$All, col= col_lines[1L], 
       main= main, xlab= xlab, ylab= ylab, xlim= xlim, ylim= ylim, 
       lty= lty[1], lwd= lwd[1], ...)

  lines(density_per_batch$others$All, col= col_lines[1L], 
        lty= lty[2], lwd= lwd[2], ...)

  stopifnot(identical(names(density_per_batch), c("repl", "others")))
  if(n_batch > 1) {
    for(ii in seq_along(density_per_batch)) {
      ea <- density_per_batch[[ii]]
      if(length(ea) < 2L) next
      for(j in 2L:length(ea)) {
        lines(ea[[j]], col= col_lines[j], 
              lty= lty[ii], lwd= lwd[ii], ...)
      }
    } 
    legend("topright", legend= names(col_lines), 
           text.col= col_lines, bty= "n",
           title= if(is.character(by_s) && length(by_s)) 
             toupper_1st_char(by_s) else NULL)
  }
  
  return(invisible(var_per_batch))
}
