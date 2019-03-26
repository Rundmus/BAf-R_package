# -----------------------------------------------------------------------------#
#' Bead count plot
#' 
#' Display bead counts across sample and binders.
#' 
#' @inheritParams plot_QC_binder_signal_boxplot
#' @param bead_count 1) a \code{data.frame} of bead counts having sample IDs at
#'   the first column like the file from LIMS, 2) a file name that contains such
#'   data, or 3) a \code{matrix} of such values only where sample IDs are in
#'   \code{rownames}
#' @param sample_order the order of sample in plots. The default is same order 
#'   as given table by \var{bead_count}.
#' @param color_tbl_s colors of samples on the plot. The information should be 
#'   provided as a \code{\link{data.frame}} that has \code{name}, \code{reg} 
#'   (regular expression to find the sample), \code{col} (color), \code{type} 
#'   (\code{negative} or \code{replicated}).\cr If this is NULL, then no 
#'   coloring is applied. The default is
#'   \code{\link{default_color_SBA_ctrl_samples}(sampleid)}, where
#'   \code{sampleid} is a vector of sample IDs that are extracted from
#'   \var{bead_count}
#' @param xlab two characters. The first one is for the plot of binders, while
#'   the 2nd is for sample plot
#' @param main,ylab,ylim,... same as those for \code{\link{boxplot}}
#' 
#' @return a \code{list} of \code{color_tbl_b} and \code{color_tbl_s}
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' 
#' @include plot_QC_0_underlying.R
#' @export
# -----------------------------------------------------------------------------#
# created  : 2014-11-25 by Mun-Gwan
# modified : 
#   2016-09-16 by Mun-Gwan : allow matrix in 'bead_count' allow a single
#     character in 'assayid'
#   2017-08-20 by Mun-Gwan : to adapt to BAf-class
# -----------------------------------------------------------------------------#

plot_QC_bead_count_boxplot <- function(
  bead_count,
  sample_order = NULL,
  color_tbl_b = default_color_SBA_bead(T),
  color_tbl_s = default_color_SBA_ctrl_samples(sampleid), 
  main = deparse(substitute(bead_count)), 
  xlab = c("Beads", "Sample (sorted)"),
  ylab = "Bead count",
  ylim,
  ...
) {
  
  stopifnot(!missing(bead_count))
  force(main)       # avoid lazy evaluation after change of 'bead_count'
  
  ##  Decompose 'bead_count' into 'sampleid' and tgBC
  if(is.data.frame(bead_count)) {
    sampleid <- bead_count[,1]
    tgBC <- as.matrix(bead_count[, -1]) 
  } else if(is.matrix(bead_count)) {
    if(is.null(rownames(bead_count))) {
      sampleid <- paste0("S", 1:nrow(bead_count))   # default sid
    } else {
      sampleid <- rownames(bead_count)
    }
    tgBC <- bead_count
  } else if(is.character(bead_count)) {
    tgBC <- read.delim(gzfile(bead_count), stringsAsFactors = FALSE, 
                       check.names = FALSE, comment.char= "#")
    sampleid <- tgBC[,1]       # first column is assumed to have sid
    tgBC <- as.matrix(tgBC[, -1])
  } else stop("'bead_count' was not given properly.")
  sampleid <- rm_repeat_indicator(sampleid)
  
  if(is.null(sample_order)) sample_order <- seq_len(nrow(tgBC))
  
  if(nrow(tgBC) != length(sample_order))
    stop("The length of 'sample_order' should be same as the number",
         " of rows of the bead count table.")
  
  ##  reorder
  sampleid <- sampleid[sample_order]
  tgBC <- tgBC[sample_order, ]
  
  # color_tbls
  color_tbl_b <- .chq_color_tbl(color_tbl_b)
  color_tbl_s <- .chq_color_tbl(color_tbl_s)
  
  ##   Plot : Bead count distribution   ##
  opar <- par(mfcol= c(2, 1), oma= c(0, 1, 4, 0), 
              mar= par()$mar + c(-2, 1, -2, 0),
              ask= dev.interactive())
  on.exit(par(opar))
  
  # upper limit = 500
  tgBC[tgBC > 500] <- 500
  
  if(missing(ylim)) ylim <- c(0, max(tgBC, na.rm= T))   
  
  boxplot_beadcount <- function(tg, ids, ylim, color_tbl, xlab, ylab, ...) {
    # plot frame for boxplot
    plot(c(1,length(ids)), ylim,  
         type= "n", axes= F, xlab= "", ylab= "", main= "")
    abline(h= BEADCOUNT_TOO_LOW, col= "skyblue2")
    
    .boxplot_color_tbl(tg, ids= ids, 
                       color_tbl= color_tbl, grep= T,
                       main= "", xlab= xlab, ylab= ylab, 
                       ..., add= T) %>%
      .$used_col
  }
  
  # *  plot for binder-wise boxplot  * #
  color_tbl_b_used <- 
    boxplot_beadcount(tgBC, colnames(tgBC), ylim, color_tbl_b, 
                      xlab[1], ylab, ...)
  
  ##   Pointers for the beads with too low bead count
  is_too_low_bc_binder <- apply(tgBC, 2, median, na.rm= T) < BEADCOUNT_TOO_LOW
  if(any(is_too_low_bc_binder)) {
    points(c(1:ncol(tgBC))[is_too_low_bc_binder], 
           rep(par()$usr[3L], sum(is_too_low_bc_binder)), 
           pch= 17, col= "deeppink")
  }
  
  # *  plot for sample-wise boxplot  * #
  color_tbl_s_used <- 
    boxplot_beadcount(t(tgBC), sampleid, ylim, color_tbl_s, 
                      xlab[2], ylab, ...)
  
  title(main= main, outer= T)
  
  return(invisible(list("color_tbl_b"= color_tbl_b_used, 
                        "color_tbl_s"= color_tbl_s_used)))
}





# -----------------------------------------------------------------------------#
#' plot for basic QC of SBA data
#' 
#' Generate plots to check basic QC of SBA data. This invokes other
#' \code{plot_QC} functions, such as 
#' \code{\link{plot_QC_sample_signal_boxplot}},
#' \code{\link{plot_QC_binder_signal_boxplot}}, \code{\link{plot_QC_repl_var}},
#' and \code{\link{plot_QC_bead_count_boxplot}}
#' 
#' @inheritParams plot_QC_sample_signal_boxplot
#' @inheritParams plot_QC_binder_signal_boxplot
#' @param baf an object of the \code{\link{BAf-class}} that contains SBA data
#' @param prefix a prefix of generating files of basic QC plots. Often, it is
#'   sample layout ID.
#' @param path path where the plots will be saved
#' @param bead_count 1) a \code{data.frame}, 2) a file name or 3) a
#'   \code{matrix}. check \code{\link{plot_QC_bead_count_boxplot}}
#' @param main same as \code{\link{boxplot}}
#' @param ylab_signal the label of y-axis of the plots that show signal
#'   distribution in boxplots.
#' 
#' @return \code{ask_save.plot_QC.SBA} 
#' 
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' 
#' @include plot_QC_0_underlying.R
#' @export
# -----------------------------------------------------------------------------#
# created  : 2014-11-24 by Mun-Gwan
# modified : 
#   2015-05-18 by Mun-Gwan : fix plot_QC_repl_var when there is only one
#     replicate per a plate
#   2017-08-20 by Mun-Gwan : to adapt to BAf-class
#   2018-04-10 by Mun-Gwan : fix
# -----------------------------------------------------------------------------#

ask_save.plot_QC.SBA <- function(
  baf,
  prefix = deparse(substitute(baf)),
  path = "../output/QC",
  bead_count = NULL,
  sample_order = order(batch(baf, "sample")),
  color_tbl_b = default_color_SBA_bead(T),
  color_tbl_s = default_color_SBA_ctrl_samples(sid(baf)),
  main = paste(unique(batch(baf, "binder")), collapse= ", "),
  ylab_signal = "Raw MFI",
  guide_line = TRUE
) {
  
  stopifnot(inherits(baf, "BAf"))
  baf <- .reorder_BAf_and_replace_0(baf, sample_order, T)
  
  if(missing(prefix) && nchar(prefix) > 15) prefix <- "SBA"
  
  answer <- readline(
    paste0("Save plots of ", prefix, " on ", path, "\n", 
           "(y : save, s : show on screen)?")) %>%
    tolower()
  
  if(! file.exists(path) && answer == "y") {
    warning("The path ", path, 
            " doesn't exist. So, it has been changed to current folder.")
    path <- "."
  }
  
  save_run_close <- function(ans, fn, path, prefix, FUN, ...) {
    if(ans == "y") {
      filename <- file.path(path, paste0(prefix, fn, ".pdf"))
      pdf(filename, width= 14, ...)
      out_r <- FUN
      dev.off() 
      cat(filename, "   [ saved ! ]\n")
      return(out_r)
    } else {
      return(FUN)
    }
  }
  
  out <- list("answer"= answer)
  if(answer %in% c("y", "s")) {
    prefix <- paste(prefix, "- ")
    
    ##  Plot : Signal distribution across samples  ##
    out <- save_run_close(
      answer, "signal distribution across samples", path, prefix, 
      FUN= plot_QC_sample_signal_boxplot(
        baf= baf, 
        sample_order= 1:nrow(baf),        # already sorted
        color_tbl_s= color_tbl_s, main= main, ylab= ylab_signal
      )
    ) %>%
    { c(out, .) }
    
    ##  Plot : Signal distribution across binders  ##
    is_neg_s <- color_tbl_s %>% 
      dplyr::filter(type == "negative") %>% 
      transmute(i = lapply(reg, grep, .fit_to_std_format(sid(baf)))) %>% 
      unlist() %>% 
      `[<-`(logical(nrow(baf)), ., T)       # unwhich
    
    out <- save_run_close(
      answer, "signal distribution across antibodies excl.ctrlSamples", 
      path, prefix,
      FUN= plot_QC_binder_signal_boxplot(
        baf= baf[ !is_neg_s, ], 
        color_tbl_b= color_tbl_b, main= main, ylab= ylab_signal, 
        guide_line= guide_line
      )
    ) %>%
    { c(out, .) }
    
    ##   Plot : CV   ##
    repl <- out$color_tbl_s %>%
      dplyr::filter(type == "replicated") %>% 
      dplyr::select(name, reg)
    if(!is.null(repl) && nrow(repl) > 0) {	
      for(ii in 1:nrow(repl)) {
        save_run_close(
          answer, paste("CV distribution of", repl$name[ii]), path, prefix,
          FUN= plot_QC_repl_var(
            baf= baf, 
            # fix on 2018-04-10
            i_repl= grepl(repl$reg[ii], .fit_to_std_format(sid(baf))), 
            main= main)
        )
      }
    }
    
    ##   Plot : Bead count distribution   ##
    if(!is.null(bead_count)) {
      save_run_close(
        answer, "bead counts of antibodies", path, prefix, height= 14,
        FUN= plot_QC_bead_count_boxplot(
          bead_count, 
          sample_order= sample_order, 
          color_tbl_b= out$color_tbl_b, 
          color_tbl_s= out$color_tbl_s, 
          main= main
        )
      )
    }
  }
  return(invisible(out))
}
