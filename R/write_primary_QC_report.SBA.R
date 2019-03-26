# -----------------------------------------------------------------------------#
#' Write primary QC report
#' 
#' Write a primary QC report for SBA data using Rmarkdown
#' 
#' @param baf an object of the \code{\link{BAf-class}} that contains SBA data
#' @param bead_count 1) a \code{data.frame}, 2) a file name or 3) a
#'   \code{matrix}. check \code{\link{plot_QC_bead_count_boxplot}}
#' @param file The name of the generating HTML file. The ".html" extension will
#'   be added if it is missing.
#' @param title The title of the generating report, which is displayed at the
#'   top
#' @param author The author of the report
#' @param date The date written as the creation date
#' @param rm_Rmd if the .Rmd (RMarkdown) file should be removed after .pdf was
#'   generated.
#' @param lowbound_mfi the lower bound of MFI. If a median of a sample is below
#'   this threshold, the data of the sample is neglected afterwords in this QC
#'   report.
#' @inheritParams plot_QC_sample_signal_boxplot
#' @inheritParams plot_QC_binder_signal_boxplot
#' 
#' @author Mun-Gwan Hong, \email{mun-gwan.hong@scilifelab.se}
#' @export
#' @importFrom whisker whisker.render
#' @import rmarkdown
#' @import knitr
#' @import Rtsne
#' 
#' @include plot_QC_0_underlying.R
# -----------------------------------------------------------------------------#
# created  : 2017-09-01 by Mun-Gwan
# modified : 
#   2017-10-27 by Mun-Gwan 
#     - Make 'render' quiet
#     - add the @param rm_Rmd
# -----------------------------------------------------------------------------#

write_primary_QC_report.SBA <- function(
  baf, 
  bead_count = NULL,
  sample_order = order(batch(baf, "sample")),
  file = "primary_QC_report.html",
  title= paste("Primary QC report of", deparse(substitute(baf))), 
  author= "Mun-Gwan",
  date= Sys.Date(),
  rm_Rmd= T,
  color_tbl_b = default_color_SBA_bead(T),
  color_tbl_s = default_color_SBA_ctrl_samples(sid(baf)),
  lowbound_mfi = 100
) {
  stopifnot(inherits(baf, "BAf"))
  if(length(unique(batch(baf, "binder"))) != 1) 
    warning("This function is designed to write a report for",
            "one SBA bead batch.")
  force(title)
  baf <- baf[sample_order, ]
  
  #  template_path = path to the template file
  template_path <- system.file("templates", "Primary_QC_report.SBA.Rmd", 
                               package= "BAf", 
                               mustWork= TRUE)
  
  ##  Substitute the entities in the templates with the given values
  data_in <- list(
    title = title,
    author = author,
    date = date
  )
  out_text <- whisker::whisker.render(readLines(template_path), data= data_in)

  ##  allow both with or without .pdf in 'file'
  base_fn <- sub("\\.htm[l]?$", "", basename(file))

  #  **  Add information into an Rmarkdown  **  #
  writeLines(out_text, file.path(dirname(file), paste0(base_fn, ".Rmd")))
  
  force(color_tbl_b)
  force(color_tbl_s)
  
  rmarkdown::render(
    file.path(dirname(file), paste0(base_fn, ".Rmd")),
    output_format= "html_document",
    quiet= T,
    params= list(
      baf = baf, 
      bead_count = bead_count,
      color_tbl_b = color_tbl_b,
      color_tbl_s = color_tbl_s,
      lowbound_mfi = lowbound_mfi
    )
  )
  
  message('# The document \"', paste0(base_fn, ".html"), '\" has been created.')
  
  if(rm_Rmd) {
    file.remove(file.path(dirname(file), paste0(base_fn, ".Rmd")))
  }
  
  invisible(TRUE)
}
