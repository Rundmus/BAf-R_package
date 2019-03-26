# -----------------------------------------------------------------------------#
#' Text files to share
#' 
#' Write text files to share the data stored in a \code{\link{BAf-class}} object
#' with collaborators
#' 
#' @param baf an object of the \code{\link{BAf-class}}, of which the data will
#'   be converted to text files
#' @param prefix prefix used in file name. An ID for sample set is recommended.
#' @param path where those files should be saved.
#' @param fn_data the characters added into the name of the BAf MFI data file to
#'   indicate it contains data (e.g '-data')
#' @param fn_sample,fn_binder the characters added into the name of the file of
#'   sample/binder information. If this is NULL as the default, the file will be
#'   not generated.
#' @param incl_sample_i,incl_binder_i the columns in BAf@sinfo/@binder that will
#'   be included in 'data'/'binder' file. The column for sample IDs at least
#'   should be included in \code{incl_sample_i}.
#' @param suffix the suffix in the file name 
#' @param fix_sid a named vector of characters. The sample IDs that match to any
#'   of those characters will be replaced with the corresponding name in this
#'   vector. An examples is \code{c(NEG_CTRL = "EMPTY-0001", REPEATED_SAMPLE =
#'   "MIX_1-0001")}
#' @param quote same as \code{quote} in \code{\link{write.table}}. Default is 
#'   TRUE
#' @param sep seperator in the output text files. Tab is the default.
#' 
#' 
#' @author Mun-Gwan Hong, \email{mun-gwan.hong@scilifelab.se}
#' @seealso \code{\link{ask_write.table}}
#' @export
# -----------------------------------------------------------------------------#
# created  : 2017-01-31 by Mun-Gwan
# modified : 
#   2017-08-23 by Mun-Gwan : to adapt to BAf-class
#   2017-10-26 by Mun-Gwan : 
#     default change of 'fn_sample' and 'fn_binder' from missing to NULL
#   2018-03-06 by Mun-Gwan :
#     rm repeat indicator of sid
#   2018-05-08 by Mun-Gwan : change default path to "_sent"
# -----------------------------------------------------------------------------#

write_textfiles_to_share <- 
  function(baf,
           prefix = "BAf",
           suffix = paste0("(", Sys.Date(), ")"),
           path = "_sent",
           fn_data = "-data",
           fn_sample = NULL,
           fn_binder = NULL,
           incl_sample_i = c("id"),
           incl_binder_i = c("id"),
           fix_sid = c(),
           quote = TRUE,
           sep = "\t") {
    
    stopifnot(inherits(baf, "BAf"))
    
    ##  make the sample IDs of "EMPTY" and "MIX_1" clearer
    for(ii in seq_along(fix_sid)) {
      mc <- sid(baf) == fix_sid[ii]
      if(any(mc)) {
        cat("The sample IDs with '", fix_sid[ii],  
            "' were replaced with '", names(fix_sid)[ii], "'.", sep= "")
      }
      sid(baf)[mc] <- names(fix_sid)[ii]
    }
    
    # fn  : 'path', 'prefix', 'suffix', and 'sep' are required
    path_fn <- function(fn) {
      file.path(path, 
                paste0(prefix, fn, suffix, 
                       if(sep == ",") ".csv" else ".txt"))
    }
    
    ##  MFI  -------------------------------------------------------------------
    
    # non-log vs log
    mfi <- round(sX(baf), if(median(sX(baf), na.rm= T) > 100) 2 else 4)
    colnames(mfi) <- rm_repeat_indicator(colnames(mfi))
    
    sI(baf) %>% 
      mutate(id = rm_repeat_indicator(id)) %>% 
      `[`(, incl_sample_i, drop= F) %>%       # sample info
      cbind(., mfi, stringsAsFactors= F) %>% 
      ask_write.table(., file= path_fn(fn_data), 
                      quote= quote, row.names= F, sep= sep)
    
    
    ##  Antibody  --------------------------------------------------------------
    
    if(!is.null(fn_binder)) {
      sB(baf)[, incl_binder_i, drop= F] %>% 
        ask_write.table(., file= path_fn(fn_binder), 
                        quote= quote, row.names= F, sep= sep)
    }
    
    
    ##  Sample Information  ----------------------------------------------------
    
    if(!is.null(fn_sample)) {
      ask_write.table(sI(baf), file= path_fn(fn_sample),
                      quote= quote, row.names= F, sep= sep)
    }
  }
# -----------------------------------------------------------------------------#
#' @describeIn write_textfiles_to_share When 'baf' has SBA data. Just different
#'   default settings.
#' @export
# -----------------------------------------------------------------------------#
# created  : 2017-08-23 by Mun-Gwan
# modified : 
# -----------------------------------------------------------------------------#

write_textfiles_to_share.SBA <- 
  function(baf,
           prefix = "SBA",
           suffix = paste0("(", Sys.Date(), ")"),
           path = "_sent",
           fn_data = "-data",
           fn_sample = NULL,
           fn_binder = NULL,
           incl_sample_i = c("id"),
           incl_binder_i = c("id", "gene_name", "gene_description", "ensg_id"),
           fix_sid = c(NEG_CTRL = "EMPTY-0001", REPEATED_SAMPLE = "MIX_1-0001"),
           quote = TRUE,
           sep = "\t") {
    
    write_textfiles_to_share(baf, prefix, suffix, path, 
                             fn_data, fn_sample, fn_binder, 
                             incl_sample_i, incl_binder_i,
                             fix_sid, quote, sep)
  }

