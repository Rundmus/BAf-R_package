# -----------------------------------------------------------------------------#
#' Extract a part of data
#' 
#' This function works in almost same way as the generic form of it.\cr
#' Note: When \code{j} is missing, this takes rows instead of columns.
#' 
#' @param x an object of the \code{\link{BAf-class}}
#' @param i row
#' @param j column
#' @param ...,drop currently not used
#' @param droplevel whether unused factor levels after extraction should be
#'   dropped out in \code{@sinfo} and \code{@binder}
#' 
#' @return an object of the \code{\link{BAf-class}} that contains a part
#' 
#' @examples
#' data(sba)
#' # the samples in 1:4 rows and the values of 3rd target
#' sba2 <- sba[1:4, 3]
#' summary(sba2)
#' 
#' # the data of "T2" and "T3"
#' sba3 <- sba[, c("T2", "T3")]
#' summary(sba3)
#' 
#' # the data of the samples on plate 2
#' sba4 <- sba[sba@sinfo$plate == "2", ]
#' summary(sba4)
#' 
#' # the data of cancer patients only
#' sba5 <- sba[!is.na(sba@sinfo$dis) & sba@sinfo$dis == "cancer", ]
#' summary(sba5)
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @seealso \code{\link{Extract}}
#' @export
# -----------------------------------------------------------------------------#
# created  : 2011-07-13 by Mun-Gwan
# modified : 
#   2011-12-19 by Mun-Gwan fix for the warning "the condition has length > 1..."
#   2012-04-04 by Mun-Gwan
#     :1. change the way to fix the factor levels in '@sinfo'
#      2. fix the problem that arrises when the 'i' is logical vector 
#         (length(i) != 1)
#   2013-06-27 by Mun-Gwan : adapt to "SBAe" class
#   2013-11-14 by Mun-Gwan : allow character matrix in @assay
#   2013-12-06 by Mun-Gwan : let the output has refreshed binder IDs in terms 
#     of "*1"
#   2015-05-28 by Mun-Gwan : When i or j has any NA, give a warning and 
#     remove it
#   2017-09-07 by Mun-Gwan : to adapt to BAf-class
# -----------------------------------------------------------------------------#

setMethod("[", signature(x = "BAf"), function(x, i, j, ..., 
                                              droplevel = TRUE, drop) {
  if(missing(j)) j <- if(ncol(x) == 0) 0 else 1:ncol(x) 
  if(missing(i)) i <- if(nrow(x) == 0) 0 else 1:nrow(x)
  
  warn_na <- function(w) {
    if(any(is.na(w))) 
      warning("There are one or more NAs in '", deparse(substitute(w)), "'!")
  }
  warn_na(i)
  warn_na(j)
  
  extract_tbl <- function(tbl, ii, drop_) {
    if(is.character(ii)) ii <- match(ii, tbl$id)
    tbl[ii, ] %>% { 
      if(drop_) droplevels(.) else .        # drop out unused factor levels
    }
  }
  
  ##  @sinfo & @binder
  subSInfo  <- extract_tbl(x@sinfo, i,  droplevel)
  subBinder <- extract_tbl(x@binder, j, droplevel)
  
  ##  @assy_s
  ## Extract the part of @assy_s removing the dropped out batches
  subAssy_s <- if(length(x@assy_s) > 0) {		# unless @assy_s is an empty list
    tmp <- colnames(x@assy_s[[1]]) %in% unique(subBinder[[x@batch_c$binder]])
    lapply(x@assy_s, . %>% .[i, tmp, drop= FALSE])
  } else x@assy_s
  
  ##  @assy_b
  ## Extract the part of @assy_b removing the dropped out batches
  subAssy_b <- if(length(x@assy_b) > 0) {		# unless @assy_b is an empty list
    tmp <- colnames(x@assy_b[[1]]) %in% unique(subSInfo[[x@batch_c$sinfo]])
    lapply(x@assy_b, . %>% .[j, tmp, drop= FALSE])
  } else x@assy_b
  
  ##  @.Data
  sub.Data <- x@.Data[i, j, drop= FALSE]
  
  BAf(sub.Data, 
      sinfo = subSInfo, 
      binder = subBinder, 
      sinfo_batch_i  = x@batch_c$sinfo, 
      binder_batch_i = x@batch_c$binder,
      assay_sinfo  = subAssy_s, 
      assay_binder = subAssy_b, 
      codebook_sinfo  = x@codebook$sinfo, 
      codebook_binder = x@codebook$binder, 
      note = x@note)
})






# -----------------------------------------------------------------------------#
#' Extract like dplyr::select, filter
#' 
#' These functions works like \code{select} and \code{filter} of \code{dplyr}
#' package.
#' 
#' @param .data an object of the \code{\link{BAf-class}}
#' @param ... Same as \code{...} in \code{\link{select}} or \code{\link{filter}}
#'   in \code{dplyr} package. This is evaluated with \code{sB(.data)} or
#'   \code{sI(.data)}, when \code{select} or \code{filter} is invoked,
#'   respectively.
#' 
#' @return an object of the \code{\link{BAf-class}} that contains a part
#' 
#' @examples
#' data(sba)
#' # the data of "T2" and "T3"
#' sba2 <- sba %>% select(id %in% c("T2", "T3"))
#' summary(sba2)
#' 
#' # the data of the samples on plate 2
#' sba3 <- sba %>% select(sba == "BA0")
#' summary(sba3)
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @seealso \code{\link{select}}
#' @export
#' @rdname select_filter
# -----------------------------------------------------------------------------#
# created  : 2018-04-11 by Mun-Gwan
# modified : 
# -----------------------------------------------------------------------------#

setMethod("select", signature(.data = "BAf"), function(.data, ...) {
  binder_new <- sB(.data) %>% 
    dplyr::filter(...)
  
  .data[, binder_new$id]
})


# -----------------------------------------------------------------------------#
#' @examples
#' data(sba)
#' # the data of the samples on plate 2
#' sba4 <- sba %>% filter(plate == "2")
#' summary(sba4)
#' 
#' # the data of cancer patients only
#' sba5 <- sba %>% filter(dis == "cancer")
#' summary(sba5)
#' 
#' @seealso \code{\link{filter}}
#' @export
#' @rdname select_filter
# -----------------------------------------------------------------------------#
# created  : 2018-04-11 by Mun-Gwan
# modified : 
# -----------------------------------------------------------------------------#

setMethod("filter", signature(.data = "BAf"), function(.data, ...) {
  sinfo_new <- sI(.data) %>% 
    dplyr::filter(...)
  
  .data[sinfo_new$id, ]
})






