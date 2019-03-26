# -----------------------------------------------------------------------------#
#' \code{replace_X} replace 0 with (the lowest measured value - 1)
#' 
#' @inheritParams replace_0
#' @param x the value to be replaced. \code{NA} is also allowed.
#' @param value the values with which the \code{'x'}s will be replaced
#' @param show_count whether the number appearance of the \code{'x'} is displayed
#'   or not
#' 
#' @examples 
#' data(sba)
#' replace_x(sba, NA, 1)
#' 
#' @include apply.R
#' @noRd
# -----------------------------------------------------------------------------#
# created  : 2011-04-17 by Mun-Gwan
# modified : 2017-08-20 by Mun-Gwan : to adapt to BAf-class
# -----------------------------------------------------------------------------#
replace_x <- function(baf,
                      x,
                      value,
                      show_count = TRUE,
                      by_s = NULL,
                      by_b = NULL) {
  stopifnot(inherits(baf, "BAf"),
            !missing(value))
  
  # required : baf
  display_count <- function(x, cond, na.rm, by_s, by_b) {
    cat("The number of ", x, " in @.Data\n")
    cond(sX(baf)) %>% sum(na.rm= na.rm) %>% print()
    if(!missing(by_s) || !missing(by_b)) {
      apply_per_group(baf, FUN= . %>% { sum(cond(.), na.rm= na.rm) }, 
                      by_s= by_s, by_b= by_b) %>% 
        print()
    }
  }
  
  if(is.na(x)) {
    if(show_count) display_count(x, is.na, F, by_s, by_b)
    ##  replacing NA with the given value
    sX(baf)[is.na(sX(baf))] <- value
    return(invisible(baf))
  }
  
  if(show_count) display_count(x, . %>% {. == x}, T, by_s, by_b)
  ##  replacing 0 with the given value
  sX(baf)[sX(baf) == x] <- value
  return(invisible(baf))
}



# -----------------------------------------------------------------------------#
#' replace 0 with given value
#' 
#' replace 0 with a given value, The default is (the lowest value - 1)
#' 
#' @param baf an object of the \code{\link{BAf-class}}
#' @param value the values with which zeros will be replaced
#' @param show_count whether the number appearance of 0s is displayed or not
#' @param by_s,by_b if show_count is TRUE, then these specify by which column the 
#'   counts are stratified, e.g. "plate", "assay". 
#' 
#' @examples 
#' data(sba)
#' replace_0(sba)
#' 
#' @include apply.R
#' @export
# -----------------------------------------------------------------------------#
# created  : 2011-04-17 by Mun-Gwan
# modified : 2017-08-20 by Mun-Gwan : to adapt to BAf-class
# -----------------------------------------------------------------------------#
replace_0 <- function(baf,
                      value = min(sX(baf)[sX(baf) > 1], na.rm = TRUE) - 1,
                      show_count = TRUE,
                      by_s = NULL,
                      by_b = NULL) {
  
  replace_x(baf, 0, value, show_count, by_s, by_b)
}

