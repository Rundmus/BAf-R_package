# -----------------------------------------------------------------------------#
#' whether every is BAf
#' 
#' Check every element in the given list is an BAf object
#' 
#' @param baf_list list of \code{\link{BAf-class}} objects
#' @return a logical value
# -----------------------------------------------------------------------------#
# created  : 2017-07-21 by Mun-Gwan
# modified :
# -----------------------------------------------------------------------------#
#' @noRd

is_every_element_BAf <- function(baf_list) {
  baf_list %>%
    vapply(inherits, 1L, "BAf") %>% 
    all()
}

# -----------------------------------------------------------------------------#
#' stop if any is not BAf
#' 
#' Stop with an alarm if any element in the given list is not a BAf object.
#' Because this function is for argument check-up, the error message is as
#' shown below.
#' 
#' @param baf_list list of \code{\link{BAf-class}} objects
#' @return same as \code{\link{stop}}
# -----------------------------------------------------------------------------#
# created  : 2017-07-21 by Mun-Gwan
# modified :
# -----------------------------------------------------------------------------#
#' @noRd

stop_if_any_is_not_BAf <- function(baf_list) {
  if(! is_every_element_BAf(baf_list) ) 
    stop("All arguments must be BAf objects!")
}

# -----------------------------------------------------------------------------#
#' rbind skipping duplicates
#' 
#' Bind in row wise and remove any duplicated rows
#' 
#' @param df_l list of \code{data.frame}
#' @return \code{data.frame}
# -----------------------------------------------------------------------------#
# created  : 2017-07-21 by Mun-Gwan
# modified :
# -----------------------------------------------------------------------------#
#' @noRd

rbind_skip_duplicates <- function(df_l) {
  do.call("rbind", df_l) %>%
    filter({
      apply(., 1, paste, collapse= "_") %>%
        duplicated() %>%
        `!`
    })
}



# -----------------------------------------------------------------------------#
#' bind the codebooks stored in multiple BAf objects
#' 
#' @param baf_list list of \code{\link{BAf-class}} objects
# -----------------------------------------------------------------------------#
# created  : 2013-07-04 by Mun-Gwan
# -----------------------------------------------------------------------------#
#' @noRd

codebook_bind <- function(baf_list, book = c("sinfo", "binder")) {
  stop_if_any_is_not_BAf(baf_list)
  book <- match.arg(book)
  
  if(length(baf_list) == 1) return(baf_list[[1]]@codebook)
  
  baf_list %>%
    lapply(codebook, wise= book) %>%
    rbind_skip_duplicates() %>%
    `rownames<-`(NULL)
}




