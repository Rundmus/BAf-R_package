# -----------------------------------------------------------------------------#
# Some functions and variables for coloring control binders and samples
# -----------------------------------------------------------------------------#
# created  : 2014-11-24 by Mun-Gwan
# modified : 
# -----------------------------------------------------------------------------#


## fit to standard format of regular expression

.fit_to_std_format <- . %>% 
  tolower() %>% 
  gsub("[^[:alnum:]]", "", .) 



##  Color of Antibodies  ---------------------------------------------------

##  Default colors of control antibodies
.color_ctrl_Ab_default <- dplyr::tribble(
  ~name,            ~col,      ~type,
  "Bare bead",      "blue",    "negative",
  "Rabbit IgG",     "green4",  "negative",
  "Anti-Albumin",   "red",     "positive",
  "Anti-Human IgG", "purple",  "positive"
) %>%
  dplyr::mutate(reg = .fit_to_std_format(name) %>% paste0("^", ., "$"))

##  Default colors of other antibodies to measure proteins
.color_measure_Ab_default <- dplyr::tribble(
  ~name,     ~col,      ~type,     ~reg,
  "HPA Ab",  "black",   "measure", "^hpa[[:digit:]]+"
)

##  Default color of various antibodies including the controls
.color_Ab_default <- dplyr::bind_rows(
  .color_ctrl_Ab_default,
  .color_measure_Ab_default
)


# -----------------------------------------------------------------------------#
#' Default color of SBA beads
#' 
#' Returns the tibble for coloring binders
#' 
#' @param only_ctrl whether the colors of control antibodies only will be
#'   returned or not
#' 
#' @return a \code{\link{tibble}} that has \code{'name'}, \code{'col'},
#'   \code{'reg'}, and \code{'type'}
#' @export
# -----------------------------------------------------------------------------#

default_color_SBA_bead <- function(only_ctrl) {
  stopifnot(is.logical(only_ctrl) && length(only_ctrl) == 1)
  if(only_ctrl) {
    .color_ctrl_Ab_default
  } else {
    .color_Ab_default
  }
}




##  Color of Samples -------------------------------------------------------


##  Default colors of control samples of SBA
.color_ctrl_samples_SBA_default <- dplyr::tribble(
  ~name,    ~col,   ~type,      ~reg,
  "EMPTY",  "red",  "negative", "^empty[[:digit:]]{4}"
)


# -----------------------------------------------------------------------------#
#' Default colors of control samples
#' 
#' The default colors of positive and negative control samples
#' 
#' @param sampleid a vector of sample IDs
#' @param rpt_id The ID of replicated samples, the data of which is analyzed to 
#'   check the stability of signals. Often pools are used, whose IDs start with 
#'   "MIX_1". If missing, the IDs with the prefix "MIX_1" are classified as
#'   replicated samples
#'   
#' @return a \code{\link{tibble}} that has 'name', 'col', and 'reg'
#' @examples
#' data(sba)
#' default_color_SBA_ctrl_samples(sampleid= sid(sba))
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @export
# -----------------------------------------------------------------------------#

default_color_SBA_ctrl_samples <-
  function(sampleid, 
           rpt_id = sampleid[grep("^MIX_1", sampleid)]) {
    sampleid <- as.character(sampleid)
    
    rpt_id <- rpt_id %>% unique() %>% sort() 
    
    #  chosen default colors for repeatedly measured samples
    def_color_rpt <- rep(c("blue", "green4", "darkgreen", 
                           "cyan", "magenta", "yellow"), 
                         ceiling(length(rpt_id)/6))
    
    .color_ctrl_samples_SBA_default %>% {
      if(length(rpt_id) > 0) {
        dplyr::add_row(
          .,
          name = rpt_id, 
          reg  = .fit_to_std_format(rpt_id) %>% paste0("^", ., "$"), 
          col  = def_color_rpt[seq_along(rpt_id)], 
          type = rep("replicated", length(rpt_id))
        )
      } else .
    }
  }
