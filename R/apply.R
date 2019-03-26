# -----------------------------------------------------------------------------#
#' Apply functions
#' 
#' \code{apply_per_group} : Apply a function to each group which is a subset of 
#' samples and binders stratified by the variables in \code{@sinfo} and 
#' \code{@binder}. The variable names are given by the argument \code{by_s} and
#' \code{by_b}.
#' 
#' @param baf an object of \code{\link{BAf-class}}
#' @param by_s,by_b by which the samples (\var{by_s}) or binders (\var{by_b})
#'   are grouped. The \var{FUN} is applied to each group and the results are
#'   combined. This can be given by 1) a \code{\link{character}} string of a
#'   column name in \code{@sinfo} or \code{@binder}, and 2) a
#'   \code{\link{factor}} vector that can be used in the stratification. If it
#'   is NULL as the default, there will be no stratification. So, the \var{FUN}
#'   will be applied to entire data.
#' @param FUN the function
#' @param ... the arguments that will be passed to FUN
#' @param passBAf if BAf object instead of matrix in @.Data is passed when
#'   applying the \var{FUN}.
#'   
#' @return an object of \code{\link{BAf-class}}
#' 
#' @examples 
#' data(sba)
#' summary(sba)
#' 
#' # mean per plate and assay
#' apply_per_group(sba, mean, by_s= "plate", by_b= "assay")
#' 
#' # calculate median after stratifying samples by "plate" and "sex"
#' apply_per_group(sba, median, na.rm= TRUE, by_s= c("plate", "sex"))
#' 
#' # verify the first median computed by above line
#' median(sba@.Data[with(sba@sinfo, plate == "1" & !is.na(sex) & sex == "female"), ])
#' 
#' # compute coefficients of variation (CV) of \code{MIX_1} wells
#' apply_per_group(sba[sba@sinfo$cohort == "MIX_1", 1:5], FUN= function(j) {
#'     apply(j, 2, function(k) sd(k, na.rm= TRUE)/mean(k, na.rm= TRUE))
#' })
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @export
# -----------------------------------------------------------------------------#
# created  : 2011-09-09 by Mun-Gwan
# modified : 
#   2012-04-23 by Mun-Gwan : upgrade 'applyEachSbaPlate' to 'apply_per_group'
#     1. accept the name of column used for stratification
#     2. mimic the arguments of other apply functions. esp. FUN
#   2012-04-25 by Mun-Gwan : 
#     1. add 'passBAf' for the case FUN handle BAf object instead of matrix in
#     @.Data
#     2. if the output is not in same dimension as input, provide list of the
#     outputs instead of BAf object
#   2012-04-26 by Mun-Gwan : 
#     If the output is a vector, then collect all outputs to generate data frame
#   2013-06-27 by Mun-Gwan : adapt to "SBAe" class
#   2017-08-17 by Mun-Gwan : to adapt to BAf-class
#   2017-11-17 by Mun-Gwan 
#     : fix the error the column names were resorted by BAf.[
# -----------------------------------------------------------------------------#

apply_per_group <-
  function(baf,
           FUN,
           ...,
           by_s = NULL,
           by_b = NULL,
           passBAf = FALSE) {
    # a primitive check up of the arguments
    stopifnot(inherits(baf, "BAf"))
    
    if(any(dim(baf) == 0)) return(NULL)
    
    ##  divide by 'by_b' and 'by_s'
    # find the columns belong to each category divide by 'by_b'
    i_by_b <- .index_grouped_by_s_or_by_b(baf, by_b)
    # find the rows belong to each category divide by 'by_s'
    i_by_s <- .index_grouped_by_s_or_by_b(baf, by_s)
    
    #* dim_chgd = if the dimension has been changed by FUN at any time 
    dim_chgd <- FALSE
    
    FUN <- match.fun(FUN)
    out <- lapply(seq_along(i_by_b), function(iB) {     # per group of binders
      ea_by_b <- if(passBAf) {
        baf[, i_by_b[[iB]] ] 
      } else {
        sX(baf)[, i_by_b[[iB]] ]
      }
      
      out_b <- lapply(seq_along(i_by_s), function(iS) {    # per group of samples
        sub_s_b <- ea_by_b[i_by_s[[iS]], , drop= F]
        attr(sub_s_b, "by_names") <- list(
          by_s = names(i_by_s)[iS], 
          by_b = names(i_by_b)[iB]
        )

        # >> execute FUN here << #
        FUN(sub_s_b, ...)
        
      }) %>% { 
        #  If any return was not an BAf object, then as.matrix
        if(passBAf && all(sapply(., inherits, "BAf"))) . 
        else lapply(., as.matrix) 
      } %>%
        do.call("rbind", .)
      
      ## if the dimension of output has changed,
      if(! identical(dim(out_b), dim(ea_by_b))) {
        assign("dim_chgd", TRUE, envir = parent.env(environment()))
        # if the number of sample is not kept, then return as it is
        if(dim(out_b)[1L] != dim(ea_by_b)[1L]) return(out_b)
      }
      
      #  restore to original order of samples unless the dimension has not been
      #  changed
      return(out_b[order(do.call("c", i_by_s)), , drop= FALSE])
    })
    
    ## if the dimension of the output has NOT been changed
    if(!dim_chgd) {
      ##  binder-wise combine and restore the original order of binders
      out <- do.call("cbind", out) %>%
        `[`(, order(do.call("c", i_by_b)) )
      
      ##  Make sure that the order of rows and columns are UNCHANGED
      errorMsg <- c(
        if(!identical(rownames(out), rownames(baf))) 
          "The order of rows have been changed.\t",
        if(!identical(colnames(out), colnames(baf)))
          "The order of columns have been changed."
      )
      if(!is.null(errorMsg)) stop(errorMsg)
      
      if(passBAf) {
        return(out)
      } else {
        sX(baf) <- as.matrix(out)
        return(baf)
      }
    } 
    
    ##  When the dimension of the output has been changed  ---------------------
    
    ### row names of return
    
    # the number of rows in 'out'
    d_1 <- vapply(out, . %>% {dim(.)[1L]}, 1L)
    
    ## combine the names of groups categorized by "by_s"
    all_i_by_s <- sapply(seq_along(i_by_s), function(j) { # per sample group
      d_i_by_s <- dim(i_by_s)
      dn_i_by_s <- dimnames(i_by_s)
      dnn_i_by_s <- names(dn_i_by_s)
      jj <- j - 1L
      c_names <- list()
      for(k in seq_along(dn_i_by_s)) {
        jjj <- jj%%d_i_by_s[k] + 1L
        jj <- jj%/%d_i_by_s[k]
        c_names <- c(c_names, paste0(dnn_i_by_s[k], ":", dn_i_by_s[[k]][jjj]))
      }
      paste(c_names, collapse= "_")
    })
    
    # if (# of rows) == (# of groups by "by_s")
    if(all(d_1 == length(all_i_by_s))) {
      out <- out %>% 
        lapply(. %>% `rownames<-`(., all_i_by_s))
    }
    
    ### column names of return
    
    d_2 <- vapply(out, . %>% {dim(.)[2L]}, 1L)
    # if all elements in "out" are vectors of same size
    if(all(d_2 == 1) && max(d_1) == min(d_1)) {		
      out <- do.call("cbind", out) %>% 
        `colnames<-`(names(out))
    }
    return(out)
  }





# -----------------------------------------------------------------------------#
#' \code{sapply_per_group} It is a variation of \code{apply_per_group} that
#' works like \code{\link{sapply}}
#' 
#' @param b_simplify if the output should be simplified like \code{simplify} for
#'   \code{sapply}
#' @param s_simplify if the output within stratified group by the \code{by_s}
#'   should be simplified.
#' 
#' @rdname apply_per_group
#' @export
# -----------------------------------------------------------------------------#
# created  : 2014-03-25 by Mun-Gwan
# modified : 
#   2017-08-17 by Mun-Gwan : to adapt to BAf-class
# -----------------------------------------------------------------------------#
sapply_per_group <- function(baf, FUN, ..., by_s = NULL, by_b = NULL, 
                             b_simplify = T, s_simplify = b_simplify) {
  
  ##  divide by 'by_b' and 'by_s'
  # find the columns belong to each category divide by 'by_b'
  i_by_b <- .index_grouped_by_s_or_by_b(baf, by_b)
  # find the rows belong to each category divide by 'by_s'
  i_by_s <- .index_grouped_by_s_or_by_b(baf, by_s)
  
  FUN <- match.fun(FUN)
  sapply(i_by_b, function(jj) {
    l_by_s_ea_b <- lapply(i_by_s, function(ii) baf[ii, jj])
    # if only one sample group
    if(length(l_by_s_ea_b) == 1) {
      # always a data frame is returned even with one column
      return(FUN(l_by_s_ea_b[[1]], ...))		
    } else {
      return(sapply(l_by_s_ea_b, FUN, ..., simplify = s_simplify))
    }
  }, simplify= b_simplify)
}

