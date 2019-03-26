# -----------------------------------------------------------------------------#
#' Probabilistic Quotient Normaliation
#' 
#' Normalize the data applying Probabilistic Quotient method introduced in
#' Dieterle et al.(2006).
#' 
#' @param X a \code{\link{matrix}} or a \code{\link{BAf-class}} object
#' @param preNorm if the integral normalization should be applied prior to the
#'   PQN. Please refer to the Dieterle \emph{et al.}'s paper.
#' @param ... not used 
#' 
#' @return \code{\link{matrix}} (or \code{\link{BAf-class}} object) after the
#'   normalization
#' 
#' @references 
#' Dieterle et al. (2006) Probabilistic quotient normalization.. - Anal.Chem
#' 
#' @seealso
#' \code{\link{apply_per_group}}
#' \code{\link{normn_MA}}
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @rdname pqn
#' @export
# -----------------------------------------------------------------------------#
# created  : 2011-08-17 by Mun-Gwan
# modified : 
#   2011-09-02 by Mun-Gwan : PQN may follows log-transformation
#   2011-09-08 by Mun-Gwan : Separate into a module (a function per step)
#   2012-04-24 by Mun-Gwan : pqn.BAf instead of 'add_pqnD'
#   2012-10-15 by Mun-Gwan : The calculation for PQN denominator was inversed.
#   2017-08-22 by Mun-Gwan : to adapt to BAf-class
# -----------------------------------------------------------------------------#

setGeneric("pqn", function(X, preNorm = FALSE, ...) {
	# X : the matrix - row:sample, col:target
	# preNorm : if global mean normalization proceeds PQN as suggested by the paper
	
	stopifnot(all(is.finite(X[!is.na(X)])))		# Stop when any infinite value included
	
	if(preNorm) {
	  # it assumes that NA appears in all elements of certain rows
		sSum <- rowSums(X)
		msSum <- mean(sSum, na.rm = T)
		X <- X / sSum * msSum
	}
	
	# reference per target
	ref <- apply(X, 2, median, na.rm = T)
	# normalized for each target
	quo <- t( t(X)/ref )
	
	# normalizing factor for each sample
	nf <- apply(quo, 1, median, na.rm = T)
	
	# any sample of which the median value is 0
	nf[nf == 0] <- NA_real_
	
	return(X / nf)
})

# -----------------------------------------------------------------------------#
#' @param by_s,by_b same as \code{by_s} and \code{by_s} in
#'   \code{\link{apply_per_group}}. If any of these are given the \code{pqn} is
#'   applied per group divided by these variables.
#' 
#' @examples
#' data(sba)
#' sba2 <- pqn(sba[sba@sinfo$cohort != "EMPTY", ])
#' plot_QC_sample_signal_boxplot(sba2)
#' 
#' @rdname pqn
#' @export
# -----------------------------------------------------------------------------#

setMethod(
  "pqn", signature(X = "BAf"), 
  function(X,
           preNorm = FALSE,
           by_s = NULL,
           by_b = NULL,
           ...) {
    
    X0 <- X
    
    X1 <- apply_per_group(X0, pqn, preNorm, by_s= by_s, by_b= by_b)

    ## compute PQN denominator and store it in @assy_s
    bat <- batch(X1, "binder")
    X1@assy_s[["pqn_denominator"]] <- 
      as.data.frame((sX(X0) / sX(X1))[, !duplicated(bat)]) %>% 
      `colnames<-`(bat[!duplicated(bat)])

    return(X1)
  })

