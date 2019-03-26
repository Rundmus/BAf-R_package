# -----------------------------------------------------------------------------#
#' Standard curve function
#' 
#' Using the data of standard dilution series, find the function that converts
#' obtained signal into concentration
#' 
#' @param std the obtained signals of standard dilution series
#' @param std_conc the known concentration of the standard series
#' @param ... not used
#' 
#' @return a list of three elements. The first, \code{fn} has a function that
#'   compute concentration when signal measure is given. Other two are
#'   \code{LOD} (limit of detection) and \code{LLOQ} (lower limit of
#'   quantification). The last is an \code{\link{drc}} object that has fitting
#'   model.
#' 
#' @examples
#' std <- c(17933, 8088, 2160, 373, 71, 9, 8, 17076, 8064, 2283, 401, 85, 8, 10)
#' std_conc <- rep(c(34589, 8647, 2161, 540, 135, 0, 0), 2)
#' scurv <- std_curve(std, std_conc)
#' scurv$fn(386)
#' 
#' @author Mun-Gwan Hong, \email{mun-gwan.hong@scilifelab.se}
#' @importFrom drc drm LL.5
#' @export
# -----------------------------------------------------------------------------#
# created  : 2018-03-23 by Mun-Gwan
# modified : 
#   2018-05-09 by Mun-Gwan : Use log-transformed values when computing LLOQ and 
#     LOD computation
#   2018-07-25 by Mun-Gwan : Use original values when computing LLOQ and LOD 
#     computation
# -----------------------------------------------------------------------------#

std_curve <- function(
  std,
  std_conc,
  ...
) {
  stopifnot(length(std) == length(std_conc))
  log_std <- log10(std)
  
  fit_model <- drc::drm(log_std ~ std_conc, fct= drc::LL.5())
  
  b <- coef(fit_model)[["b:(Intercept)"]]
  c <- coef(fit_model)[["c:(Intercept)"]]
  d <- coef(fit_model)[["d:(Intercept)"]]
  e <- coef(fit_model)[["e:(Intercept)"]]
  f <- coef(fit_model)[["f:(Intercept)"]]
  
  fn <- function(x) {
    stopifnot(is.numeric(x))
    log_x <- log10(x)
    (( ((d - c)/(log_x - c))^(1/f) ) - 1)^(1 / b) * e
  }
  
  
  ##  Compute LOD and LLOQ
  
  i_min_conc <- which(std_conc == min(std_conc, na.rm= T))
  
  # m_min <- mean(log_std[i_min_conc], na.rm= T) # out @2018-07-25
  # sd_min <- sd(log_std[i_min_conc], na.rm= T)  # out @2018-07-25 
  
  m_min <- mean(std[i_min_conc], na.rm= T)
  sd_min <- sd(std[i_min_conc], na.rm= T)
  
  if(is.na(sd_min)) {
    lod <- NA
    lloq <- NA
  } else {
    lod <- m_min + 3 * sd_min
    lloq <- m_min + 10 * sd_min
  }
  
  # lod <- 10^lod   # out @2018-07-25
  # lloq <- 10^lloq # out @2018-07-25
  
  return(
    list(fn= fn, LOD= lod, LLOQ= lloq, fit_model= fit_model)
  )
}

