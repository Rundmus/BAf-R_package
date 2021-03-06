% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trim_outlier.R
\name{trim_outlier}
\alias{trim_outlier}
\alias{trim_outlier,matrix-method}
\alias{trim_outlier,BAf-method}
\title{Trim outlier out by Robust PCA}
\usage{
trim_outlier(X, ...)

\S4method{trim_outlier}{matrix}(
  X,
  ...,
  alpha = 0.9,
  cutoff.od.P = 0.025,
  cutoff.sd.P = 0.025,
  coord = c("o&s", "o", "s", "os"),
  plotit = FALSE,
  main = "",
  kmax = 10
)

\S4method{trim_outlier}{BAf}(
  X,
  ...,
  by_s = NULL,
  alpha = 0.9,
  cutoff.od.P = 0.025,
  cutoff.sd.P = 0.025,
  coord = c("o&s", "o", "s", "os"),
  plotit = FALSE,
  kmax = 10
)
}
\arguments{
\item{X}{a \code{\link{matrix}} or a \code{\link{BAf-class}} object}

\item{...}{in order to encourage to use all parameter names}

\item{alpha, kmax}{the parameter for \code{rrcov::\link{PcaHubert}}}

\item{cutoff.od.P, cutoff.sd.P}{the probability threshold for the orthogonal
and score distances}

\item{coord}{the distance coordinates to be used in outlier classification.
"\code{o&s}" indicates the points beyond the cutoffs on both coordinates
are labelled as outliers. In other options, the points over any cutoff on
any coordinates are marked as outliers.}

\item{plotit}{if plots were wanted}

\item{main}{title of plot}

\item{by_s}{Robust PCA per sample set divided by this. If it is a character,
then the column named as it in \code{@sinfo} wil be used for
stratification. When a factor vector is given, then it is used as it is in
dividing into groups. If it is NULL as the default, there will be no
stratification.}
}
\value{
The BAf object after outlier removal
}
\description{
Trim outlier out based on the orthogonal and score distances computed by 
robust principal components analysis (PCA). After log-transformation, like
ordinary PCA, the values are scaled, but using robust statistics such as
median instead of mean and MAD instead SD. Please note this function is
applied to each binder batch separately and omits the data of the samples
having any NA first.
}
\examples{
data(sba)
B <- trim_outlier(sba, applyBy= "plate", plotit = FALSE)
summary(B)
}
\references{
Hubert, M., Rousseeuw, P. J., Branden, K. V., (2005) ROBPCA: A New Approach 
to Robust Principal Component Analysis. Technometrics 47, 64-79
}
\seealso{
\code{\link{apply_per_group}}
}
\author{
Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
}
\keyword{Array}
\keyword{BAf}
\keyword{Bead}
\keyword{Suspension}
