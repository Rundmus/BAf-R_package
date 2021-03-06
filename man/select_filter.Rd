% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3.extract.R
\name{select,BAf-method}
\alias{select,BAf-method}
\alias{filter,BAf-method}
\title{Extract like dplyr::select, filter}
\usage{
\S4method{select}{BAf}(.data, ...)

\S4method{filter}{BAf}(.data, ..., .preserve = FALSE)
}
\arguments{
\item{.data}{an object of the \code{\link{BAf-class}}}

\item{...}{Same as \code{...} in \code{\link{select}} or \code{\link{filter}}
in \code{dplyr} package. This is evaluated with \code{sB(.data)} or
\code{sI(.data)}, when \code{select} or \code{filter} is invoked,
respectively.}
}
\value{
an object of the \code{\link{BAf-class}} that contains a part
}
\description{
These functions works like \code{select} and \code{filter} of \code{dplyr}
package.
}
\examples{
data(sba)
# the data of "T2" and "T3"
sba2 <- sba \%>\% select(id \%in\% c("T2", "T3"))
summary(sba2)

# the data of the samples on plate 2
sba3 <- sba \%>\% select(sba == "BA0")
summary(sba3)

data(sba)
# the data of the samples on plate 2
sba4 <- sba \%>\% filter(plate == "2")
summary(sba4)

# the data of cancer patients only
sba5 <- sba \%>\% filter(dis == "cancer")
summary(sba5)

}
\seealso{
\code{\link{select}}

\code{\link{filter}}
}
\author{
Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
}
