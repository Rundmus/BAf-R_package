% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_QC_0_underlying.R
\name{.boxplot_color_tbl}
\alias{.boxplot_color_tbl}
\title{boxplot using color table}
\usage{
.boxplot_color_tbl(x, ids, color_tbl, grep = FALSE, main, xlab, ylab, ...)
}
\arguments{
\item{x}{a matrix}

\item{ids}{a vector of IDs}

\item{color_tbl}{a \code{tbl_df} having coloring variables. All
\code{'reg'}, \code{col}, and \code{name} are required.}

\item{grep}{same as \code{\link{farg}}}

\item{main, xlab, ylab, ...}{same as \code{\link{boxplot}}}
}
\value{
a list of \code{fg_col} and \code{used_col} (a subset of
  \var{color_tbl} of the colors used here)
}
\description{
Find colors from a given color table and draw a box-and-whisker plot
}
