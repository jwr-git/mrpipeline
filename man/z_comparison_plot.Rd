% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.r
\name{z_comparison_plot}
\alias{z_comparison_plot}
\title{Z score comparison plot}
\usage{
z_comparison_plot(dat1, dat2, z_col = "z", p_col = "pvalues", verbose = TRUE)
}
\arguments{
\item{dat1}{A list of data}

\item{dat2}{A list of data}

\item{z_col}{Column name for Z scores (Optional)}

\item{verbose}{Display verbose information (Optional, boolean)}

\item{force_static}{True for forcing the plot to be returned as a static plot (Optional)}
}
\value{
Plot
}
\description{
If `plotly` is installed, the plot will be returned interactive.
}
