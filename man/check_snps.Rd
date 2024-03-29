% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_loading.r
\name{check_snps}
\alias{check_snps}
\title{Check if SNPs are good for use in analyses and mark them as such.}
\usage{
check_snps(dat, analyses = c("mr", "coloc"), drop = T)
}
\arguments{
\item{dat}{A data.frame of formatted data (exposure or outcome)}

\item{analyses}{Which analyses should be checked?}

\item{drop}{Whether to drop SNPs if they failed the check}
}
\value{
Data.frame
}
\description{
Check if SNPs are good for use in analyses and mark them as such.
}
\details{
List of analyses and what data are checked for:
\itemize{
 \item{"MR"}{beta, SE, P value}
 \item{"coloc"}{chromosome, position, P value}
}
}
