% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_processing.R
\name{qc_report}
\alias{qc_report}
\title{qc_report}
\usage{
qc_report(x, column, samples)
}
\arguments{
\item{x}{Seurat object}

\item{column}{a metadata column with values of "Pass" or "Fail"}

\item{samples}{a metadata column of sample names}
}
\description{
Report cells remaining from each sample after quality control
}
