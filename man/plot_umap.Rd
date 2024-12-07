% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_plotting.R
\name{plot_umap}
\alias{plot_umap}
\title{plot_umap}
\usage{
plot_umap(
  x,
  reduction = "umap",
  group.by,
  split.by = NULL,
  cols = NULL,
  count = T,
  pt.size = 0.5,
  legend.ncol = 1,
  shuffle = F
)
}
\arguments{
\item{x}{Seurat object}

\item{reduction}{reduction name, defaults to "umap"}

\item{group.by}{metadata column to group by}

\item{split.by}{metadata column to split by}

\item{cols}{a vector of colors, defaults to NULL}

\item{count}{display no. of cells if TRUE}

\item{pt.size}{geom_point size}

\item{legend.ncol}{no. of columns for group.by keys in legend}

\item{shuffle}{shuffle points}
}
\description{
Plot UMAP from Seurat object
}
