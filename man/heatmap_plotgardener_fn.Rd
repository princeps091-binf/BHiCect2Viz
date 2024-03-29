% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gardener_plot.R
\name{heatmap_plotgardener_fn}
\alias{heatmap_plotgardener_fn}
\title{Heatmap plotting with plotgardener}
\usage{
heatmap_plotgardener_fn(
  cl_col_dat,
  col_specs,
  chromo,
  genome,
  figure_height,
  figure_width,
  res_num
)
}
\arguments{
\item{cl_col_dat}{Tibble collecting all heatmap entries to plot with corresponding coordinates}

\item{col_specs}{List collecting the necessary components to build color-map}

\item{chromo}{String indicating the chromosome containing the considered cluster}

\item{genome}{String indicating the genome symbol containing the considered cluster}

\item{figure_height}{Height in inches of the heatmap to plot}

\item{figure_width}{Width in inches of the heatmap to plot}

\item{res_num}{Named vector with the resolutions and their labels present in the HiC data}
}
\value{
Heatmap plot produced by the plotgardener backend
}
\description{
Heatmap plotting with plotgardener
}
