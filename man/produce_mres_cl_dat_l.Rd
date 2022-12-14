% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hub_2_mres_tbl.R
\name{produce_mres_cl_dat_l}
\alias{produce_mres_cl_dat_l}
\title{Produce cluster specific HiC data-object}
\usage{
produce_mres_cl_dat_l(
  tmp_cl,
  chr_spec_res,
  tmp_cl_res,
  tmp_cl_res_set,
  res_num,
  chr_dat_l
)
}
\arguments{
\item{tmp_cl}{target cluster to visualise}

\item{chr_spec_res}{BHiCect result object for target cluster chromosome}

\item{tmp_cl_res}{Resolution at which target cluster was detected}

\item{tmp_cl_res_set}{Set of resolutions contained within target cluster}

\item{res_num}{Complete set of resolutions present in HiC data}

\item{chr_dat_l}{List of chromosome-wide HiC data where each element correspond to a particular resolution}
}
\value{
List where each element correspond to cluster HiC data at a particular resolution
}
\description{
Produce cluster specific HiC data-object
}
