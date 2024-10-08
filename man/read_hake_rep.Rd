% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/om.R
\name{read_hake_rep}
\alias{read_hake_rep}
\alias{read_hake_dat}
\title{Read hake ADMB assessment}
\usage{
read_hake_rep(
  filename = "admb-assessment/modelo_140.rep",
  headings = "^[[:alpha:]]"
)

read_hake_dat(filename = "admb-assessment/datos_14021.dat", headings = "^#")
}
\arguments{
\item{filename}{Directory for location and name of file}

\item{headings}{Regex for identifying lines in the text file that contain variable names}
}
\value{
Named list
}
\description{
Reads in text files for the ADMB report \code{read_hake_rep()} and data \code{read_hake_dat()} files for the hake assessment. Code is adapted from \code{R2admb::read_rep()}.
}
