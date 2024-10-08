% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/proj-figures.R
\name{plot_Kobe_time}
\alias{plot_Kobe_time}
\title{Kobe time figure}
\usage{
plot_Kobe_time(MSE, output, type = c("B", "F"))
}
\arguments{
\item{MSE}{\link[MSEtool:MSE-class]{MSEtool::MSE} object}

\item{output}{Data frame for plotting. Used for averaging across operating models. (Optional)}

\item{type}{Character for biomass or fishing mortality figure}
}
\value{
A ggplot2 object
}
\description{
Plots proportion of simulations in each Kobe region by year and MP for a single operating model
}
