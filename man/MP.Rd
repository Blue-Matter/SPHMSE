% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mp.R
\name{MP}
\alias{MP}
\alias{Sin_Pesca}
\alias{Manejo_Perfecto}
\alias{PM_A}
\alias{PM_B}
\alias{PM_C}
\alias{I3}
\alias{I5}
\alias{I3_lambda1}
\alias{I5_lambda1}
\title{Hake management procedures}
\usage{
Sin_Pesca(x, Data, reps = 1, ...)

Manejo_Perfecto(x, Data, reps = 1, ...)

PM_A(x = 1, Data, reps = 1, CBA_previous)

PM_B(x = 1, Data, reps = 1, CBA_previous)

PM_C(x = 1, Data, reps = 1, CBA_previous)

I3(
  x = 1,
  Data,
  y = 3,
  lambda = c(2, 1),
  delta = c(0.9, 1.1),
  CBA_previous,
  ...
)

I5(
  x = 1,
  Data,
  y = 5,
  lambda = c(2, 1),
  delta = c(0.9, 1.1),
  CBA_previous,
  ...
)

I3_lambda1(
  x = 1,
  Data,
  y = 3,
  lambda = c(1, 1),
  delta = c(0.9, 1.1),
  CBA_previous,
  ...
)

I5_lambda1(
  x = 1,
  Data,
  y = 5,
  lambda = c(1, 1),
  delta = c(0.9, 1.1),
  CBA_previous,
  ...
)
}
\arguments{
\item{x}{Integer, simulation number}

\item{Data}{\link[MSEtool:Data-class]{MSEtool::Data} object}

\item{reps}{Integer, number of stochastic replicates for generating the advice. Not defined for index-based MPs.}

\item{...}{Additional arguments (not used)}

\item{CBA_previous}{Numeric, the CBA in the previous year, used to calculate the CBA for the following year. Only used if \code{max(Data@Year) == Data@LHYear}.
If not provided, uses a value of 41.4 t (2022 CBA).}

\item{y}{Integer, length of years for calculating the trend, i.e., slope, in the index}

\item{lambda}{Vector length 2, change in CBA relative to slope. First value if slope < 0, second value if lambda > 0}

\item{delta}{Vector length 2, minimum and maximum change in CBA (applied after lambda)}
}
\value{
\link[MSEtool:Rec-class]{MSEtool::Rec} object
}
\description{
Functions of the various management procedures evaluated in the hake MSE (2023)
}
\details{
Sin_Pesca: No fishing reference MP, effort-based and thus will not be subject to CBA implementation error

Manejo_Perfecto: Perfect management is an effort-based MP and thus will not be subject to CBA implementation error
MP specifies the fishing effort that corresponds to 95 percent FMSY
In openMSE, the effort is relative to the last historical year of the operating model.
Effort can also be more efficient/powerful in the projection
This code ensures that the correct fishing mortality is specified regardless of the
improvement in the fishing effort, based on \code{MSEtool::FMSYref()}

PM_A fits an assessment model and applies the low-compliance control rule for the CBA

PM_B fits an assessment model and applies the high-compliance control rule for the CBA

PM_C fits an assessment model and applies the ramped harvest control rule for the CBA

I3 is an index-based management procedure proposed by SSPA (slope over 3 years)

I5 is an index-based management procedure proposed by SSPA (slope over 5 years)

I3_lambda1 is an index-based management procedure proposed by SSPA (slope over 3 years with symmetric lambda)

I5_lambda1 is an index-based management procedure proposed by SSPA (slope over 5 years with symmetric lambda)
}
\section{Model-based MPs}{

Model-based MPs use the RCM model (coded in TMB). By default, the RCM object developed from data to 2023 is located in \code{SPHenv$RCM_hake}.
To update the model for the MP, replace \code{SPHenv$RCM_hake} with a new \code{RCModel-class} object.
}

