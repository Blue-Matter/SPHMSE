% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pm.R
\name{ZV_CP}
\alias{ZV_CP}
\alias{PM}
\alias{ZV_MP}
\alias{ZV_LP}
\alias{NZR_CP}
\alias{NZR_MP}
\alias{NZR_LP}
\alias{NSP_CP}
\alias{NSP_MP}
\alias{NSP_LP}
\alias{CBAmin_CP}
\alias{CBAmin_MP}
\alias{CBAmin_LP}
\alias{CBAprom_CP}
\alias{CBAprom_MP}
\alias{CBAprom_LP}
\alias{CBAv_CP}
\alias{CBAv_MP}
\alias{CBAv_LP}
\title{Performance metrics}
\usage{
ZV_CP(MSEobj)

ZV_MP(MSEobj)

ZV_LP(MSEobj)

NZR_CP(MSEobj)

NZR_MP(MSEobj)

NZR_LP(MSEobj)

NSP_CP(MSEobj)

NSP_MP(MSEobj)

NSP_LP(MSEobj)

CBAmin_CP(MSEobj)

CBAmin_MP(MSEobj)

CBAmin_LP(MSEobj)

CBAprom_CP(MSEobj)

CBAprom_MP(MSEobj)

CBAprom_LP(MSEobj)

CBAv_CP(MSEobj)

CBAv_MP(MSEobj)

CBAv_LP(MSEobj)
}
\arguments{
\item{MSEobj}{\link[MSEtool:MSE-class]{MSEtool::MSE} object}
}
\value{
\link[MSEtool:PMobj-class]{MSEtool::PMobj} object
}
\description{
Functions to calculate hake performance metrics
}
\details{
ZV_CP - Probability Kobe Green in Years 1 - 5

ZV_MP - Probability Kobe Green in Years 12 - 15

ZV_LP - Probability Kobe Green in Years 24 - 36

NZR_CP - Probability of avoiding Kobe Red in Years 1 - 5

NZR_MP - Probability of avoiding Kobe Red in Years 12 - 15

NZR_LP - Probability of avoiding Kobe Red in Years 24 - 36

NSP_CP - Probability of avoiding overfishing in Years 1 - 5

NSP_MP - Probability of avoiding overfishing in Years 12 - 15

NSP_LP - Probability of avoiding overfishing in Years 24 - 36

CBAmin_CP - Probability CBA > 20 kt in Years 1-5

CBAmin_MP - Probability CBA > 20 kt in Years 12-15

CBAmin_LP - Probability CBA > 20 kt in Years 24-36

CBAprom_CP - Mean CBA in Years 1-5

CBAprom_MP - Mean CBA in Years 12-15

CBAprom_LP - Mean CBA in Years 24-36

CBAv_CP - Average variability in CBA in Years 1-5

CBAv_MP - Average variability in CBA in Years 12-15

CBAv_LP - Average variability in CBA in Years 24-36
}
