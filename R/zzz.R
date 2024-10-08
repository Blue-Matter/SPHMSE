
#' @importFrom utils packageVersion
.onLoad <- function(libname, pkgname) {
  if (packageVersion("SAMtool") != "1.6.1.999") {
    warning("SPH MSE was designed for SAMtool version 1.6.1.999 (custom TMB model for the assessment model). Version ",
            packageVersion("SAMtool"), " is detected.")
  }
}

#' @importFrom utils globalVariables
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("value", "MP", "name", "Year", "Simulation", "OM", "med", "lwr2", "upr2", "Sim", "type", "txt", ".", "scenario"))
}
