
#' @importFrom utils packageVersion
.onLoad <- function(libname, pkgname) {
  if (packageVersion("SAMtool") != "1.6.1.999") {
    warning("SPH MSE was designed for SAMtool version 1.6.1.999 (custom TMB model for the assessment model). Version ",
            packageVersion("SAMtool"), " is detected.")
  }
}
