
#' Read hake ADMB assessment
#'
#' @description Reads in text files for the ADMB report `read_hake_rep()` and data `read_hake_dat()` files for the hake assessment. Code is adapted from `R2admb::read_rep()`.
#' @param filename Directory for location and name of file
#' @param headings Regex for identifying lines in the text file that contain variable names
#' @return Named list
#' @export
read_hake_rep <- function(filename = "admb-assessment/modelo_140.rep", headings = "^[[:alpha:]]") {
  rr <- readLines(filename)

  #numLines <- grepl("^[0-9. e+-]*$", rr)
  #commentLines <- grepl("^#", rr)
  #stringLines <- grepl("^[[:alpha:]]*$", rr)
  #stringLines <- grepl("^[[:alpha:]]", rr)

  stringLines <- grepl(headings, rr)

  vars <- rr[stringLines]

  vars_list <- lapply(1:length(vars), function(x) {
    line_start <- which(rr == vars[x])
    if (x == length(vars)) {
      line_end <- length(rr) + 1
    } else {
      line_end <- which(rr == vars[x+1])
    }

    rr_char <- rr[seq(line_start + 1, line_end - 1)]
    rr_num <- lapply(strsplit(rr_char, " "), as.numeric)

    rr_array <- do.call(rbind, rr_num)

    if (nrow(rr_array) == 1) { # Convert to vector
      rr_array <- as.numeric(rr_array)
      if (is.na(rr_array[1])) rr_array <- rr_array[-1]
    } else {
      if (all(is.na(rr_array[, 1]))) rr_array <- rr_array[, -1]
    }
    return(rr_array)
  })

  names(vars_list) <- vars
  return(vars_list)
}

#' @rdname read_hake_rep
#' @export
read_hake_dat <- function(filename = "admb-assessment/datos_14021.dat", headings = "^#") {

  rr <- readLines(filename)

  #numLines <- grepl("^[0-9. e+-]*$", rr)
  commentLines <- grepl(headings, rr)
  #stringLines <- grepl("^[[:alpha:]]*$", rr)

  vars <- rr[commentLines]

  vars_list <- lapply(1:length(vars), function(x) {

    line_start <- which(rr == vars[x])
    if (x == length(vars)) {
      line_end <- length(rr) + 1
    } else {
      line_end <- which(rr == vars[x+1])
    }

    rr_char <- rr[seq(line_start + 1, line_end - 1)]
    rr_num <- lapply(strsplit(rr_char, "[[:space:]]"), as.numeric)

    rr_array <- do.call(rbind, rr_num)

    if (nrow(rr_array) == 1) { # Convert to vector
      rr_array <- as.numeric(rr_array)
      if (is.na(rr_array[1])) rr_array <- rr_array[-1]
    } else {
      if (all(is.na(rr_array[, 1]))) rr_array <- rr_array[, -1]
    }
    return(rr_array)
  })

  names(vars_list) <- sub("# ", "", vars)
  return(vars_list)

}

#' Create hake operating model from ADMB estimation model
#'
#' Reads in ADMB data and report files to generate an openMSE operating model
#'
#' @param dat Filename of ADMB data file
#' @param rep Filename of ADMB report file
#' @param Name Character, name of operating model to insert to `OM@Name`
#' @param nsim Integer, number of simulations in the operating model
#' @param proyears Integer, number of projection years in the operating model
#' @param AC Numeric, lag-1 autocorrelation for the simulated recruitment deviations in the projection
#' @param Perr Numeric, lognormal standard deviation for the simulated recruitment deviations in the projection
#' @param tv_maturity Logical, whether to use time-varying maturity in the operating model (specified in the ADMB data file)
#' @return [MSEtool::OM-class] object
#' @export
make_hake_OM <- function(
    dat = "admb-assessment/datos_14021.dat",
    rep = "admb-assessment/modelo_140.rep",
    Name = "Merluza comun - OM 1",
    nsim = 100,
    proyears = 36,
    AC = 0.65,
    Perr = 0.5,
    tv_maturity = FALSE # Time-varying maturity
) {

  # Read ADMB report ----
  admb_rep <- read_hake_rep(rep)

  # Read ADMB data file ----
  admb_dat <- read_hake_dat(dat)

  # Historical parameters ----
  nyears <- admb_dat$nanos
  years <- admb_dat$indices[, 1] #1940:2022

  # Operating model parameters ----
  age <- admb_dat$edades #2:13
  nage <- length(age)
  nage_om <- max(age) + 1

  spawn_time_frac <- 0.5833 # Fraction of year when spawning occurs (7/12)
  SRrel <- 2 # (1 = Beverton-Holt, 2 = Ricker)

  R0 <- rep(admb_rep$ro, nsim)
  SB0 <- rep(admb_rep$bo, nsim)
  phi0 <- SB0/R0
  h <- rep(admb_rep$Steepness, nsim)

  # Double check stock recruit parameters
  #alfa <- admb_rep$alfa
  #bbeta <- admb_rep$bbeta
  #MSEtool::SRalphaconv(h, phi0, SR = 2)
  #admb_rep$alfa

  #MSEtool::SRbetaconv(h, R0, phi0, SR = 2)
  #admb_rep$bbeta

  # Natural mortality ----
  M <- local({
    maa <- array(NA, c(nsim, nage_om, nyears))
    for(x in 1:nsim) {
      for(a in 3:nage_om) maa[x, a, ] <- admb_rep[["M"]]
    }
    maa[, 1:2, ] <- 1e-15 + .Machine$double.eps
    maa
  })

  # Fishing mortality ----
  FM <- local({
    faa <- array(NA, c(nsim, nage_om, nyears))
    for(x in 1:nsim) {
      faa[x, -c(1:2), ] <- t(admb_rep[["F"]] * admb_rep[["Sel.flota"]])
    }
    faa[, 1:2, ] <- 0
    faa
  })
  #matplot(years, t(FM[1, , ]), typ = "l", ylab = "F", col = 1:6, lty = 1:5)
  #legend("topleft", legend = 0:max(age), col = 1:6, lty = 1:5)

  #persp(x = 0:max(age), y = years, z = FM[1, , ], xlab = "Age", ylab = "Year", zlab = "F", theta = 45 + 180)

  # Abundance at age ----
  N <- local({
    naa <- array(NA, c(nsim, nage_om, nyears))
    for(x in 1:nsim) {
      naa[x, age + 1, ] <- t(admb_rep[["N"]])
    }

    # Backwards calculation for age 0 and 1 abundance
    for(a in 2:1) {
      naa[, a, 2:nyears - 1] <- naa[, a+1, 2:nyears] *
        exp(FM[, a, 2:nyears - 1] + M[, a, 2:nyears - 1]) # F+M should be numerically zero!
    }

    naa
  })

  # Maturity ogive ----
  Maturity <- local({
    mat <- array(NA, c(nsim, nage_om, nyears))
    for(x in 1:nsim) {
      mat[x, age + 1, ] <- admb_dat[["msex1"]] # No cambio ojiva
      if (tv_maturity) {
        mat[x, age + 1, -c(1:64)] <- t(admb_dat[["msex2"]])
      }
    }
    mat[, 1:2, ] <- 0
    mat
  })
  #matplot(years, t(Maturity[1, , ]), typ = "l", xlab = "Year", ylab = "Maturity", col = 1:6, lty = 1:5)
  #legend("topleft", legend = 0:max(age), col = 1:6, lty = 1:5)

  #plot(age, Maturity[1, age + 1, 1], typ = 'o', xlab = "Age", ylab = "Maturity")

  # Stock weight at age ----
  Weight_age <- local({
    wt <- array(NA, c(nsim, nage_om, nyears))
    for(x in 1:nsim) {
      wt[x, age + 1, ] <- t(admb_dat[["Wm"]])
    }
    wt[, 1:2, ] <- 0
    wt
  })
  #matplot(years, t(Weight_age[1, , ]), typ = "l", ylab = "Weight", col = 1:6, lty = 1:5)
  #legend("topleft", legend = 0:max(age), col = 1:6, lty = 1:5)

  # Length at age ----
  # The current ADMB model does not use length data
  # We need a placeholder for the operating model but we will have to ignore the length modeling
  Length_age <- local({
    laa <- array(NA, c(nsim, nage_om, nyears))
    for(x in 1:nsim) {
      laa[x, , ] <- 1:nage_om
    }
    laa
  })

  # Make operating model ----
  OM <- MSEtool::Assess2OM(
    Name = Name,
    proyears = proyears,
    CurrentYr = max(years),
    h = h,
    naa = N,
    faa = FM,
    waa = Weight_age,
    Mataa = Maturity,
    Maa = M,
    laa = Length_age,
    nyr_par_mu = 1,
    LowerTri = 2,
    R0 = R0,
    phi0 = phi0,
    SRrel = SRrel,
    spawn_time_frac = spawn_time_frac,
    Perr = Perr,
    AC = AC
  )

  # This is a placeholder for growth parameters that we will ignore ----
  OM@cpars$Linf <- OM@cpars$K <- OM@cpars$t0 <- rep(1, nsim)

  # Add acoustic index and its selectivity to operating model ----
  # See class?Data
  RealData <- new("Data")
  RealData@AddInd <- ifelse(admb_dat$indices[, 4] > 0, admb_dat$indices[, 4], NA) %>%
    array(c(1, 1, OM@nyears))
  RealData@CV_AddInd <- ifelse(admb_dat$indices[, 4] > 0, 0.15, NA) %>%
    array(c(1, 1, OM@nyears))
  RealData@AddIndV <- c(rep(0, 2), admb_rep$Sel.crucero[1, ]) %>%
    array(c(1, 1, OM@maxage + 1))
  RealData@AddIunits <- 1  # Biomass
  RealData@AddIndType <- 1

  # Add catch and catch at age to operating model ----
  RealData@Cat <- matrix(admb_dat$indices[, 3], 1, OM@nyears) # Adjusted catch
  RealData@CV_Cat <- matrix(0.05, 1, OM@nyears)

  RealData@CAA <- array(admb_dat$cap.arr, c(1, OM@nyears, OM@maxage + 1))

  # Pass data to OM ----
  OM@cpars$Data <- RealData

  # Specify that the index is proportional to the vulnerable biomass ----
  OM@cpars$AddIbeta <- matrix(1, OM@nsim, 1)

  # Specify the multinomial sample size of the fishery catch at age
  OM@CAA_ESS <- c(250, 250)

  # Specify the nominal sample size (expands values returned by dmultinom function)
  OM@CAA_nsamp <- c(3500, 3500)

  # Add adjustment factor for implementation error
  OM@TACFrac <- c(1.85, 1.85)

  # Save model ----
  return(OM)
}




#' @name make_hake_OM
#' @param jjdir Directory of the JJM assessment
#' @param model Character string for name of model
#' @export
make_hake_OM_AMAC <- function(
    jjdir = "admb-assessment/om2",
    model = "mc_0.0",
    Name = "Merluza comun - OM 2",
    nsim = 100,
    proyears = 36,
    AC = 0.65,
    Perr = 0.5
) {

  # Read ADMB report ----
  admb_rep <- read_hake_rep("admb-assessment/om1/modelo_140.rep")

  # Read ADMB data file ----
  admb_dat <- read_hake_dat("admb-assessment/om1/datos_14022.dat")

  # Read JJM report
  if (file.exists(file.path(jjdir, "jjmR_hake.rds"))) {
    res <- readRDS(file.path(jjdir, "jjmR_hake.rds"))
  } else {

    if (!"jjmR" %in% installed.packages()) {
      message("Installing jjmr package from Github..")
      remotes::install_github("SPRFMO/jjmr")
    }
    dir_cur <- getwd()
    on.exit(setwd(dir_cur))
    setwd(jjdir)
    jjm_results <- jjmR::readJJM(model, output = getwd(), input = getwd())
    res <- jjm_results[[1]]
    saveRDS(res, file = "jjmR_hake.rds")

  }

  # Historical parameters ----
  #nyears <- admb_dat$nanos
  #years <- admb_dat$indices[, 1] #1940:2022
  years <- seq(res$info$data$year[1], res$info$data$year[2])
  nyears <- length(years)

  # Operating model parameters ----
  age <- seq(res$info$data$age[1], res$info$data$age[2]) #2:13
  nage <- length(age)
  nage_om <- max(age) + 1

  spawn_time_frac <- 0.5833 # Fraction of year when spawning occurs (7/12)
  SRrel <- 2 # (1 = Beverton-Holt, 2 = Ricker)

  # Do not use these stock recruit parameters! ----
  #R0 <- exp(6.00663621455) # file jjms.par
  #SB0 <- res$output$Stock_1$msy_mt[, "bzero"][nyears]
  #phi0 <- SB0/R0
  #h <- 0.7

  SRR <- res$output$Stock_1$Stock_Rec
  #plot(SRR[, 1], SRR[, 2], xlab = "Year", ylab = "SSB", typ = 'o')
  #plot(SRR[, 1], SRR[, 4], xlab = "Year", ylab = "Recruitment", typ = 'o')
  #plot(SRR[, 2], SRR[, 4], xlab = "SSB", ylab = "Recruitment", typ = 'o', xlim = c(0, 700), ylim = c(0, 2000))
  #phi0 <- SRR[1, 2]/SRR[1, 4]
  #abline(a = 0, b = 1/phi0, col = 2)

  # Estimate stock recruit parameters from annual values
  R0_start <- log(mean(SRR[, 4]))
  h_start <- log(0.7 - 0.2)
  opt <- nlminb(c(R0_start, h_start), SAMtool:::get_SR, E = SRR[, 2], R = SRR[, 4], EPR0 = SRR[1, 2]/SRR[1, 4], type = "Ricker")
  SR_par <- SAMtool:::get_SR(opt$par, E = SRR[, 2], R = SRR[, 4], EPR0 = SRR[1, 2]/SRR[1, 4], opt = FALSE, figure = FALSE, type = "Ricker")

  R0 <- SR_par$R0
  SB0 <- SR_par$E0
  phi0 <- SR_par$E0/SR_par$R0
  h <- SR_par$h

  # Natural mortality ----
  M <- local({
    maa <- array(NA, c(nsim, nage_om, nyears))
    for(x in 1:nsim) {
      maa[x, -c(1, 2), ] <- t(res$output[[1]]$M)
    }
    maa[, 1:2, ] <- 1e-15 + .Machine$double.eps
    maa
  })

  # Fishing mortality ----
  F_fleet <- sapply(1:3, function(i) res$output[[1]][[paste0("F_age_", i)]][, -1], simplify = "array")
  F_age <- apply(F_fleet, 1:2, sum)

  FM <- local({
    faa <- array(NA, c(nsim, nage_om, nyears))
    for(x in 1:nsim) {
      faa[x, -c(1:2), ] <- t(F_age)
    }
    faa[, 1:2, ] <- 0
    faa
  })
  #matplot(years, t(FM[1, , ]), typ = "l", ylab = "F", col = 1:6, lty = 1:5)
  #legend("topleft", legend = 0:max(age), col = 1:6, lty = 1:5)

  #persp(x = 0:max(age), y = years, z = FM[1, , ], xlab = "Age", ylab = "Year", zlab = "F", theta = 45 + 180)

  # Sample random walk in selectivity in projection
  V_fleet <- apply(F_fleet, c(1, 3), function(x) x/max(x))
  V_all <- apply(F_age, 1, function(x) x/max(x))

  set.seed(324)
  Vproj <- sapply(1:nsim, function(...) {
    delta <- rnorm(proyears, 0, 0.025)
    delta2 <- rnorm(proyears * length(age), 0, 0.005) %>% matrix(proyears, length(age))

    Vlast <- F_age[nyears, ]/max(F_age[nyears, ])

    Vout <- array(NA, c(length(age), proyears))
    Vout[, 1] <- log(Vlast) + delta[1] + delta2[1, ]
    for(y in 2:proyears) Vout[, y] <- Vout[, y-1] + delta[y] + delta2[y, ]
    V <- exp(Vout)
    V[V > 1] <- 1
    V[11:12, ] <- 1
    return(V)
  }, simplify = "array")

  #matplot(Vproj[, , 5], typ = 'l')
  #persp(x = age, y = 2022 + 1:proyears, z = Vproj[, , 25], phi = 25, theta = -45)

  # Abundance at age ----
  N <- local({
    naa <- array(NA, c(nsim, nage_om, nyears))
    for(x in 1:nsim) {
      naa[x, age + 1, ] <- t(res$output[[1]]$N[, -1])
    }

    # Backwards calculation for age 0 and 1 abundance
    for(a in 2:1) {
      naa[, a, 2:nyears - 1] <- naa[, a+1, 2:nyears] *
        exp(FM[, a, 2:nyears - 1] + M[, a, 2:nyears - 1]) # F+M should be numerically zero!
    }

    naa
  })

  # Maturity ogive ----
  Pmat <- res$output$Stock_1$mature_a
  Wmat <- res$output$Stock_1$wt_a_pop
  #Pmat <- c(0, 0.07, 0.651, 0.967, 0.997, 1, 1, 1, 1, 1, 1, 1) # mc_0.0.ctl, line 171
  #Wmat <- c(0.182, 0.284,	0.420,	0.537,	0.647,	0.778,	0.861,	0.993,	1.063,	1.209,	1.363,	1.829) # mc_0.0.ctl, line 168
  Maturity <- local({
    mat <- array(NA, c(nsim, nage_om, nyears))
    for(x in 1:nsim) {
      mat[x, age + 1, ] <- Pmat # No cambio ojiva
    }
    mat[, 1:2, ] <- 0
    mat
  })
  #matplot(years, t(Maturity[1, , ]), typ = "l", xlab = "Year", ylab = "Maturity", col = 1:6, lty = 1:5)
  #legend("topleft", legend = 0:max(age), col = 1:6, lty = 1:5)

  #plot(age, Maturity[1, age + 1, 1], typ = 'o', xlab = "Age", ylab = "Maturity")

  # Stock weight at age ----
  Weight_age <- local({
    wt <- array(NA, c(nsim, nage_om, nyears))
    for(x in 1:nsim) {
      wt[x, age + 1, ] <- matrix(Wmat, nage, nyears)
    }
    wt[, 1:2, ] <- 0
    wt
  })
  #matplot(years, t(Weight_age[1, , ]), typ = "l", ylab = "Weight", col = 1:6, lty = 1:5)
  #legend("topleft", legend = 0:max(age), col = 1:6, lty = 1:5)

  # Length at age ----
  # The current ADMB model does not use length data
  # We need a placeholder for the operating model but we will have to ignore the length modeling
  Length_age <- local({
    laa <- array(NA, c(nsim, nage_om, nyears))
    for(x in 1:nsim) {
      laa[x, , ] <- 1:nage_om
    }
    laa
  })

  # Make operating model ----
  OM <- MSEtool::Assess2OM(
    Name = Name,
    proyears = proyears,
    CurrentYr = max(years),
    h = h,
    naa = N,
    faa = FM,
    waa = Weight_age,
    Mataa = Maturity,
    Maa = M,
    laa = Length_age,
    nyr_par_mu = 1,
    LowerTri = 2,
    R0 = R0,
    phi0 = phi0,
    SRrel = SRrel,
    spawn_time_frac = spawn_time_frac,
    Perr = Perr,
    AC = AC
  )

  # This is a placeholder for growth parameters that we will ignore ----
  OM@cpars$Linf <- OM@cpars$K <- OM@cpars$t0 <- rep(1, nsim)

  # Update selectivity in projection
  OM@cpars$V[, -c(1:2), nyears + 1:proyears] <- aperm(Vproj, c(3, 1, 2))

  # Weight at age of fishery catch
  Cobs <- sapply(1:3, function(i) res$output[[1]][[paste0("Obs_catch_", i)]])
  Z <- F_age + 0.33
  CN <- sapply(1:3, function(f) F_fleet[, , f]/Z * (1 - exp(-Z)) * t(N[1, -c(1:2), ]), simplify = "array")
  CB <- CN * res$data$Fwtatage
  Wt_age_C <- apply(CB, 1:2, sum)/apply(CN, 1:2, sum)
  OM@cpars$Wt_age_C <- local({
    out <- array(0, c(nsim, nage_om, nyears + proyears))
    out[, -c(1:2), 1:nyears] <- array(Wt_age_C, c(nyears, nage, nsim)) %>% aperm(3:1)
    out[, -c(1:2), nyears + 1:proyears] <- array(Wt_age_C[nyears, ], c(nage, proyears, nsim)) %>% aperm(c(3, 1, 2))
    out
  })

  # Add acoustic index and its selectivity to operating model ----
  # See class?Data
  RealData <- new("Data")
  RealData@AddInd <- ifelse(admb_dat$indices[, 4] > 0, admb_dat$indices[, 4], NA) %>%
    array(c(1, 1, OM@nyears))
  RealData@CV_AddInd <- ifelse(admb_dat$indices[, 4] > 0, 0.15, NA) %>%
    array(c(1, 1, OM@nyears))
  RealData@AddIndV <- c(rep(0, 2), admb_rep$Sel.crucero[1, ]) %>%
    array(c(1, 1, OM@maxage + 1))
  RealData@AddIunits <- 1  # Biomass
  RealData@AddIndType <- 1

  # Add catch and catch at age to operating model ----
  #RealData@Cat <- matrix(admb_dat$indices[, 3], 1, OM@nyears) # Adjusted catch
  #RealData@CV_Cat <- matrix(0.05, 1, OM@nyears)
  OM@Cbiascv <- 1e-8
  OM@Cobs <- c(0.01, 0.01)

  RealData@CAA <- array(admb_dat$cap.arr, c(1, OM@nyears, OM@maxage + 1))

  # Pass data to OM ----
  OM@cpars$Data <- RealData

  # Specify that the index is proportional to the vulnerable biomass ----
  OM@cpars$AddIbeta <- matrix(1, OM@nsim, 1)

  # Specify the multinomial sample size of the fishery catch at age
  OM@CAA_ESS <- c(250, 250)

  # Specify the nominal sample size (expands values returned by dmultinom function)
  OM@CAA_nsamp <- c(3500, 3500)

  # Add adjustment factor for implementation error
  OM@TACFrac <- c(1.85, 1.85)

  # Save model ----
  return(OM)
}
