
# Simulate new index age composition
# Later on, we must pass to Rec@Misc in the MP to save sampled values
simIAA <- function(x, Data, N_OM) {
  #browser(expr = x == 3 && max(Data@Year) == 2030)

  if (is.null(Data@Misc[[x]]$IAA)) {
    proyears <- dim(Data@Misc$StockPars$N_P)[3]
    IAA <- matrix(NA, Data@MaxAge + 1, proyears) # n_age x projection years
  } else {
    IAA <- Data@Misc[[x]]$IAA
  }

  # Assumes annual survey
  styear <- which(!colSums(IAA, na.rm = TRUE))[1]
  endyear <- max(Data@Year) - Data@LHYear

  # Product of survey selectivity and abundance
  N_P <- Data@Misc$StockPars$N_P[x, , , ]
  if (is.null(N_P)) {
    n_y_rcm <- sum(Data@Year <= Data@LHYear)
    n_y_Data <- length(Data@Year)
    N_P <- N_OM[x, , seq(n_y_rcm + 1, n_y_Data) - n_y_rcm, ]
  }
  N_P <- apply(N_P, 1:2, sum)
  IAApred <- array(Data@AddIndV[x, 1, ] * N_P[, seq(styear, endyear), drop = FALSE],
                   c(Data@MaxAge + 1, length(styear:endyear), 1)) %>%
    aperm(c(3, 1, 2))

  # Update the master array
  IAA[, seq(styear, endyear)] <- MSEtool:::simCAA(
    nsim = 1, yrs = length(styear:endyear), n_age = Data@MaxAge + 1,
    Cret = IAApred, CAA_ESS = 500, CAA_nsamp = 60000
  )[1, , ] %>%
    matrix(length(styear:endyear), Data@MaxAge + 1) %>%
    t()

  return(IAA)
}


RCM2Assess <- function(RCModel) {

  .data <- RCModel@mean_fit$obj$env$data
  .parameters <- lapply(RCModel@mean_fit$obj$env$parameters, function(x) {
    y <- attr(x, "shape")
    if(is.null(y)) y <- x
    return(y)
  })

  .map <- RCModel@mean_fit$obj$env$map

  .stockpars <- MSEtool::SampleStockPars(RCModel@OM, msg = FALSE)
  .fleetpars <- MSEtool::SampleFleetPars(RCModel@OM, Stock = .stockpars, msg = FALSE)

  stockpar_names <- c("M_ageArray", "Len_age", "Linf", "LatASD",
                      "Wt_age", "Mat_age", "Fec_Age", "ageMarray",
                      "spawn_time_frac",
                      "SRrel", "hs", "R0", "procsd")
  .stockpars <- .stockpars[stockpar_names]

  fleetpar_names <- c("L5_y", "LFS_y", "Vmaxlen_y")
  .fleetpars <- .fleetpars[fleetpar_names]

  rcm_data <- RCModel@data

  force(.data)
  force(.parameters)
  force(.map)
  force(rcm_data)
  
  force(.stockpars)
  force(.fleetpars)

  fn <- function(x = 1, Data, StockPars, FleetPars, ...) {
    dots <- list(...)

    n_y_rcm <- rcm_data@Misc$nyears # sum(Data@Year <= Data@LHYear)
    n_y_Data <- length(Data@Year)

    new_data <- rcm_data

    #### Update estimation model with data and parameter variables that are time-dependent ----
    if (n_y_Data > n_y_rcm) {

      ## RCM data variables ----
      y_new <- Data@Year[seq(n_y_rcm + 1, n_y_Data)]
      n_y_new <- length(y_new)

      nfleet <- rcm_data@Misc$nfleet
      nsurvey <- rcm_data@Misc$nsurvey
      n_age <- rcm_data@Misc$maxage + 1

      # Append catch series from Data object
      Chist_new <- matrix(NA, n_y_new, nfleet)
      Chist_new[, 1] <- Data@Cat[x, seq(n_y_rcm + 1, n_y_Data)]
      new_data@Chist <- rbind(rcm_data@Chist, Chist_new)

      # Append SD of catch - extend time series from RCM
      C_sd_new <- matrix(rcm_data@C_sd[n_y_rcm, ], n_y_new, nfleet, byrow = TRUE)
      new_data@C_sd <- rbind(rcm_data@C_sd, C_sd_new)

      # Append effort (currently add zeros: MPs are naive to future squid mortality)
      Ehist_new <- matrix(0, n_y_new, nfleet)
      new_data@Ehist <- rbind(rcm_data@Ehist, Ehist_new)

      # Append catch at age from Data object
      CAA_new <- array(NA, c(n_y_new, n_age, nfleet))
      CAA_new[, , 1] <- matrix(Data@CAA[x, seq(n_y_rcm + 1, n_y_Data), ], n_y_new, n_age)
      new_data@CAA <- abind::abind(rcm_data@CAA, CAA_new, along = 1)

      # Append sample size - catch at age from RCM
      CAA_ESS_new <- matrix(rcm_data@CAA_ESS[n_y_rcm, ], n_y_new, nfleet, byrow = TRUE)
      new_data@CAA_ESS <- rbind(rcm_data@CAA_ESS, CAA_ESS_new)

      # Ignore catch at length, and CAL_ESS, length, mean size inputs

      # Append index series from Data object
      Index_new <- matrix(Data@AddInd[x, , seq(n_y_rcm + 1, n_y_Data)], nsurvey, n_y_new) %>% t()
      new_data@Index <- rbind(rcm_data@Index, Index_new)

      # Append SD of index from RCM
      I_sd_new <- matrix(rcm_data@I_sd[n_y_rcm, ], n_y_new, nsurvey, byrow = TRUE)
      new_data@I_sd <- rbind(rcm_data@I_sd, I_sd_new)

      # Index age composition
      IAA <- Data@Misc[[x]]$IAA
      if (is.null(IAA)) stop("Need index at age matrix, but cannot find it. Simulaton ", x, " , Year ", max(Data@Year))

      # Update the array for estimation model
      IAA_new <- IAA[, seq(1, n_y_new)] %>%
        array(c(Data@MaxAge + 1, n_y_new, 1)) %>%
        aperm(c(2, 1, 3))
      new_data@IAA <- abind::abind(rcm_data@IAA, IAA_new, along = 1)

      # Append sample size - index age from RCM
      IAA_ESS_new <- matrix(rcm_data@IAA_ESS[n_y_rcm, ], n_y_new, nsurvey, byrow = TRUE)
      new_data@IAA_ESS <- rbind(rcm_data@IAA_ESS, IAA_ESS_new)

      # Ignore index length composition

      # Append selectivity block - same as terminal RCM year
      sel_block_new <- matrix(rcm_data@sel_block[n_y_rcm, ], n_y_new, nfleet, byrow = TRUE)
      new_data@sel_block <- rbind(rcm_data@sel_block, sel_block_new)
    }

    ## Update the number of years for the RCM model in the MP ----
    new_data@Misc$nyears <- nyears <- n_y_Data
    new_data@Misc$CurrentYr <- max(Data@Year)

    ## Biological parameters in the MP will be identical to the operating model ----
    .new_data <- SAMtool::check_RCMdata(new_data, condition = new_data@Misc$condition, silent = TRUE) %>%
      getElement("RCMdata")

    # Update for map, start parameters, recruitment deviates
    new_start <- .parameters[c("vul_par", "ivul_par")]
    new_map <- .map[c("vul_par", "ivul_par")]
    new_map <- Map(function(x, y) array(as.integer(x), dim(y)), x = new_map, y = new_start)
    new_map$log_rec_dev <- local({
      m <- ifelse(1:n_y_Data <= 27 | 1:n_y_Data > n_y_Data - 2, NA, TRUE)
      m[!is.na(m)] <- 1:sum(m, na.rm = TRUE)
      m
    })

    ## Fit assessment model ----
    fit_rcm <- SAMtool:::RCM_est(
      x = 1, RCMdata = .new_data, selectivity = .data$vul_type, s_selectivity = .data$ivul_type,
      LWT = .new_data@Misc$LWT, comp_like = .data$comp_like, prior = .new_data@Misc$prior,
      max_F = .data$max_F,
      StockPars = .stockpars,
      FleetPars = .fleetpars,
      start = new_start, map = new_map,
      dots = list(vul_transform = FALSE)
    )
    #return(fit_rcm)

    ## Generate some assessment output ----
    obj <- fit_rcm$obj
    opt <- fit_rcm$opt
    SD <- fit_rcm$SD
    report <- fit_rcm$report
    conv <- SD$pdHess

    Year <- Data@Year
    Yearplusone <- YearR <- c(Year, max(Year) + 1)

    Assessment <- new("Assessment", Model = "RCM", Name = Data@Name, conv = conv,
                      h = report$h,
                      FMort = structure(apply(report$F_at_age, 1, max), names = Year), # Annual apical F
                      B = structure(report$B[1:nyears], names = Year),
                      SSB = structure(report$E[1:nyears], names = Year),
                      VB = structure(report$VB[1:nyears, 1], names = Year),
                      R = structure(report$R, names = YearR),
                      N = structure(rowSums(report$N), names = Yearplusone),
                      N_at_age = report$N,
                      Selectivity = report$F_at_age/apply(report$F_at_age, 1, max), # Aggregate selectivity
                      Dev = structure(report$log_rec_dev, names = Year),
                      Dev_type = "log-Recruitment deviations",
                      NLL = ifelse(is.character(opt), NA_real_, opt$objective),
                      obj = obj, opt = opt, SD = SD, TMB_report = report)

    ## Calculate reference points ----
    # Note that squid mortality is estimated as F and should assigned to natural mortality
    if (conv) {
      # For B0, need SRR, M = 0.33, fec from maturity
      # wt = Wcru, but we use weight-at-age from fishery
      #ref_pt <- SAMtool:::RCM_assess_ref(obj, report, yref = 1:n_y_Data)
      ref_pt_unfished <- SAMtool:::RCM_assess_ref(obj, report, yref = 1)
      ref_pt_MSY <- SAMtool:::RCM_assess_ref(obj, report, yref = n_y_Data)

      report$FMSY <- sapply(ref_pt_MSY, getElement, "FMSY")
      #tv_ref_pt <- length(unique(report$FMSY)) > 1

      report$MSY <- sapply(ref_pt_MSY, getElement, "MSY")
      report$VBMSY <- sapply(ref_pt_MSY, getElement, "VBMSY")
      report$RMSY <- sapply(ref_pt_MSY, getElement, "RMSY")
      report$BMSY <- sapply(ref_pt_MSY, getElement, "BMSY")
      report$EMSY <- sapply(ref_pt_MSY, getElement, "EMSY")

      refyear <- 1 #n_y_Data # Year of reference points for reporting

      report$new_B0 <- Assessment@B0 <- ref_pt_unfished[[refyear]]$new_B0
      report$new_E0 <- Assessment@SSB0 <- ref_pt_unfished[[refyear]]$new_E0
      report$new_VB0 <- Assessment@VB0 <- ref_pt_unfished[[refyear]]$new_VB0
      report$new_R0 <- Assessment@R0 <- ref_pt_unfished[[refyear]]$new_R0
      report$new_h <- Assessment@h <- ref_pt_unfished[[refyear]]$new_h

      Assessment@B_B0 <- Assessment@B/Assessment@B0
      Assessment@SSB_SSB0 <- Assessment@SSB/Assessment@SSB0
      Assessment@VB_VB0 <- Assessment@VB/Assessment@VB0

      Assessment@FMSY <- report$FMSY[refyear]
      Assessment@F_FMSY <- structure(Assessment@FMort/Assessment@FMSY, names = Year)

      Assessment@MSY <- report$MSY[refyear]
      Assessment@BMSY <- report$BMSY[refyear]
      Assessment@SSBMSY <- report$EMSY[refyear]
      Assessment@VBMSY <- report$VBMSY[refyear]
      Assessment@B_BMSY <- structure(Assessment@B/Assessment@BMSY, names = Year)
      Assessment@SSB_SSBMSY <- structure(Assessment@SSB/Assessment@SSBMSY, names = Year)
      #Assessment@VB_VBMSY <- structure(Assessment@VB/Assessment@VBMSY, names = Yearplusone)
      Assessment@TMB_report <- report

      catch_eq_fn <- function(Ftarget) {
        SAMtool:::catch_equation(method = "Baranov", Ftarget = Ftarget,
                                 M = obj$env$data$M_data[nyears, ],
                                 wt = obj$env$data$wt[nyears, ],
                                 N = report$N[nyears + 1, ],
                                 sel = report$F_at_age[nyears, ]/max(report$F_at_age[nyears, ]))
      }
      Assessment@forecast <- list(per_recruit = ref_pt_MSY[[1]][["per_recruit"]],
                                  catch_eq = catch_eq_fn)

    }
    return(Assessment)
  }
  class(fn) <- "Assessment"

  return(fn)
}

#hake_RCM <- readRDS("RCM/RCM_hake1_recdev.rds")
#MSE <- readRDS("OM_ADMB/hake_MSE1.rds")
#Hist <- readRDS("OM_ADMB/hake_Hist1.rds")

#m <- 3
#hake_assess <- RCM2Assess(hake_RCM)
#do_assess <- hake_assess(Data = MSE@PPD[[m]],
#                         N_OM = MSE@N[, , m, , ],
#                         StockPars = Hist@SampPars$Stock,
#                         FleetPars = Hist@SampPars$Fleet)
#
#do_assess <- hake_assess(Data = Hist@Data,
#                         N_OM = MSE@N[, , m, , ],
#                         StockPars = Hist@SampPars$Stock,
#                         FleetPars = Hist@SampPars$Fleet)
#do_assess@forecast$catch_eq(0.31)
