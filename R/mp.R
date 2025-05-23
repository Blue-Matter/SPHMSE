
#' @name MP
#' @title Hake management procedures
#' @description Functions of the various management procedures evaluated in the hake MSE (2023)
#'
#' @param x Integer, simulation number
#' @param Data [MSEtool::Data-class] object
#' @param ... Additional arguments (not used)
#' @returns [MSEtool::Rec-class] object
NULL

#' @name MP
#' @details
#' Sin_Pesca: No fishing reference MP, effort-based and thus will not be subject to CBA implementation error
#' @export
Sin_Pesca <- function(x, Data, reps = 1, ...) {
  rec <- new("Rec") # create recommendation object
  rec@Effort <- rep(1e-15, reps)
  rec
}
class(Sin_Pesca) <- "MP"
#MSEtool::NFref()

#' @name MP
#' @details
#' Manejo_Perfecto: Perfect management is an effort-based MP and thus will not be subject to CBA implementation error
#' MP specifies the fishing effort that corresponds to 95 percent FMSY
#' In openMSE, the effort is relative to the last historical year of the operating model.
#' Effort can also be more efficient/powerful in the projection
#' This code ensures that the correct fishing mortality is specified regardless of the
#' improvement in the fishing effort, based on `MSEtool::FMSYref()`
#' @export
Manejo_Perfecto <- function(x, Data, reps = 1, ...) {
  y <- max(Data@Year) - Data@LHYear+1
  nyears <- length(Data@Misc$FleetPars$Find[x,])
  FMSY <- Data@Misc$ReferencePoints$ByYear$FMSY[x,nyears+y]
  q <- Data@Misc$FleetPars$qs[x]
  qvar <- Data@Misc$FleetPars$qvar[x,y] # future only
  if (length(qvar)<1) qvar <- 1
  qinc <- Data@Misc$FleetPars$qinc[x] # future only
  qcur <- qvar * q*(1+qinc/100)^y # catchability this year

  HistE <- Data@OM$FinF[x] # Last historical fishing effort
  MSYE <- FMSY/qcur # effort for this year's FMSY

  Rec <- new('Rec')
  Rec@Effort <- rep(0.95 * MSYE/HistE, reps)
  Rec
}
class(Manejo_Perfecto) <- "MP"

# Examples of harvest control rules
SPH_A <- SAMtool::HCR_segment
formals(SPH_A)$OCP_type <- "SSB_SSB0"
formals(SPH_A)$Ftarget_type <- "FSPR"
formals(SPH_A)$SPR_targ <- 0.4
formals(SPH_A)$OCP <- c(0.2, 0.2)
formals(SPH_A)$relF <- c(0.6, 1)

SPH_C <- SPH_B <- SPH_A
formals(SPH_B)$relF <- c(0.75, 1)

formals(SPH_C)$OCP <- c(0.25, 0.5, 0.95)
formals(SPH_C)$relF <- c(0.1, 0.6, 0.95)

# Additional functions for stochastic MP code
catch_eq <- function(Ftarget, M, N, wt, sel) {
  Fage <- sel * Ftarget
  Z <- Fage + M
  CAA <- Fage/Z * (1 - exp(-Z)) * N
  sum(CAA * wt)
}

# See SAMtool:::linesegment to verify harvest control rule with low compliance
#BMSY <- 0.4 # BMSY = 0.4 B0
#B_B0 <- seq(0, 1, 0.01)
#B_BMSY <- B_B0/0.4

#F_FMSY <- SAMtool:::HCRlinesegment(B_B0, OCP = c(0.2, 0.2), relF = c(0.6, 1))
#plot(B_BMSY, F_FMSY, typ = 'l', xlab = expression(B/B[MSY]), ylab = expression(F/F[MSY]))
#SAMtool:::HCRlinesegment(0.2, OCP = c(0.2, 0.2), relF = c(0.6, 1))
#SAMtool:::HCRlinesegment(0.1999, OCP = c(0.2, 0.2), relF = c(0.6, 1))

# Subpesca proposed control rule
# Additional constraint in the annual change in catch to 15% if below estimated B_lim
#F_FMSY <- SAMtool:::HCRlinesegment(
#  B_BMSY,
#  OCP = c(0.25, 0.5, 0.95),
#  relF = c(0.1, 0.6, 0.95)
#)
#plot(B_BMSY, F_FMSY, typ = 'l',
#     lwd = 2,
#     ylim = c(0, 1),
#     xlab = expression(B/B[MSY]), ylab = expression(F/F[MSY]))
#abline(h = c(0.1, 0.6, 0.95), lty = 2)
#abline(v = c(0.25, 0.5, 0.95), lty = 2)


#source("98-RCM2Assess.R")
#hake_RCM <- readRDS("RCM/RCM_hake1.rds")
#hake_assess <- RCM2Assess(hake_RCM)

#' Make model-based management procedure
#'
#' Creates a management procedure from an assessment and harvest control rule
#'
#' @param HCR_fn A harvest control rule function that converts assessment output into catch advice, excluding hyper-rule
#' @returns A MP function that generates the CBA
#' @import MSEtool SAMtool
#' @importFrom mvtnorm rmvnorm
#' @export
make_hake_MP <- function(HCR_fn) {
  HCR_fn <- substitute(HCR_fn)
  delta <- substitute(delta)
  hake_assess <- substitute(hake_assess)

  fn_body <- bquote({

    # Sample IAA if needed
    if (max(Data@Year) > Data@LHYear) {
      Data@Misc[[x]]$IAA <- simIAA(x, Data)
    }

    # Fit assessment
    if (is.null(SPHenv$RCM_hake)) stop("No RCM model was found in SPHenv$RCM_hake")
    hake_assess <- RCM2Assess(SPHenv$RCM_hake)
    Mod <- hake_assess(x, Data)

    Rec <- new("Rec")
    Rec@Misc <- list(IAA = Data@Misc[[x]]$IAA)
    Rec@Misc <- c(Rec@Misc, SAMtool:::Assess_diagnostic(x, Data, Mod, include_assessment = TRUE))

    if (!Mod@conv) {
      Rec@TAC <- NA_real_
      return(Rec)
    }

    # Run harvest control rule to calculate the TAC
    if (reps > 1) {
      # Sample depletion from covariance matrix
      cov_matrix <- Mod@SD$cov.fixed

      set.seed(seed)
      samp_par <- mvtnorm::rmvnorm(reps, mean = Mod@SD$par.fixed, sigma = cov_matrix)
      samp_report <- lapply(1:reps, function(i) {
        Mod@obj$report(samp_par[i, ]) %>%
          SAMtool:::RCM_posthoc_adjust(obj = Mod@obj, par = samp_par[i, ])
      })

      nyears <- length(Mod@SSB)
      SSB <- sapply(samp_report, function(x) x$E[nyears])
      SSB0 <- sapply(samp_report, getElement, "E0_SR")
      B_B0 <- SSB/SSB0

      N <- sapply(samp_report, function(x) x$N[nrow(x$N), ])
      sel <- sapply(samp_report, function(x) {
        F_at_age <- x$F_at_age[nrow(x$F_at_age), ]
        F_at_age/max(F_at_age)
      })

      wt <- Mod@obj$env$data$wt
      wt <- wt[nrow(wt), ]

      M <- Mod@obj$env$data$M_data
      M <- M[nrow(M), ]

      F40 <- SAMtool:::get_FSPR(Mod@forecast$per_recruit$FM, Mod@forecast$per_recruit$SPR, target = 0.4)

      alpha <- SAMtool:::HCRlinesegment(B_B0, OCP = formals(.(HCR_fn))$OCP, relF = formals(.(HCR_fn))$relF)

      CBA <- sapply(1:reps, function(i) catch_eq(alpha[i] * F40, M = M, N = N[, i], wt = wt, sel = sel[, i]))

    } else if (reps == 1) {
      B_B0 <- Mod@SSB_SSB0[length(Mod@SSB_SSB0)]
      CBA <- .(HCR_fn)(Mod)@TAC
    }

    # Retrieve the previous CBA
    if (max(Data@Year) == Data@LHYear) {

      if (missing(CBA_previous)) {
        # If we are at the beginning of the projection, we specify the real 2022 CBA
        CBA_previous <- 41.4
      }
    } else {
      CBA_previous <- Data@MPrec[x]
    }

    CBA_range <- delta * CBA_previous

    CBA[CBA < min(CBA_range)] <- min(CBA_range)
    CBA[CBA > max(CBA_range)] <- max(CBA_range)

    if (verbose) {
      message("Range of B/B0: ", paste(round(range(B_B0), 2), collapse = " - "))
      message("Maximum CBA range (delta * CBA_previous): ", paste(round(range(CBA_range), 2), collapse = " - "))
    }

    Rec@TAC <- CBA
    return(Rec)
  })

  hake_MP <- eval(
    call("function",
         as.pairlist(alist(x = 1, Data = , reps = 1, delta = c(0.85, 1.15), CBA_previous = , seed = 1, verbose = FALSE)),
         fn_body)
  )
  return(structure(hake_MP, class = "MP"))
}



#' @importFrom stats loess sd
make_hake_interim_MP <- function(HCR_fn) { # Compare with SAMtool::make_interim_MP
  HCR_fn <- substitute(HCR_fn)
  delta <- substitute(delta)
  hake_assess <- substitute(hake_assess)

  fn_body <- bquote({

    type <- match.arg(type)

    # Sample IAA if needed
    if (max(Data@Year) > Data@LHYear) {
      Data@Misc[[x]]$IAA <- simIAA(x, Data)
    }
    Rec <- new("Rec")
    Rec@TAC <- rep(NA_real_, reps)
    Rec@Misc <- list(IAA = Data@Misc[[x]]$IAA)

    # Setup interim procedure
    ny <- length(Data@Year)
    Current_Yr <- Data@Year[ny]
    first_assess_yr <- Current_Yr == Data@LHYear
    run_assessment <- first_assess_yr || Current_Yr == Data@Misc[[x]]$interim$next_assess_yr
    run_interim <- !run_assessment

    if (run_assessment) {
      # Fit assessment
      if (is.null(SPHenv$RCM_hake)) stop("No RCM model was found in SPHenv$RCM_hake")
      hake_assess <- RCM2Assess(SPHenv$RCM_hake)
      Mod <- hake_assess(x, Data)

      Rec@Misc <- c(Rec@Misc, SAMtool:::Assess_diagnostic(x, Data, Mod, include_assessment = TRUE))

      if (Mod@conv) {  # Assessment converged. Run the HCR and report parameters for future interim procedure

        # Run harvest control rule to calculate the TAC
        if (reps > 1) {
          # Sample depletion from covariance matrix
          cov_matrix <- Mod@SD$cov.fixed

          set.seed(seed)
          samp_par <- mvtnorm::rmvnorm(reps, mean = Mod@SD$par.fixed, sigma = cov_matrix)
          samp_report <- lapply(1:reps, function(i) {
            Mod@obj$report(samp_par[i, ]) %>%
              SAMtool:::RCM_posthoc_adjust(obj = Mod@obj, par = samp_par[i, ])
          })

          nyears <- length(Mod@SSB)
          SSB <- sapply(samp_report, function(x) x$E[nyears])
          SSB0 <- sapply(samp_report, getElement, "E0_SR")
          B_B0 <- SSB/SSB0

          N <- sapply(samp_report, function(x) x$N[nrow(x$N), ])
          sel <- sapply(samp_report, function(x) {
            F_at_age <- x$F_at_age[nrow(x$F_at_age), ]
            F_at_age/max(F_at_age)
          })

          wt <- Mod@obj$env$data$wt
          wt <- wt[nrow(wt), ]

          M <- Mod@obj$env$data$M_data
          M <- M[nrow(M), ]

          F40 <- SAMtool:::get_FSPR(Mod@forecast$per_recruit$FM, Mod@forecast$per_recruit$SPR, target = 0.4)

          alpha <- SAMtool:::HCRlinesegment(B_B0, OCP = formals(.(HCR_fn))$OCP, relF = formals(.(HCR_fn))$relF)

          CBA <- sapply(1:reps, function(i) catch_eq(alpha[i] * F40, M = M, N = N[, i], wt = wt, sel = sel[, i]))

        } else if (reps == 1) {
          B_B0 <- Mod@SSB_SSB0[length(Mod@SSB_SSB0)]
          CBA <- .(HCR_fn)(Mod)@TAC
        }

        Rec@TAC <- CBA
        Rec@Misc$interim <- list(Cref = CBA, # A vector of reps long
                                 Iref = Mod@Index[nrow(Mod@Index), 1],
                                 next_assess_yr = Current_Yr + assessment_interval)
        if (type == "buffer") {
          Rec@Misc$interim$s <- log(Mod@Obs_Index[, 1]/Mod@Index[, 1]) %>% sd(na.rm = TRUE)
        } else {
          Rec@Misc$interim$s <- 0
        }

      } else {
        # Assessment did not converge. Try the assessment again next year
        next_assess_yr <- Current_Yr + 1

        if (first_assess_yr || is.null(Data@Misc[[x]]$interim)) { # No assessment has been run in the past, can't do interim procedure, return TAC = NA

          Rec@TAC <- rep(NA_real_, reps)

          if (!is.null(Mod@info$Misc)) Rec@Misc <- c(Rec@Misc, Mod@info$Misc)

          Rec@Misc$interim <- list(Cref = NA_real_, Iref = NA_real_, next_assess_yr = next_assess_yr, s = 0)
          return(Rec)

        } else { # There should be a previous assessment to continue the interim procedure for the current year
          run_interim <- TRUE
        }
      }

    }

    if (run_interim) {
      Cref <- Data@Misc[[x]]$interim$Cref
      Iref <- Data@Misc[[x]]$interim$Iref
      I_y <- switch(
        type,
        "buffer" = Data@AddInd[x, 1, ny],
        "mean" =  local({
          if (is.null(type_par)) type_par <- 3
          mean(Data@AddInd[x, 1, seq(ny - type_par + 1, ny)], na.rm = TRUE)
        }),
        "loess" = local({
          I_df <- data.frame(Year = Data@Year, Ind = Data@AddInd[x, 1, 1:ny])
          if (is.null(type_par)) type_par <- formals(loess)$span
          fit <- loess(Ind ~ Year, data = I_df, span = type_par)
          fit$fitted[length(fit$fitted)]
        }),
        "none" = Data@AddInd[x, 1, ny]
      )

      if (type == "buffer") {
        if (is.null(type_par)) {
          b <- 1
        } else {
          b <- type_par
        }
        s <- Data@Misc[[x]]$interim$s
      } else {
        b <- s <- 0
      }

      #Rec <- new("Rec")
      TAC <- Cref * (I_y + b * s)/(Iref + b * s)
      if (is.null(TAC)) TAC <- NA_real_
      Rec@TAC <- TAC

      Rec@Misc <- c(Rec@Misc, Data@Misc[[x]][setdiff(names(Data@Misc[[x]]), names(Rec@Misc))]) # Only add items that are missing from Rec@Misc
      if (exists("next_assess_yr", inherits = FALSE)) Rec@Misc$interim$next_assess_yr <- next_assess_yr
    }

    # Apply the hyper rule last
    if (any(!is.na(Rec@TAC))) {
      CBA_hyperrule <- Rec@TAC

      # Retrieve the previous CBA
      if (max(Data@Year) == Data@LHYear) {

        if (missing(CBA_previous)) {
          # If we are at the beginning of the projection, we specify the real 2022 CBA
          CBA_previous <- 41.4
        }
      } else {
        CBA_previous <- Data@MPrec[x]
      }

      CBA_range <- delta * CBA_previous

      CBA_hyperrule[Rec@TAC < min(CBA_range)] <- min(CBA_range)
      CBA_hyperrule[Rec@TAC > max(CBA_range)] <- max(CBA_range)

      Rec@TAC <- CBA_hyperrule
    }
    return(Rec)
  })

  hake_MP <- eval(
    call("function",
         as.pairlist(alist(x = 1, Data = , reps = 1, delta = c(0.85, 1.15), CBA_previous = , seed = 1,
                           assessment_interval = 3,
                           type = c("buffer", "mean", "loess", "none"),
                           type_par = NULL,
                           verbose = FALSE)),
         fn_body)
  )
  return(structure(hake_MP, class = "MP"))
}



#' Environment for storing hake assessment object
#'
#' An environment that stores variables for the hake assessment.
#' @export
SPHenv <- new.env(parent = emptyenv())
load("data/RCM_hake_2023.rda", envir = SPHenv)
SPHenv$RCM_hake <- SPHenv$RCM_hake_2023

#' @name RCM_hake_2023
#' @title Hake assessment 2023
#' @description TMB model fitted to common hake in 2023 to replicate ADMB assessment
#' @examples
#' data(RCM_hake_2023)
NULL

#' @name MP
#' @details PM_A fits an assessment model and applies the low-compliance control rule for the CBA
#' @param reps Integer, number of stochastic replicates for generating the advice. Not defined for index-based MPs.
#' @section Model-based MPs:
#' Model-based MPs use the RCM model (coded in TMB). By default, the RCM object developed from data to 2023 is located in `SPHenv$RCM_hake`.
#' To update the model for the MP, replace `SPHenv$RCM_hake` with a new `RCModel-class` object.
#' @export
PM_A <- make_hake_MP(SPH_A)

#' @name MP
#' @details PM_B fits an assessment model and applies the high-compliance control rule for the CBA
#' @export
PM_B <- make_hake_MP(SPH_B)

#' @name MP
#' @details PM_C fits an assessment model and applies the ramped harvest control rule for the CBA
#' @export
PM_C <- make_hake_MP(SPH_C)


#' @name MP
#' @details PM_A_interim fits an assessment model once every three years with the low-compliance control rule for the CBA.
#' In the intervening years, an adjustment is made based on new index values. See details below. In each year, the hyperrule is
#' applied at the last step.
#' @param assessment_interval Integer, frequency of assessment for the interim MPs (3 years is the default)
#' @param type Character, how the index is used to calculate the CBA in the interim procedure. See details.
#' @param type_par A control parameter for the interim procedure. See details.
#' @section Interim MP:
#' The interim procedure updates the CBA according to index observations between periodic
#' assessment intervals. **Always ensure to set:** `OM@@interval <- 1`. The actual assessment frequency is specified in function argument for each
#' management procedure (argument `assessment_interval`).
#'
#' In the year when the assessment is applied, the CBA is set by fitting the model and then running the harvest control rule.
#'
#' Between assessments, the CBA is updated as
#' \deqn{
#' \textrm{CBA}_{y+1} = C_{\textrm{ref}} (I_y + b \times s)/(I_{\textrm{ref}} + b \times s)
#' }
#' where `Cref` is the CBA calculated from the most recent assessment, `Iref` is the value of the index when `Cref` was calculated
#' (see Equations 6 and 7 of Huynh et al. 2020). The value of `I_y` depends on `type`, with `b` and `s` equal zero unless
#' `type = "buffer"`:
#'
#' \itemize{
#' \item `"buffer"` - `I_y` is the most recent index with `b` is specifed by `type_par` (default = 1), and `s` is
#' the standard deviation of index residuals from the most recent assessment.
#' \item `"mean"` - `I_y` is the mean value of the index over the most recent `type_par` years (default = 3).
#' \item `"loess"` - `I_y` is the most recent index predicted by a [loess][stats::loess] smoother applied over the entire time series of the index.
#' Use `type_par` to adjust the `span` parameter (default = 0.75).
#' \item `"none"` - `I_y` is the most recent index. Index values are not adjusted in the interim procedure.
#' }
#'
#' @references
#' Huynh et al. 2020. The interim management procedure approach for assessed stocks: Responsive management advice and lower assessment
#' frequency. Fish Fish. 21:663–679. \doi{10.1111/faf.12453}
#' @export
PM_A_interim <- make_hake_interim_MP(SPH_A)

#' @name MP
#' @details PM_B_interim fits an assessment model once every three years with the high-compliance control rule for the CBA.
#' In the intervening years, an adjustment is made based on new index values. See details below.
#' @export
PM_B_interim <- make_hake_interim_MP(SPH_B)

#' @name MP
#' @details PM_A_interim fits an assessment model once every three years with the ramped control rule for the CBA.
#' In the intervening years, an adjustment is made based on new index values. See details below.
#' @export
PM_C_interim <- make_hake_interim_MP(SPH_C)

#PM_Actual <- function(x, Data, reps = 1, delta = c(0.85, 1.15), HCR_fn = SPH_A, ...) {
#
#  # Fit assessment
#  Mod <- hake_assess(x, Data)
#
#  # Run harvest control rule to calculate the TAC
#  CBA <- HCR_fn(Mod, ...)@TAC
#
#  # Retrieve the previous CBA
#  if (max(Data@Year) == Data@LHYear) {
#    # If we are at the beginning of the projection, we need to specify the real CBA
#    # This is just a placeholder (software uses the last historical catch as the real CBA)
#    CBA_previous <- Data@MPrec[x]
#  } else {
#    CBA_previous <- Data@MPrec[x]
#  }
#  CBA_range <- delta * CBA_previous
#
#  if (CBA < min(CBA_range)) CBA <- min(CBA_range)
#  if (CBA > max(CBA_range)) CBA <- max(CBA_range)
#
#  Rec <- new("Rec")
#  Rec@TAC <- CBA
#  Rec@Misc <- list(IAA = Mod@info$IAA)
#  Rec@Misc <- c(Rec@Misc, SAMtool:::Assess_diagnostic(x, Data, Mod, include_assessment = TRUE))
#  return(Rec)
#}
#class(PM_Actual) <- "MP"

#' @name MP
#' @details I3 is an index-based management procedure proposed by SSPA (slope over 3 years)
#' @param y Integer, length of years for calculating the trend, i.e., slope, in the index
#' @param lambda Vector length 2, change in CBA relative to slope. First value if slope < 0, second value if lambda > 0
#' @param delta Vector length 2, minimum and maximum change in CBA (applied at the last step of the CBA calculation). Set to `c(0, Inf)` to remove the hyperrule.
#' @param CBA_previous Numeric, the CBA in the previous year, used to calculate the CBA for the following year. Only used if `max(Data@Year) == Data@LHYear`.
#' If not provided, uses a value of 41.4 t (2022 CBA).
#' @param verbose Logical, whether to report intermediary steps for the CBA calculation to console
#' @param seed Integer to set the random number generator when sampling the covariance matrix (when `reps` > 1)
#' @importFrom stats lm coef
#' @export
I3 <- function(x = 1, Data, y = 3, lambda = c(2, 1), delta = c(0.9, 1.1), CBA_previous, ...) {

  # Index series
  Ind <- Data@AddInd[x, 1, ]
  n_y <- length(Ind)

  # Calculate the slope of the natural logarithm
  Ind_slope <- Ind[seq(n_y - y + 1, n_y)]
  anos <- seq(1, y)
  reg <- stats::lm(log(Ind_slope) ~ anos)
  slope <- coef(reg)["anos"]

  # Calculate lambda and the relative change in CBA
  lam <- ifelse(slope < 0, lambda[1], lambda[2])

  rate <- 1 + lam * slope
  if (rate < delta[1]) rate <- delta[1]
  if (rate > delta[2]) rate <- delta[2]

  # Calculate the new CBA
  if (max(Data@Year) == Data@LHYear) {

    if (missing(CBA_previous)) {
      # If we are at the beginning of the projection, we specify the real 2022 CBA
      CBA_previous <- 41.4
    }

  } else {
    CBA_previous <- Data@MPrec[x]
  }

  CBA <- CBA_previous * rate

  Rec <- new("Rec")
  Rec@TAC <- CBA
  return(Rec)
}
class(I3) <- "MP"

#' @name MP
#' @details I5 is an index-based management procedure proposed by SSPA (slope over 5 years)
#' @export
I5 <- I3
formals(I5)$y <- 5
class(I5) <- "MP"

#' @name MP
#' @details I3_lambda1 is an index-based management procedure proposed by SSPA (slope over 3 years with symmetric lambda)
#' @export
I3_lambda1 <- I3
formals(I3_lambda1)$lambda <- c(1, 1)
class(I3_lambda1) <- "MP"

#' @name MP
#' @details I5_lambda1 is an index-based management procedure proposed by SSPA (slope over 5 years with symmetric lambda)
#' @export
I5_lambda1 <- I3
formals(I5_lambda1)$lambda <- c(1, 1)
formals(I5_lambda1)$y <- 5
class(I5_lambda1) <- "MP"
