

.get_array <- function(.mse, .name, type = c("SB", "B_BMSY", "F_FMSY", "CBA", "Index", "Catch", "R")) {
  type <- match.arg(type)

  x <- switch(
    type,
    "SB" = .mse@SSB,
    "B_BMSY" = .mse@SSB/.mse@RefPoint$ByYear$SSBMSY[, .mse@nyears],
    "F_FMSY" = .mse@FM/.mse@RefPoint$ByYear$FMSY[, .mse@nyears],
    "CBA" = .mse@TAC,
    "Catch" = .mse@Catch,
    "Index" = NULL,
    "R" = apply(.mse@N[, 3, , , ], 1:3, sum), # age 2
  )

  if (!is.null(x)) {
    x <- x %>%
      structure(dimnames = list(Simulation = 1:.mse@nsim, MP = .mse@MPs, Year = .mse@OM$CurrentYr[1] + seq(1, .mse@proyears))) %>%
      reshape2::melt()
  }

  xhist <- switch(
    type,
    "SB" = .mse@SSB_hist,
    "B_BMSY" = .mse@SSB_hist/.mse@RefPoint$ByYear$SSBMSY[, .mse@nyears],
    "F_FMSY" = .mse@FM_hist/.mse@RefPoint$ByYear$FMSY[, .mse@nyears],
    "Catch" = .mse@CB_hist,
    "CBA" = NULL,
    "Index" = NULL,
    "R" = apply(.mse@Hist@AtAge$Number[, 3, , ], 1:2, sum) # age 2
  )
  if (!is.null(xhist)) {
    xhist <- xhist %>%
      array(c(.mse@nsim, .mse@nyears, .mse@nMPs)) %>%
      aperm(c(1, 3, 2)) %>%
      structure(dimnames = list(Simulation = 1:.mse@nsim, MP = .mse@MPs, Year = .mse@OM$CurrentYr[1] - seq(.mse@nyears, 1) + 1)) %>%
      reshape2::melt()
  }

  if (type == "Index") {
    index_list <- lapply(1:.mse@nMPs, function(i) {
      Data <- .mse@PPD[[i]]

      Data@AddInd[, 1, ] %>%
        structure(dimnames = list(Simulation = 1:.mse@nsim, Year = Data@Year)) %>%
        reshape2::melt() %>%
        mutate(MP = .mse@MPs[i])
    })

    x <- do.call(rbind, index_list) %>%
      dplyr::filter(!is.na(value)) %>%
      mutate(MP = factor(MP, levels = .mse@MPs))
  }

  output <- rbind(x, xhist) %>%
    dplyr::mutate(OM = .name)

  return(output)
}

#' Plot time series figures
#'
#' Plots state variables for each operating model and management procedure
#'
#' @param MSE Either a list of [MSEtool::MSE-class] objects or single MSE object
#' @param names Character vector for the name of the operating model
#' @param sims Integer vector for a subset of simulations to plot
#' @param type Character for the state variable to plot
#' @returns A ggplot2 object
#' @importFrom stats median quantile
#' @importFrom methods is
#' @export
plot_array <- function(MSE, names = "", sims, type = c("SB", "B_BMSY", "F_FMSY", "CBA", "Index", "Catch", "R")) {
  type <- match.arg(type)

  if (is.list(MSE)) {
    #BMSY <- lapply(MSE_list, function(x) x@RefPoint$ByYear$SSBMSY[, x@nyears + 0:x@proyears])
    #matplot(t(BMSY[[1]]), typ = 'l')
    #BMSY <- lapply(MSE, function(x) c(x@RefPoint$ByYear$SSBMSY[1, x@nyears], x@RefPoint$SSBMSY[1]))

    #FMSY <- lapply(MSE_list, function(x) x@RefPoint$ByYear$FMSY[, x@nyears + 0:x@proyears])
    #FMSY <- lapply(MSE, function(x) c(x@RefPoint$ByYear$SSBMSY[1, x@nyears], x@RefPoint$SSBMSY[1]))

    array_df <- Map(.get_array, .mse = MSE, .name = names, type = type) %>%
      dplyr::bind_rows()

  } else if (is(MSE, "MSE")) {
    array_df <- .get_array(MSE, names, type)
  } else {
    stop("MSE list or object not found.")
  }

  ylab <- switch(type,
                 "SB" = "Biomasa desovante",
                 "B_BMSY" = expression(B/B[RMS]),
                 "F_FMSY" = expression(F/F[RMS]),
                 "CBA" = "CBA",
                 "Catch" = "Adjusted catch",
                 "Index" = "Crucero acústico",
                 "R" = "Reclutamiento")

  if (!missing(sims)) {

    array_sims <- filter(array_df, Simulation %in% sims) %>%
      mutate(Simulation = factor(Simulation))
    g <- ggplot(array_sims, aes(Year, value, group = Simulation, linetype = Simulation)) +
      facet_grid(vars(MP), vars(OM)) +
      geom_line() +
      #geom_hline(yintercept = 1, linetype = 3) +
      expand_limits(y = 0) +
      theme(panel.spacing = unit(0, "in"),
            legend.position = "bottom",
            axis.text.x = element_text(angle = 45, vjust = 0.5),
            strip.background = element_blank()) +
      coord_cartesian(expand = FALSE) +
      labs(x = "Año", y = ylab, linetype = "Simulación")

  } else {

    array_band <- array_df %>%
      group_by(Year, MP, OM) %>%
      summarise(med = median(value),
                lwr = quantile(value, 0.05),
                upr = quantile(value, 0.95),
                lwr2 = quantile(value, 0.25),
                upr2 = quantile(value, 0.75))


    g <- ggplot(array_band, aes(Year, med)) +
      facet_wrap(vars(MP)) +
      #geom_ribbon(fill = "grey90", aes(ymin = lwr, ymax = upr)) +
      #geom_ribbon(fill = "grey70", aes(ymin = lwr2, ymax = upr2)) +
      geom_ribbon(alpha = 0.1, aes(ymin = lwr2, ymax = upr2, fill = OM)) + # Plot only the interquartile range
      geom_line(aes(colour = OM)) +
      expand_limits(y = 0) +
      theme(panel.spacing = unit(0, "in"),
            axis.text.x = element_text(angle = 45, vjust = 0.5),
            legend.position = "bottom",
            strip.background = element_blank()) +
      coord_cartesian(expand = FALSE) +
      scale_fill_brewer(palette = "Dark2") +
      scale_colour_brewer(palette = "Dark2") +
      labs(x = "Año", y = ylab, fill = "MO", colour = "MO") +
      guides(fill = guide_legend(ncol = 2),
             colour = guide_legend(ncol = 2))
  }

  if (any(array_df$Year < 2022)) { # 2022 is the historical period
    g <- g + geom_vline(xintercept = 2022, linetype = 2)
  }

  g
}

#' @importFrom graphics matpoints legend abline
plot_HCR <- function(B_B0, F_FMSY) {
  Brel <- seq(0, 1, length.out = 200)
  Frel <- HCRlin(Brel, 0.2, 0.4, 0.75, 1)
  plot(
    Brel, Frel,
    xlab = expression(B/B[0]),
    ylab = expression(F/F[MSY]),
    ylim = c(0, 1.25), type = "l",
    lwd = 2,
    panel.first = {
      abline(v = c(0.2, 0.4), col = "red", lty = 3)
      abline(h = c(0.75, 1), col = "red", lty = 3)
    }
  )
  if (!missing(B_B0) && !missing(F_FMSY)) {
    matpoints(B_B0, F_FMSY, pch = 1, type = "p")
    legend(
      "bottomright",
      c("HCR", "Operating model"),
      lwd = c(2, 0),
      pch = c(NA, 1),
      col = c(1, 4)
    )
  }

  invisible()
}



#' Kobe time figure
#'
#' Plots proportion of simulations in each Kobe region by year and MP for a single operating model
#'
#' @param MSE [MSEtool::MSE-class] object
#' @param output Data frame for plotting. Used for averaging across operating models. (Optional)
#' @param type Character for biomass or fishing mortality figure
#' @returns A ggplot2 object
#' @export
plot_Kobe_time <- function(MSE, output, type = c("B", "F")) {
  type <- match.arg(type)

  if (!missing(MSE) && missing(output)) {

    if (type == "B") {

      B_BMSY <- .get_array(MSE, "OM", type = "B_BMSY") %>%
        filter(Year > MSE@OM$CurrentYr[1])

      Kgreen <- summarise(B_BMSY, value = mean(value >= 1.05), .by = c(MP, Year)) %>%
        mutate(val = "green", name = "Subexplotado")

      Kdarkgreen <- summarise(B_BMSY, value = mean(value < 1.05 & value >= 0.95), .by = c(MP, Year)) %>%
        mutate(val = "darkgreen", name = "Plena explotado")

      Kyellow <- summarise(B_BMSY, value = mean(value < 0.95 & value >= 0.5), .by = c(MP, Year)) %>%
        mutate(val = "yellow", name = "Sobreexplotado")

      Kred <- summarise(B_BMSY, value = mean(value < 0.5), .by = c(MP, Year)) %>%
        mutate(val = "red", name = "Agotado")

      output <- rbind(Kgreen, Kdarkgreen, Kyellow, Kred)

    } else {

      F_FMSY <- .get_array(MSE, "OM", type = "F_FMSY") %>%
        filter(Year > MSE@OM$CurrentYr[1])

      output <- summarise(F_FMSY, value = mean(value <= 1), .by = c(MP, Year)) %>%
        mutate(val = "grey60", name = "No sobre pesca")

    }
  }

  if (type == "B") {
    fill_values = c("Subexplotado" = "green",
                    "Plena explotado" = "darkgreen",
                    "Sobreexplotado" = "yellow",
                    "Agotado" = "red")

  } else {
    fill_values = c("No sobre pesca" = "grey60")
  }

  g <- output %>%
    mutate(name = factor(name, levels = names(fill_values))) %>%
    ggplot(aes(Year, value)) +
    geom_col(aes(fill = name), width = 1, colour = NA) +
    facet_wrap(vars(MP)) +
    scale_fill_manual(values = fill_values) +
    labs(x = "Año", y = "Probabilidad", fill = NULL) +
    theme(legend.position = "bottom",
          strip.background = element_blank()) +
    guides(fill = guide_legend(ncol = 2)) +
    coord_cartesian(expand = FALSE)

  g
}





.yield_curve <- function(mse, Fvec = seq(0.01, 1, 0.01), name = "OM", yr.ind = mse@nyears) {

  nyears <- mse@nyears
  #yr.ind <- mse@nyears
  x <- 1 # First simulation, yield curve should be identical in all simulations

  StockPars <- mse@Hist@SampPars$Stock
  Wt_age_C <- mse@Hist@SampPars$Fleet$Wt_age_C
  V <- mse@Hist@SampPars$Fleet$V

  maxage <- StockPars$maxage
  plusgroup <- StockPars$plusgroup
  hs <- StockPars$hs
  SSBpR <- StockPars$SSBpR
  spawn_time_frac <- StockPars$spawn_time_frac
  SRrel <- StockPars$SRrel
  R0 <- StockPars$R0

  if (length(yr.ind)==1) {
    M_at_Age <- StockPars$M_ageArray[x,,yr.ind]
    Wt_at_Age <- StockPars$Wt_age[x,, yr.ind]
    Fec_at_Age <- StockPars$Fec_Age[x,, yr.ind]
    Mat_at_Age <- StockPars$Mat_age[x,, yr.ind]
    V_at_Age <- V[x,, yr.ind]
    Wt_at_Age_C <- Wt_age_C[x,, yr.ind]
  } else {
    M_at_Age <- apply(StockPars$M_ageArray[x,,yr.ind], 1, mean)
    Wt_at_Age <- apply(StockPars$Wt_age[x,, yr.ind], 1, mean)
    Fec_at_Age <- apply(StockPars$Fec_Age[x,, yr.ind], 1, mean)
    Mat_at_Age <- apply(StockPars$Mat_age[x,, yr.ind], 1, mean)
    V_at_Age <- apply(V[x,, yr.ind], 1, mean)
    Wt_at_Age_C <- apply(Wt_age_C[x,, yr.ind], 1, mean)
  }

  # check for M = 0 in MOMs where maxage isn't the same for each stock
  if (max(which(M_at_Age!=0)) != (maxage+1)) {
    ind <- which(M_at_Age>0)
    M_at_Age <- M_at_Age[ind]
    Wt_at_Age <- Wt_at_Age[ind]
    Mat_at_Age <- Mat_at_Age[ind]
    Fec_at_Age <- Fec_at_Age[ind]
    V_at_Age <- V_at_Age[ind]
    Wt_at_Age_C <- Wt_at_Age_C[ind]

    maxage <- length(ind)-1
  }

  YC <- sapply(
    log(Fvec),
    MSEtool:::MSYCalcs,
    M_at_Age,
    Wt_at_Age,
    Mat_at_Age,
    Fec_at_Age,
    V_at_Age,
    Wt_at_Age_C,
    maxage,
    relRfun = StockPars$relRfun,
    SRRpars=StockPars$SRRpars[[x]],
    R0[x], SRrel[x], hs[x], SSBpR[x, 1], opt=2,
    plusgroup=plusgroup,
    spawn_time_frac=spawn_time_frac[x]
  )

  out <- as.data.frame(t(YC))
  out$OM <- name

  return(out)
}


.per_recruit <- function(mse, Fvec = seq(0.01, 1, 0.01), name = "OM", yr.ind = mse@nyears) {

  nyears <- mse@nyears
  #yr.ind <- mse@nyears
  x <- 1 # First simulation, yield curve should be identical in all simulations

  StockPars <- mse@Hist@SampPars$Stock
  Wt_age_C <- mse@Hist@SampPars$Fleet$Wt_age_C
  V <- mse@Hist@SampPars$Fleet$V

  maxage <- StockPars$maxage
  plusgroup <- StockPars$plusgroup
  hs <- StockPars$hs
  SSBpR <- StockPars$SSBpR
  spawn_time_frac <- StockPars$spawn_time_frac
  SRrel <- StockPars$SRrel
  R0 <- StockPars$R0

  if (length(yr.ind)==1) {
    M_at_Age <- StockPars$M_ageArray[x,,yr.ind]
    Wt_at_Age <- StockPars$Wt_age[x,, yr.ind]
    Fec_at_Age <- StockPars$Fec_Age[x,, yr.ind]
    Mat_at_Age <- StockPars$Mat_age[x,, yr.ind]
    V_at_Age <- V[x,, yr.ind]
    Wt_at_Age_C <- Wt_age_C[x,, yr.ind]
  } else {
    M_at_Age <- apply(StockPars$M_ageArray[x,,yr.ind], 1, mean)
    Wt_at_Age <- apply(StockPars$Wt_age[x,, yr.ind], 1, mean)
    Fec_at_Age <- apply(StockPars$Fec_Age[x,, yr.ind], 1, mean)
    Mat_at_Age <- apply(StockPars$Mat_age[x,, yr.ind], 1, mean)
    V_at_Age <- apply(V[x,, yr.ind], 1, mean)
    Wt_at_Age_C <- apply(Wt_age_C[x,, yr.ind], 1, mean)
  }

  # check for M = 0 in MOMs where maxage isn't the same for each stock
  if (max(which(M_at_Age!=0)) != (maxage+1)) {
    ind <- which(M_at_Age>0)
    M_at_Age <- M_at_Age[ind]
    Wt_at_Age <- Wt_at_Age[ind]
    Mat_at_Age <- Mat_at_Age[ind]
    Fec_at_Age <- Fec_at_Age[ind]
    V_at_Age <- V_at_Age[ind]
    Wt_at_Age_C <- Wt_at_Age_C[ind]

    maxage <- length(ind)-1
  }

  Ref_search <- MSEtool:::Ref_int_cpp(
    Fvec,
    M_at_Age = M_at_Age,
    Wt_at_Age = Wt_at_Age,
    Mat_at_Age = Mat_at_Age,
    Fec_at_Age=Fec_at_Age,
    V_at_Age = V_at_Age,
    Wt_at_Age_C = Wt_at_Age_C,
    StockPars$relRfun,
    StockPars$SRRpars[[x]],
    maxage = maxage,
    plusgroup = plusgroup,
    spawn_time_frac=StockPars$spawn_time_frac[x]
  )

  data.frame(
    F = Fvec,
    YPR = Ref_search[1, ],
    SPR = Ref_search[2, ],
    OM = name
  )
}

#' Yield curve
#'
#' Compare yield curve from set of operating models
#'
#' @param MSE_list List of [MSEtool::MSE-class] objects
#' @param names Character vector for the name of each operating model
#' @param Fvec Numeric vector for values of fishing mortality to calculate the yield curve
#' @param Fmaxplot Numeric, axis limit for F in the ggplot figure
#' @param yr.ind Integer, year of which to obtain the biological parameters and selectivity for the yield curve
#' @param by_OM Logical, whether to plot individual panels for each operating model
#' @return ggplot object
#' @export
#' @importFrom ggpubr ggarrange
plot_yield_curve <- function(MSE_list, names, Fvec = c(1e-8, seq(0.01, 3, 0.01)), Fmaxplot = 0.8,
                             yr.ind = MSE_list[[1]]@nyears, by_OM = FALSE) {

  RelRec <- SPR <- SB_SB0 <- NULL

  YC <- Map(.yield_curve, mse = MSE_list, name = names, MoreArgs = list(Fvec = Fvec, yr.ind = yr.ind)) %>%
    bind_rows() %>%
    filter(RelRec > 0)

  PR <- Map(.per_recruit, mse = MSE_list, name = names, MoreArgs = list(Fvec = Fvec, yr.ind = yr.ind)) %>%
    bind_rows()

  if (by_OM) {

    g <- ggplot(YC, aes(F, Yield, colour = OM)) +
      geom_line() +
      coord_cartesian(xlim = c(0, Fmaxplot)) +
      expand_limits(y = 0) +
      labs(colour = "MO") +
      facet_wrap(vars(OM))
    return(g)

  } else {

    g1 <- ggplot(YC, aes(F, Yield, colour = OM)) +
      geom_line() +
      coord_cartesian(xlim = c(0, Fmaxplot)) +
      expand_limits(y = 0) +
      labs(colour = "MO") +
      guides(colour = guide_legend(ncol = 2))

    g2 <- ggplot(PR, aes(F, SPR, colour = OM)) +
      geom_line() +
      coord_cartesian(xlim = c(0, Fmaxplot)) +
      expand_limits(y = 0) +
      labs(colour = "MO") +
      guides(colour = guide_legend(ncol = 2))

    g3 <- ggplot(YC, aes(SB_SB0, Yield, colour = OM)) +
      geom_line() +
      labs(x = expression(BD/BD[0]), colour = "MO") +
      guides(colour = guide_legend(ncol = 2))

    g <- ggpubr::ggarrange(g1, g3, g2, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

    return(g)
  }
}

