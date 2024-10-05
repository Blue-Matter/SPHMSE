

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
    "R" = apply(.mse@N[, 1, , , ], 1:3, sum),
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
    "R" = NULL
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
    matpoints(B_B0, F_FMSY, pch = 1, typ = "p")
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

      Kgreen <- structure(
        MSE@SB_SBMSY >= 1.05,
        dimnames = list(Simulation = 1:MSE@nsim, MP = MSE@MPs,
                        Year = 2022 + seq(1, MSE@proyears))) %>%
        apply(2:3, mean) %>%
        reshape2::melt() %>%
        mutate(val = "green", name = "Subexplotado")

      Kdarkgreen <- structure(
        MSE@SB_SBMSY < 1.05 & MSE@SB_SBMSY >= 0.95,
        dimnames = list(Simulation = 1:MSE@nsim, MP = MSE@MPs,
                        Year = 2022 + seq(1, MSE@proyears))) %>%
        apply(2:3, mean) %>%
        reshape2::melt() %>%
        mutate(val = "darkgreen", name = "Plena explotado")

      Kyellow <- structure(
        MSE@SB_SBMSY < 0.95 & MSE@SB_SBMSY >= 0.5,
        dimnames = list(Simulation = 1:MSE@nsim, MP = MSE@MPs,
                        Year = 2022 + seq(1, MSE@proyears))) %>%
        apply(2:3, mean) %>%
        reshape2::melt() %>%
        mutate(val = "yellow", name = "Sobreexplotado")

      Kred <- structure(
        MSE@SB_SBMSY < 0.5,
        dimnames = list(Simulation = 1:MSE@nsim, MP = MSE@MPs,
                        Year = 2022 + seq(1, MSE@proyears))) %>%
        apply(2:3, mean) %>%
        reshape2::melt() %>%
        mutate(val = "red", name = "Agotado")

      output <- rbind(Kgreen, Kdarkgreen, Kyellow, Kred)

    } else {

      output <- structure(
        MSE@F_FMSY <= 1,
        dimnames = list(Simulation = 1:MSE@nsim, MP = MSE@MPs,
                        Year = 2022 + seq(1, MSE@proyears))) %>%
        apply(2:3, mean) %>%
        reshape2::melt() %>%
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
