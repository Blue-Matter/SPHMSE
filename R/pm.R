
#' @import dplyr ggplot2
NULL

.Kobe <- function(MSEobj, Yrs = c(1, 5), Ref = 0.95) {
  Yrs <- ChkYrs(Yrs, MSEobj)
  y_vector <- seq(Yrs[1], Yrs[2])

  PMobj <- new("PMobj")
  #PMobj@Name <- Name
  #PMobj@Caption <- Caption

  # See articles on reference points: https://openmse.com/tutorial-reference-points/
  #BMSY <- MSEobj@RefPoint$ByYear$SSBMSY
  #BMSY <- MSEobj@RefPoint$SSBMSY
  #B_BMSY <- MSEobj@SSB/BMSY
  B_BMSY <- MSEobj@SB_SBMSY

  PMobj@Ref <- Ref
  tt <- B_BMSY[, , y_vector] >= Ref
  if (is.null(dim(tt)))
    tt <- matrix(tt, nrow=MSEobj@nsim, ncol=1)
  PMobj@Stat <- tt
  PMobj@Prob <- calcProb(PMobj@Stat, MSEobj)
  PMobj@Mean <- calcMean(PMobj@Prob)
  PMobj@MPs <- MSEobj@MPs

  return(PMobj)
}

#' @rdname PM
#' @title Performance metrics
#'
#' @description Functions to calculate hake performance metrics
#' @param MSEobj [MSEtool::MSE-class] object
#' @returns [MSEtool::PMobj-class] object
#' @details
#' ZV_CP - Probability Kobe Green in Years 1 - 5
#' @export
ZV_CP <- function(MSEobj) {
  if(!inherits(MSEobj, 'MSE')) stop('This PM method is designed for objects of class `MSE`')
  PM <- .Kobe(MSEobj, Ref = 0.95, Yrs = c(1, 5))
  PM@Name <- PM@Caption <- "Probability Kobe Green in Years 1 - 5"
  PM
}
class(ZV_CP) <- "PM"

#' @name PM
#' @details
#' ZV_MP - Probability Kobe Green in Years 12 - 15
#' @export
ZV_MP <- function(MSEobj) {
  if(!inherits(MSEobj, 'MSE')) stop('This PM method is designed for objects of class `MSE`')
  PM <- .Kobe(MSEobj, Ref = 0.95, Yrs = c(12, 15))
  PM@Name <- PM@Caption <- "Probability Kobe Green in Years 12 - 15"
  PM
}
class(ZV_MP) <- "PM"

#' @name PM
#' @details
#' ZV_LP - Probability Kobe Green in Years 24 - 36
#' @export
ZV_LP <- function(MSEobj) {
  if(!inherits(MSEobj, 'MSE')) stop('This PM method is designed for objects of class `MSE`')
  PM <- .Kobe(MSEobj, Ref = 0.95, Yrs = c(24, 36))
  PM@Name <- PM@Caption <- "Probability Kobe Green in Years 24 - 36"
  PM
}
class(ZV_LP) <- "PM"

#' @name PM
#' @details
#' NZR_CP - Probability of avoiding Kobe Red in Years 1 - 5
#' @export
NZR_CP <- function(MSEobj) {
  if(!inherits(MSEobj, 'MSE')) stop('This PM method is designed for objects of class `MSE`')
  PM <- .Kobe(MSEobj, Ref = 0.5, Yrs = c(1, 5))
  PM@Name <- PM@Caption <- "Probability of avoiding Kobe Red in Years 1 - 5"
  PM
}
class(NZR_CP) <- "PM"

#' @name PM
#' @details
#' NZR_MP - Probability of avoiding Kobe Red in Years 12 - 15
#' @export
NZR_MP <- function(MSEobj) {
  if(!inherits(MSEobj, 'MSE')) stop('This PM method is designed for objects of class `MSE`')
  PM <- .Kobe(MSEobj, Ref = 0.5, Yrs = c(12, 15))
  PM@Name <- PM@Caption <- "Probability of avoiding Kobe Red in Years 12 - 15"
  PM
}
class(NZR_CP) <- "PM"

#' @name PM
#' @details
#' NZR_LP - Probability of avoiding Kobe Red in Years 24 - 36
#' @export
NZR_LP <- function(MSEobj) {
  if(!inherits(MSEobj, 'MSE')) stop('This PM method is designed for objects of class `MSE`')
  PM <- .Kobe(MSEobj, Ref = 0.5, Yrs = c(24, 36))
  PM@Name <- PM@Caption <- "Probability of avoiding Kobe Red in Years 24 - 36"
  PM
}
class(NZR_LP) <- "PM"


.NSP <- function(MSEobj, Yrs = c(1, 5), Ref = 1) {
  Yrs <- ChkYrs(Yrs, MSEobj)
  y_vector <- seq(Yrs[1], Yrs[2])

  PMobj <- new("PMobj")
  #PMobj@Name <- Name
  #PMobj@Caption <- Caption

  # See articles on reference points: https://openmse.com/tutorial-reference-points/
  #FMSY <- MSEobj@RefPoint$ByYear$FMSY
  #FMSY <- MSEobj@RefPoint$FMSY
  #F_FMSY <- MSEobj@FM/FMSY
  F_FMSY <- MSEobj@F_FMSY

  PMobj@Ref <- Ref
  tt <- F_FMSY[, , y_vector] <= Ref
  if (is.null(dim(tt)))
    tt <- matrix(tt, nrow=MSEobj@nsim, ncol=1)
  PMobj@Stat <- tt
  PMobj@Prob <- calcProb(PMobj@Stat, MSEobj)
  PMobj@Mean <- calcMean(PMobj@Prob)
  PMobj@MPs <- MSEobj@MPs

  return(PMobj)
}

#' @name PM
#' @details
#' NSP_CP - Probability of avoiding overfishing in Years 1 - 5
#' @export
NSP_CP <- function(MSEobj) {
  if(!inherits(MSEobj, 'MSE')) stop('This PM method is designed for objects of class `MSE`')
  PM <- .NSP(MSEobj, Yrs = c(1, 5))
  PM@Name <- PM@Caption <- "Probability of avoiding overfishing in Years 1 - 5"
  PM
}
class(NSP_CP) <- "PM"


#' @name PM
#' @details
#' NSP_MP - Probability of avoiding overfishing in Years 12 - 15
#' @export
NSP_MP <- function(MSEobj) {
  if(!inherits(MSEobj, 'MSE')) stop('This PM method is designed for objects of class `MSE`')
  PM <- .NSP(MSEobj, Yrs = c(12, 15))
  PM@Name <- PM@Caption <- "Probability of avoiding overfishing in Years 12 - 15"
  PM
}
class(NSP_MP) <- "PM"

#' @name PM
#' @details
#' NSP_LP - Probability of avoiding overfishing in Years 24 - 36
#' @export
NSP_LP <- function(MSEobj) {
  if(!inherits(MSEobj, 'MSE')) stop('This PM method is designed for objects of class `MSE`')
  PM <- .NSP(MSEobj, Yrs = c(24, 36))
  PM@Name <- PM@Caption <- "Probability of avoiding overfishing in Years 24 - 36"
  PM
}
class(NSP_LP) <- "PM"

#' @importFrom methods new
.CBA <- function(MSEobj, Yrs = c(1, 5), Ref = 20, type = c("min", "mean", "var")) {
  type <- match.arg(type)

  Yrs <- ChkYrs(Yrs, MSEobj)
  y_vector <- seq(Yrs[1], Yrs[2])

  PMobj <- new("PMobj")
  #PMobj@Name <- Name
  #PMobj@Caption <- Caption

  if (type == "min") {
    PMobj@Ref <- Ref
    tt <- MSEobj@TAC[, , y_vector] >= Ref
  }

  if (type == "mean") {
    PMobj@Ref <- 0
    tt <- MSEobj@TAC[, , y_vector]
  }

  if (type == "var") {
    PMobj@Ref <- 0
    y2 <- y_vector[-1]
    y1 <- y_vector[-length(y_vector)]
    tt <- abs(MSEobj@TAC[, , y2]/MSEobj@TAC[, , y1] - 1)
  }

  if (is.null(dim(tt)))
    tt <- matrix(tt, nrow=MSEobj@nsim, ncol=1)
  PMobj@Stat <- tt
  PMobj@Prob <- calcProb(PMobj@Stat, MSEobj)
  PMobj@Mean <- calcMean(PMobj@Prob)
  PMobj@MPs <- MSEobj@MPs

  return(PMobj)
}


#' @name PM
#' @details
#' CBAmin_CP - Probability CBA > 20 kt in Years 1-5
#' @export
CBAmin_CP <- function(MSEobj) {
  if(!inherits(MSEobj, 'MSE')) stop('This PM method is designed for objects of class `MSE`')
  PM <- .CBA(MSEobj, Yrs = c(1, 5))
  PM@Name <- PM@Caption <- "Probability CBA > 20 kt in Years 1-5"
  PM
}
class(CBAmin_CP) <- "PM"

#' @name PM
#' @details
#' CBAmin_MP - Probability CBA > 20 kt in Years 12-15
#' @export
CBAmin_MP <- function(MSEobj) {
  if(!inherits(MSEobj, 'MSE')) stop('This PM method is designed for objects of class `MSE`')
  PM <- .CBA(MSEobj, Yrs = c(12, 15))
  PM@Name <- PM@Caption <- "Probability CBA > 20 kt in Years 12-15"
  PM
}
class(CBAmin_MP) <- "PM"

#' @name PM
#' @details
#' CBAmin_LP - Probability CBA > 20 kt in Years 24-36
#' @export
CBAmin_LP <- function(MSEobj) {
  if(!inherits(MSEobj, 'MSE')) stop('This PM method is designed for objects of class `MSE`')
  PM <- .CBA(MSEobj, Yrs = c(24, 36))
  PM@Name <- PM@Caption <- "Probability CBA > 20 kt in Years 24-36"
  PM
}
class(CBAmin_LP) <- "PM"


#' @name PM
#' @details
#' CBAprom_CP - Mean CBA in Years 1-5
#' @export
CBAprom_CP <- function(MSEobj) {
  if(!inherits(MSEobj, 'MSE')) stop('This PM method is designed for objects of class `MSE`')
  PM <- .CBA(MSEobj, Yrs = c(1, 5), type = "mean")
  PM@Name <- PM@Caption <- "Mean CBA in Years 1-5"
  PM
}
class(CBAprom_CP) <- "PM"



#' @name PM
#' @details
#' CBAprom_MP - Mean CBA in Years 12-15
#' @export
CBAprom_MP <- function(MSEobj) {
  if(!inherits(MSEobj, 'MSE')) stop('This PM method is designed for objects of class `MSE`')
  PM <- .CBA(MSEobj, Yrs = c(12, 15), type = "mean")
  PM@Name <- PM@Caption <- "Mean CBA in Years 12-15"
  PM
}
class(CBAprom_MP) <- "PM"


#' @name PM
#' @details
#' CBAprom_LP - Mean CBA in Years 24-36
#' @export
CBAprom_LP <- function(MSEobj) {
  if(!inherits(MSEobj, 'MSE')) stop('This PM method is designed for objects of class `MSE`')
  PM <- .CBA(MSEobj, Yrs = c(24, 36), type = "mean")
  PM@Name <- PM@Caption <- "Mean CBA in Years 24-36"
  PM
}
class(CBAprom_LP) <- "PM"


#' @name PM
#' @details
#' CBAv_CP - Average variability in CBA in Years 1-5
#' @export
CBAv_CP <- function(MSEobj) {
  if(!inherits(MSEobj, 'MSE')) stop('This PM method is designed for objects of class `MSE`')
  PM <- .CBA(MSEobj, Yrs = c(1, 5), type = "var")
  PM@Name <- PM@Caption <- "Average variability in CBA in Years 1-5"
  PM
}
class(CBAv_CP) <- "PM"

#' @name PM
#' @details
#' CBAv_MP - Average variability in CBA in Years 12-15
#' @export
CBAv_MP <- function(MSEobj) {
  if(!inherits(MSEobj, 'MSE')) stop('This PM method is designed for objects of class `MSE`')
  PM <- .CBA(MSEobj, Yrs = c(12, 15), type = "var")
  PM@Name <- PM@Caption <- "Average variability in CBA in Years 12-15"
  PM
}
class(CBAv_MP) <- "PM"


#' @name PM
#' @details
#' CBAv_LP - Average variability in CBA in Years 24-36
#' @export
CBAv_LP <- function(MSEobj) {
  if(!inherits(MSEobj, 'MSE')) stop('This PM method is designed for objects of class `MSE`')
  PM <- .CBA(MSEobj, Yrs = c(24, 36), type = "var")
  PM@Name <- PM@Caption <- "Average variability in CBA in Years 24-36"
  PM
}
class(CBAv_LP) <- "PM"




# Copied code from https://github.com/pbs-assess/ggmse/blob/master/R/figures-tigure.R#L95
make_gg_df <- function(.probs_dat, .name, mp_order = NULL, do_approx = TRUE, sort_by = "decreasing",
                       relative_max = FALSE, scale_0_1 = FALSE) {
  df <- .probs_dat

  if (is.null(mp_order)) {
    if (sort_by == "decreasing") {
      df$MP <- factor(df$MP, levels = df$MP[do.call(order, df[-1])])
    } else if (sort_by == "increasing") {
      df$MP <- factor(df$MP, levels = df$MP[rev(do.call(order, df[-1]))])
    } else {
      stop("sort_by must be either 'increasing' or 'decreasing'",
           call. = FALSE
      )
    }
  } else {
    df$MP <- factor(df$MP, levels = mp_order)
  }

  df <- reshape2::melt(
    df,
    id.vars = "MP",
    variable.name = "type",
    value.name = "value"
  )

  df$txt <- vapply(df$value, function(x) {
    digits <- ifelse(x <= 1, 2, 1)
    format(round(x, digits), big.mark = ",", decimal.mark = ".", nsmall = digits)
  }, character(1L))

  if (do_approx) {
    OutDec <- options()$OutDec # Decimal is a point or comma?
    df$txt <- gsub(paste0("1\\", OutDec, "00"), paste0(">0", OutDec , "99"), df$txt)
    df$txt <- gsub(paste0("0\\", OutDec, "00"), paste0("<0", OutDec , "01"), df$txt)
  }

  if (relative_max) {
    df <- group_by(df, type) %>%
      mutate(value = value / max(value)) %>%
      ungroup()
  }
  if (scale_0_1) {
    df <- group_by(df, type) %>%
      mutate(value = (value - min(value)) / (max(value) - min(value))) %>%
      ungroup()
  }
  df$MP <- as.factor(df$MP)
  df$OM <- .name
  return(df)
}

#' Performance metric table
#'
#' Create coloured performance metric table in ggplot2 for a set of operating models. Code
#' taken from \url{https://github.com/pbs-assess/ggmse/blob/master/R/figures-tigure.R}
#'
#' @param probs_dat List of data frames containing performance measures
#' @param names Character vector of operating model names corresponding to `probs_dat`
#' @param ncol Integer, optional number of columns.
#' @param relative_max Make the plot have each column use a relative maximum. If
#'   `scale_0_1` is used, this will be ignored
#' @param scale_0_1 Scale each column from 0 to 1, so that the colours in each
#'   column are fully represented
#' @param sort_by show values in decreasing or increasing format
#' @param mp_order Optional hardcoded MP order
#' @param satisficed An optional named numeric vector. The names correspond to
#'   the performance metrics and the values correspond to the values above which
#'   (`>`) the cells will be outlined as "satisficed".
#' @param alpha Transparency of underlying colour.
#' @param do_approx Logical. If `TRUE`, values greater than 0.99 are replaced with ">0.99" and less than 0.01 are replaced
#' "<0.01".
#'
#' @importFrom purrr map_df
#' @return A ggplot2 object
#' @export
plot_table <- function(
    probs_dat,
    names,
    ncol = NULL,
    relative_max = FALSE,
    scale_0_1 = FALSE,
    sort_by = "decreasing",
    mp_order = NULL,
    satisficed = NULL,
    alpha = 0.6,
    do_approx = FALSE
  ) {

  if (!is.list(probs_dat)) stop("probs_dat should be a list of data frames")

  df <- Map(make_gg_df, .probs_dat = probs_dat, .name = names,
            MoreArgs = list(mp_order = mp_order, do_approx = do_approx, sort_by = sort_by,
                            relative_max = relative_max, scale_0_1 = scale_0_1)) %>%
    dplyr::bind_rows()

  padding <- 0.52

  g <- ggplot(df, aes(type, MP)) +
    geom_tile(aes(fill = value), alpha = alpha, color = "white") +
    geom_text(aes(x = type, label = txt), size = ggplot2::rel(3)) +
    scale_fill_gradient2(limits = c(0, 1), midpoint = 0.5) +
    #scale_fill_viridis_c(limits = c(0, 1), begin = 0.15, end = 1, alpha = alpha, option = "D", direction = 1) +
    facet_wrap(vars(OM), scales = "free_x", ncol = ncol) +
    guides(fill = "none") +
    labs(x = NULL, y = NULL) +
    coord_cartesian(
      expand = FALSE,
      xlim = range(as.numeric(df$type)) + c(-padding, padding),
      ylim = range(as.numeric(df$MP)) + c(-padding - 0.01, padding + 0.01)
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(color = "grey10"),
      strip.placement = "outside",
      strip.background = element_blank()
    ) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(labels = levels(df$MP)) +

  if (!is.null(satisficed)) {
    h <- purrr::map_df(
      seq_along(satisficed),
      ~ dplyr::filter(df, value > satisficed[[.x]] & type == names(satisficed)[.x])
    )
    g <- g + geom_tile(data = h, color = "grey30", lwd = 0.45, fill = NA)
  }

  g
}


#' Plot tradeoff figure
#'
#' Create figure to compare management procedures across 2 performance metrics. Code is copied from
#' \url{https://github.com/pbs-assess/ggmse/blob/master/R/figures-bivariate.R#L26}
#'
#' @param pm_df_list A named list of performance metric data frames. The names will be used as the plot labels.
#' @param xvar The performance metric for the x axis (as character).
#' @param yvar The performance metric for the y axis (as character).
#' @param custom_pal An optional custom color palette. Should be a named character vector.
#' @param mp An optional character vector of MPs to include. By default includes all.
#' @param nudge_x How far to nudge the labels left/right from the x-value for \link[ggplot2]{geom_text}
#' @param nudge_y How far to nudge the labels up/down from the y-value for \link[ggplot2]{geom_text}
#' @param mp_ref Character vector of reference management procedures
#' @return A ggplot2 object
#' @import reshape2
#' @export
plot_tradeoff <- function(
    pm_df_list, xvar, yvar, custom_pal = NULL,
    mp = NULL,
    mp_ref = NULL,
    nudge_x = 0, nudge_y = 0.05
  ) {

  df <- lapply(names(pm_df_list), function(x) {
    pm_df_list[[x]] %>%
      mutate(scenario = x, MP_label = 1:nrow(.))
  }) %>%
    bind_rows()

  if (!is.null(mp)) {
    df <- dplyr::filter(df, MP %in% mp) %>% mutate(MP = factor(MP, levels = mp))
  }
  df_long <- reshape2::melt(
    df,
    id.vars = c("MP", "scenario", "MP_label"),
    value.name = "prob",
    variable.name = "pm"
  )
  df_wide <- df_long %>%
    reshape2::dcast(MP + MP_label + scenario ~ pm, value.var = "prob") %>%
    dplyr::mutate(`Reference` = MP %in% mp_ref)

  xmin <- pull(df_wide, !!xvar) %>% min()
  ymin <- pull(df_wide, !!yvar) %>% min()
  #xvar <- paste0("`", xvar, "`")
  #yvar <- paste0("`", yvar, "`")

  n_mp <- length(unique(df_wide$MP))
  ref_or_not <- dplyr::select(df_wide, .data$MP, .data$Reference) %>% dplyr::distinct()
  #mp_shapes <- vector(mode = "numeric", length = n_mp)
  mp_shapes <- ifelse(ref_or_not$Reference, 1, 16) %>%
    structure(names = as.character(ref_or_not$MP))
  mp_label <- filter(df_wide, scenario == names(pm_df_list)[[1]]) %>% pull("MP_label")

  g <- ggplot(
    df_wide,
    aes(.data[[xvar]], .data[[yvar]], colour = MP, shape = MP)
  ) +
    geom_text(aes(label = .data$MP_label), nudge_x = nudge_x, nudge_y = nudge_y) +
    geom_point(show.legend = FALSE) +
    guides(colour = guide_legend(override.aes = list(label = mp_label))) +
    facet_wrap(vars(scenario), nrow = 2) +
    scale_shape_manual(values = mp_shapes) +
    theme(
      strip.background = element_blank(),
      panel.grid.major.y = element_line(colour = "grey85"),
      panel.grid.major.x = element_line(colour = "grey85")
    ) +
    labs(colour = "PM", fill = "PM")

  if (!is.null(custom_pal)) {
    g <- g + scale_color_manual(values = custom_pal)
  }

  g
}

#' Calculate performance metric
#'
#' Wrapper function that conveniently calculates the values for a set of performance metrics
#'
#' @param .mse MSE object
#' @param PMs Character vector of the performance metric functions
#' @param all_sims Logical, whether to report the value for individual simulations (TRUE), or averaged over all simulations (FALSE)
#' @return If `all_sims = TRUE` a matrix, otherwise a data frame
#' @importFrom methods slot
#' @export
calculate_PM <- function(.mse, PMs, all_sims = FALSE) {
  if (all_sims) {
    PMobj <- sapply(1:length(PMs), function(i) {
      x <- get(PMs[i])(.mse)
      slot(x, 'Prob')
    }, simplify = "array") %>%
      structure(dimnames = list(Sim = 1:.mse@nsim, MP = .mse@MPs, PM = PMs))
  } else {
    PMobj <- sapply(1:length(PMs), function(i) {
      x <- get(PMs[i])(.mse)
      slot(x, 'Mean')
    }) %>%
      structure(dimnames = list(MP = .mse@MPs, PM = PMs)) %>%
      as.data.frame() %>%
      mutate(MP = .mse@MPs)
  }

  return(PMobj)
}

update_mean <- function(old_mean, xnext, n) (n-1)/n * old_mean + xnext/n

cumulative_mean <- function(PMobj) {
  x <- array(NA_real_, dim(PMobj)) %>% structure(dimnames = dimnames(PMobj))
  x[1, , ] <- PMobj[1, , ]
  for(i in 2:dim(x)[1]) {
    x[i, , ] <- update_mean(x[i-1, , ], PMobj[i, , ], i)
  }
  return(x)
}

#' @importFrom rlang .env
plot_cumulative_PM <- function(x, MP = dimnames(x)$MP, PM = "NZR_MP") {
  df <- reshape2::melt(x) %>%
    filter(MP %in% .env$MP, PM %in% .env$PM, !is.na(value))

  g <- ggplot(df, aes(Sim, value, colour = MP)) +
    geom_line() +
    facet_wrap(vars(PM), scales = "free_y") +
    theme(
      strip.placement = "outside",
      strip.background = element_blank(),
      legend.position = "bottom"
    ) +
    labs(x = "Simulación de cada MO", y = "Valor acumulado", colour = "PM")
  g
}
