
mask_too_far <- function(model, x, y, n = 100, dist = 0.05) {
  m <- model.frame(model)

  # check presence
  miss <- setdiff(c(x, y), names(m))
  if (length(miss) > 0) {
    stop("Variable(s) not found in model frame: ", paste(miss, collapse = ", "))
  }

  # helper to coerce to numeric vector, handling matrices/list-columns
  as_num_vec <- function(z) {
    if (is.matrix(z) || is.data.frame(z)) return(as.numeric(z))
    if (is.list(z)) return(unlist(lapply(z, function(el) if (length(el) == 1) as.numeric(el) else as.numeric(el))))
    as.numeric(z)
  }

  x_vals <- as_num_vec(m[[x]])
  y_vals <- as_num_vec(m[[y]])

  # keep only finite observed values
  ok <- is.finite(x_vals) & is.finite(y_vals)
  if (!any(ok)) stop("No finite observed pairs found for variables '", x, "' and '", y, "'. Inspect your data.")
  x_obs <- x_vals[ok]
  y_obs <- y_vals[ok]

  # sequences for grid (now guaranteed finite)
  x_seq <- seq(min(x_obs), max(x_obs), length.out = n)
  y_seq <- seq(min(y_obs), max(y_obs), length.out = n)

  grd <- expand.grid(x = x_seq, y = y_seq)
  names(grd) <- c(x, y)

  # compute mask: TRUE where grid point is too far from any observed pair
  tf_mask <- mgcv::exclude.too.far(grd[[x]], grd[[y]], x_obs, y_obs, dist = dist)

  # return logical vector aligned to grid rows (length = n*n)
  tf_mask
}
#' get lagged smooth estimates from a GAM model
#' @param model, a fitted GAM model
#' @param lag_var, name of the lag variable in the model
#' @param smooth_var, name of the smooth variable in the model
#' @param weights, name of the weights variable in the model
#' @param n, number of points to evaluate the smooth over
#' @return a list with two dataframes: smooth_df (all estimates) and mask (significant estimates)

get_lag_estimates <- function(model, lag_var = "L", smooth_var = "vpd", weights = "weights_btl", n = 100) {
  
  if (!is.character(lag_var)) stop("lag_var must be a string, the name of the lag variable in the model")
  if (!is.character(smooth_var)) stop("smooth_var must be a string, the name of the smooth variable in the model")
  if (!is.character(weights)) stop("weights must be a string, the name of the weights variable in the model")
  
  # Create the smooth specification
  smooth = paste0("te(", smooth_var, ",", lag_var,"):",  weights)


  tf <- mask_too_far(model = model, x = lag_var, y = smooth_var, n = 100, dist = 0.1)

 
  # get intercept from gam object
  intercept <- model$coefficients[[1]]


  smooth_df <- gratia::smooth_estimates(
    object = model,
    smooth = smooth,
    partial_match = TRUE,
    unconditional = TRUE,
    n = 100
  ) |>
    gratia::add_confint(type = "confidence") |> 
# Flag significant locations where CI excludes 0 (null effect)
    tibble::add_column(tf) |> 
    #dplyr::mutate(across(c(.estimate, .lower_ci, .upper_ci), ~ if_else(tf, NA_real_, .))) |>
    dplyr::mutate(sig = ifelse(.lower_ci > 0 | .upper_ci < 0, TRUE, FALSE)) 

return(smooth_df)
  }


#' plot heatmap of lagged smooth effects from a GAM model
#' @param smooth_df, dataframe of smooth estimates with confidence intervals
#' @param mask, dataframe of significant smooth estimates (CI excludes 0)



plot_lag_smooth <- function(model, 
  smooth_df, 
  smooth_var, 
  lag_var, 
  weights) {
  ylabel <- switch(smooth_var,
                   vpd = "Vapour Pressure Deficit (z-score)",
                   cmi = "Climate Moisture Index (z-score)",
                   smooth_var)
  response <- as.character(formula(model)[[2]])
  response_lab <- if (response == "median") "Median\nBurn Severity" else if (response == "pct_95") "Extreme\nBurn Severity" else response
  #if response = "median" set breaks  = 1:30, lmits = c(1,30) but if response = "pct_90" set breaks = 1:5, limits = c(1,5)
  breaks <- if (smooth_var == "vpd") seq(1,30,5) else if (smooth_var == "cmi_wy") 1:5 else 1:4
  limits <- if (smooth_var == "vpd") c(1,30) else if (smooth_var == "cmi_wy") c(1,5) else c(1,4)
  # base heatmap mapping columns by name
  xlab<- if (smooth_var == "cmi_wy") "Lag (years before fire start)" else if (smooth_var == "vpd") "Lag (days before fire start)" else "Lag (weeks before fire start)"
  ylab <- if (smooth_var == "cmi_wy") "CMI" else if (smooth_var == "vpd") "VPD"
  p <- ggplot(smooth_df, aes(x = .data[[lag_var]], y = .data[[smooth_var]],z =.estimate, fill = .estimate)) +
    geom_raster() +
    geom_contour(color = "grey50", size = 0.4) +
      geom_raster(data = filter(smooth_df, tf== TRUE), fill = "white", alpha = 0.4) +
      scale_fill_viridis_c(
        name = response_lab,
        option = "A",
        guide = guide_colorbar(
          direction = "vertical",
          title.position = "top",
          title.hjust = 0.5
        )
      ) +
      theme(legend.position = "bottom") +
      scale_x_continuous(
        name = xlab,
        breaks = breaks,
        expand = c(0, 0)
      ) +
    scale_y_continuous(name = ylab, expand = c(0, 0)) +
    theme_bw() +
    theme(
      panel.spacing = grid::unit(0, "pt"),
      plot.margin = grid::unit(c(0, 0, 0, 0), "pt"),
      aspect.ratio = 0.5
    )

  # outline significant cells (sig == TRUE)
  mask <- dplyr::filter(smooth_df, sig == TRUE)

  if (nrow(mask) > 0) {
     p <- p +
      geom_tile(
        data = mask,
        aes(x = .data[[lag_var]], y = .data[[smooth_var]]),
        color = "black",
        size = 0.75,
        linejoin = "round",
        fill = NA,
        inherit.aes = FALSE
      ) +  
      geom_tile(data = mask, aes(x = .data[[lag_var]], y = .data[[smooth_var]], fill = .estimate), inherit.aes = FALSE) + geom_contour(color = "grey50", size = 0.4) 

  }
  
  # apply x-limits via coord_cartesian (does not drop data)
  if (!is.null(limits)) {
    p <- p + coord_cartesian(xlim = limits, expand = FALSE)
  } else {
    p <- p + coord_cartesian(expand = FALSE)
  }

  p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "grey") 
  p
}

###


#' place label inside ggplot at specified position
#' @param plot, ggplot object to modify
#' @param label, text label to place
#' @param x_pos, x position for label (default -Inf for left side)
#' @param y_pos, y position for label (default 100)
#' @param hjust, horizontal justification (default -0.5 to nudge left
#' 
place_label_inside <- function(plot, label, x_pos = -Inf, y_pos = 100, hjust = -0.5, vjust = 1, size = 5, padding = unit(5, "pt")) {
  plot +
    annotate("text",
             label = label, x = x_pos, y = y_pos,
             hjust = hjust, vjust = vjust, size = size,
             gp = gpar(fontsize = size) # Use grid::gpar for more control if needed
    ) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "pt")) # Optional: remove margins
}


#' get min max estimates
#' 
#' 
get_steepest <- function(mod.estimates){
  max_est <- mod.estimates |> 
    dplyr::filter(tf == FALSE) |>
    dplyr::arrange(dplyr::desc(.estimate)) |>
    dplyr::slice_head(n = 5)

  min_est <- mod.estimates |> 
    dplyr::filter(tf == FALSE) |>
    dplyr::arrange(.estimate) |>
    dplyr::slice_head(n = 5)


  sig_est <- mod.estimates |> 
    dplyr::filter(sig == TRUE)

  return(list(max = max_est, min = min_est, sig = sig_est))
}
