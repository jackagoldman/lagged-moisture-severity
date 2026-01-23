
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

get_lag_estimates <- function(model, lag_var = "L", smooth_var = "vpd", n = 100) {
  
  if (!is.character(lag_var)) stop("lag_var must be a string, the name of the lag variable in the model")
  if (!is.character(smooth_var)) stop("smooth_var must be a string, the name of the smooth variable in the model")
  
  # Create the smooth specification
  smooth = paste0("te(", smooth_var, ",", lag_var,")")


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
    dplyr::mutate(sig = ifelse(.lower_ci > 0 | .upper_ci < 0 , TRUE, FALSE)) 

return(smooth_df)
  }


#' plot heatmap of lagged smooth effects from a GAM model
#' @param smooth_df, dataframe of smooth estimates with confidence intervals
#' @param mask, dataframe of significant smooth estimates (CI excludes 0)



plot_lag_smooth <- function(model, 
  smooth_df, 
  smooth_var, 
  lag_var) {
  ylabel <- switch(smooth_var,
                   vpd = "Vapour Pressure Deficit",
                   cmi = "Climate Moisture Index",
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




#' draws function 
#' @param model, a fitted GAM model
#' @param smooth_id, index of the tensor smooth in the model
#' @param q, quantile to fix the other variable at
#' @param moisture_data, dataframe with moisture data to get quantiles and medians
#' @param seed, random seed for reproducibility
#' @return a list with curves dataframe, median dataframe, smooth variable name, and plot
#' 
#' Example usage: 
#' test <- one_d_slice_draws(model = m_median, smooth_id = 1, q = "5%", moisture_data = moisture_data)
#' p <- test$plot # get plot
#' p # display plot
#' curves_df <- test$curves # get curves dataframe
#' median_df <- test$median # get median dataframe

one_d_slice_draws <- function(model, smooth_id = 1, q = "5%" , moisture_data, seed = 42) {
    # choose model (m or m1)
    model <- model      # replace with m or m1 as appropriate
    smooth_id <- smooth_id     # index of the tensor smooth (check by names(model$smooth))
  
  get_smooth_cols <- function(mod, smooth_id = 1) {
    s <- mod$smooth[[smooth_id]] 
    seq(s$first.para, s$last.para)
  }

  

    if (smooth_id == 1){

    # build grid/newdata for a particular CMI quartile (example: 50%)
    quartiles <- quantile(moisture_data$cmi_wy, probs = c(0.05,0.25,0.5,0.75), na.rm = TRUE)

    # Sequence for lags (LL)
    ll_seq <- seq(min(moisture_data$LL, na.rm = TRUE), max(moisture_data$LL, na.rm = TRUE), length.out = 100)

    # Fix other covariates at medians
    pct_conf_med <- median(moisture_data$pct_conf, na.rm = TRUE)
    dc_90th_med <- median(moisture_data$dc_90th, na.rm = TRUE)
    vpd_med <- median(moisture_data$vpd, na.rm = TRUE)
    l_seq <- min(moisture_data$L, na.rm = TRUE)
    fire_size_class <- "large"
    ecoregion_lvl <- 'LN'  # Representative ecoregion


    newdat <- data.frame(
      cmi_wy = rep(quartiles[q], length(ll_seq)),
      LL = ll_seq,
      pct_conf = pct_conf_med,
      dc_90th = dc_90th_med,
      fire_size_class = fire_size_class,
      vpd = vpd_med,
      L = l_seq,
      ecoregion_lvl = ecoregion_lvl
    )
  }else {
    # build grid/newdata for a particular CMI quartile (example: 50%)
    quartiles <- quantile(moisture_data$vpd, probs = c(0.05,0.25,0.5,0.75), na.rm = TRUE)

    # Sequence for lags (LL)
    l_seq <- seq(min(moisture_data$L, na.rm = TRUE), max(moisture_data$L, na.rm = TRUE), length.out = 100)

    # Fix other covariates at medians
    pct_conf_med <- median(moisture_data$pct_conf, na.rm = TRUE)
    dc_90th_med <- median(moisture_data$dc_90th, na.rm = TRUE)
    cmi_wy_med <- median(moisture_data$cmi_wy, na.rm = TRUE)
    ll_seq <- min(moisture_data$LL, na.rm = TRUE)
    fire_size_class <- "large"
    ecoregion_lvl <- 'LN'  # Representative ecoregion


    newdat <- data.frame(
      vpd = rep(quartiles[q], length(l_seq)),
      LL = ll_seq,
      pct_conf = pct_conf_med,
      dc_90th = dc_90th_med,
      fire_size_class = fire_size_class,
      cmi_wy = cmi_wy_med,
      L = l_seq,
      ecoregion_lvl = ecoregion_lvl
    )
  }

    B <- 200
    beta_hat <- coef(model)
    Vb <- vcov(model)
    df_t <- 1000   # choose degrees of freedom

    # get Xp_partial (one-time)
    Xp <- predict(model, newdata = newdat, type = "lpmatrix")
    cols <- get_smooth_cols(model, smooth_id)
    Xp_part <- matrix(0, nrow = nrow(Xp), ncol = ncol(Xp))
    Xp_part[, cols] <- Xp[, cols, drop = FALSE]

    # simulate B coefficient draws from multivariate-t and compute curves
    set.seed(seed)
    beta_draws <- mvtnorm::rmvt(B, sigma = Vb, df = df_t, delta = beta_hat, type = "shifted")
    # beta_draws: B x p, convert to p x B if needed
    curve_mat <- t(apply(beta_draws, 1, function(b) as.numeric(Xp_part %*% b))) # B x ngrid
    curve_mat <- t(curve_mat)   # ngrid x B


    B_plot <- min(50, ncol(curve_mat))
    cols_to_plot <- sample(ncol(curve_mat), B_plot)

    if (smooth_id == 1){
    smooth_var <- "cmi_wy"
    df_long <- as_tibble(curve_mat[, cols_to_plot]) %>%
      mutate(LL = newdat$LL) %>%
      pivot_longer(-LL, names_to = "sim", values_to = "partial") |> mutate(quartile  = q)

    # median curve across refits
    median_curve <- apply(curve_mat, 1, mean, na.rm = TRUE)
    median_df <- tibble(LL = newdat$LL, fit = median_curve,  quartile = q)
     plot <- ggplot(df_long, aes(x = LL, y = partial, group = sim)) +
        geom_line(alpha = 0.12, color = "black") +
        geom_line(data = median_df, aes(x = LL, y = fit), color = "firebrick", linewidth = 1.1, inherit.aes = FALSE) +
        scale_x_reverse() +
        labs(x = "Lag", y = "Partial effect", title = "Refit subsample curves + median") +
        theme_bw()+ geom_hline(yintercept = 0, linetype = "dashed", color = "grey")
    } else {
      smooth_var <- "vpd"
      df_long <- as_tibble(curve_mat[, cols_to_plot]) %>%
        mutate(L = newdat$L) %>%
        pivot_longer(-L, names_to = "sim", values_to = "partial")|> mutate(quartile  = q)

      # median curve across refits
      median_curve <- apply(curve_mat, 1, mean, na.rm = TRUE)
      median_df <- tibble(L = newdat$L, fit = median_curve, quartile = q)
       plot <- ggplot(df_long, aes(x = L, y = partial, group = sim)) +
          geom_line(alpha = 0.12, color = "black") +
          geom_line(data = median_df, aes(x = L, y = fit), color = "firebrick", linewidth = 1.1, inherit.aes = FALSE) +
          scale_x_reverse() +
          labs(x = "Lag", y = "Partial effect", title = "Refit subsample curves + median") +
          theme_bw() +  geom_hline(yintercept = 0, linetype = "dashed", color = "grey")

    }
 
    return(list(curves = df_long, median = median_df, smooth_var = smooth_var, plot = plot, curve_mat = curve_mat))
}


