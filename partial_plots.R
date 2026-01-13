library(mgcv)
library(gratia)
library(tidyverse)
library(ggplot2)

# Load model and data (adjust paths as needed)
m_cmi_median <- readRDS("path/to/m_cmi_median.rds")
moisture_data <- readRDS("path/to/moisture_data.rds")

#number of observations per ecoregin lvl
moisture_data %>%
  group_by(ecoregion_lvl) %>%
  summarise(n = n())

m <- gam(median~ 
       te(cmi_wy, LL, k = c(30,5))+# cmi with 5 year lag
        s(pct_conf, k = 10) +  
         s(mean_evi, k = 10) + 
         s(fwi_90th, k = 10) +
  s(ecoregion_lvl, bs = "re") ,
      data = moisture_data,
      family = scat(link = "identity"),
      method = "REML")
summary(m)
gam.check(m)

draw(m)

# Get coefficients for the te smooth
smooth_obj <- m$smooth[[1]]
coefs <- smooth_obj$first.para : smooth_obj$last.para
beta <- coef(m)



# Assuming cmi_quartiles is defined as above
cmi_quartiles <- quantile(moisture_data$cmi_wy, probs = c(0.05, 0.25, 0.5, 0.75), na.rm = TRUE)

# Sequence for lags (LL)
ll_seq <- seq(min(moisture_data$LL, na.rm = TRUE), max(moisture_data$LL, na.rm = TRUE), length.out = 100)

# Fix other covariates at medians
pct_conf_med <- median(moisture_data$pct_conf, na.rm = TRUE)
mean_evi_med <- median(moisture_data$mean_evi, na.rm = TRUE)
fwi_90th_med <- median(moisture_data$fwi_90th, na.rm = TRUE)
ecoregion_lvl <- 'BTL'  # Representative ecoregion

# Get coefficients for the te smooth
smooth_obj <- m$smooth[[1]]
coefs <- smooth_obj$first.para : smooth_obj$last.para
beta <- coef(m)

# Loop over quartiles and compute partial effects
plot_dat_list <- lapply(names(cmi_quartiles), function(q) {
  cmi_fixed <- cmi_quartiles[q]
  newdat <- data.frame(cmi_wy = cmi_fixed, LL = ll_seq, pct_conf = pct_conf_med, 
                       mean_evi = mean_evi_med, fwi_90th = fwi_90th_med, 
                       ecoregion_lvl = ecoregion_lvl)
  newXp <- predict(m, type = 'lpmatrix', newdata = newdat)
  newXp_adj <- matrix(0, nrow = nrow(newXp), ncol = ncol(newXp))
  newXp_adj[, coefs] <- newXp[, coefs]
  partial_effect <- newXp_adj[, coefs] %*% beta[coefs]
  data.frame(LL = ll_seq, partial_effect = as.vector(partial_effect), quartile = q)
})

plot_dat <- do.call(rbind, plot_dat_list)

plot_dat$quartile <- factor(plot_dat$quartile, levels = c("5%", "25%", "50%", "75%"))

# Plot (unchanged, but legend will now be in this order)
ggplot(plot_dat, aes(x = LL, y = partial_effect, color = quartile)) +
  geom_line(linewidth = 1) +
  labs(x = "Lag (months)", y = "Partial Effect of CMI on Median Severity", 
       title = "Partial Effects at CMI Quartiles", color = "CMI Quartile") +
  scale_x_reverse(breaks = 1:5, labels = paste0("-", 1:5)) +
  scale_color_viridis_d() +  # Viridis colors for lines
  scale_fill_viridis_d() +   # Viridis colors for ribbons
  theme_bw()

#how many cmi_wy values below 0
sum(moisture_data$cmi_wy < 0, na.rm = TRUE)


m1<- gam(pct_95~ 
       te(cmi_wy, LL, k = c(30,5))+# cmi with 5 year lag
        s(pct_conf, k = 10) +  
         s(mean_evi, k = 10) + 
         s(fwi_90th, k = 10) +
  s(ecoregion_lvl, bs = "re") ,
      data = moisture_data,
      family = scat(link = "identity"),
      method = "REML")
summary(m1)
gam.check(m1)



summary(m)


# Assuming cmi_quartiles is defined as above
cmi_quartiles <- quantile(moisture_data$cmi_wy, probs = c(0.05, 0.25, 0.5, 0.75), na.rm = TRUE)

# Sequence for lags (LL)
ll_seq <- seq(min(moisture_data$LL, na.rm = TRUE), max(moisture_data$LL, na.rm = TRUE), length.out = 100)

# Fix other covariates at medians
pct_conf_med <- median(moisture_data$pct_conf, na.rm = TRUE)
mean_evi_med <- median(moisture_data$mean_evi, na.rm = TRUE)
fwi_90th_med <- median(moisture_data$fwi_90th, na.rm = TRUE)
ecoregion_lvl <- 'BTL'  # Representative ecoregion

# Get coefficients for the te smooth
smooth_obj <- m1$smooth[[1]]
coefs <- smooth_obj$first.para : smooth_obj$last.para
beta <- coef(m1)

# Loop over quartiles and compute partial effects
plot_dat_list <- lapply(names(cmi_quartiles), function(q) {
  cmi_fixed <- cmi_quartiles[q]
  newdat <- data.frame(cmi_wy = cmi_fixed, LL = ll_seq, pct_conf = pct_conf_med, 
                       mean_evi = mean_evi_med, fwi_90th = fwi_90th_med, 
                       ecoregion_lvl = ecoregion_lvl)
  newXp <- predict(m1, type = 'lpmatrix', newdata = newdat)
  newXp_adj <- matrix(0, nrow = nrow(newXp), ncol = ncol(newXp))
  newXp_adj[, coefs] <- newXp[, coefs]
  partial_effect <- newXp_adj[, coefs] %*% beta[coefs]
  data.frame(LL = ll_seq, partial_effect = as.vector(partial_effect), quartile = q)
})

plot_dat <- do.call(rbind, plot_dat_list)

plot_dat$quartile <- factor(plot_dat$quartile, levels = c("5%", "25%", "50%", "75%"))

# Plot (unchanged, but legend will now be in this order)
ggplot(plot_dat, aes(x = LL, y = partial_effect, color = quartile)) +
  geom_line(linewidth = 1) +
  labs(x = "Lag (months)", y = "Partial Effect of CMI on Median Severity", 
       title = "Partial Effects at CMI Quartiles", color = "CMI Quartile") +
  scale_x_reverse(breaks = 1:5, labels = paste0("-", 1:5)) +
  theme_minimal()

### VPD
# 
m2 <- gam(median~ 
       te(vpd, L, k = c(30,10))+# cmi with 5 year lag
        s(pct_conf, k = 10) +  
         s(mean_evi, k = 10) + 
         s(fwi_90th, k = 10) +
  s(ecoregion_lvl, bs = "re") ,
      data = moisture_data,
      family = scat(link = "identity"),
      method = "REML")

summary(m2)

vpd_quartiles <- quantile(moisture_data$vpd, probs = c(0.05, 0.25, 0.5, 0.75), na.rm = TRUE)

# Get coefficients for the te smooth
smooth_obj <- m2$smooth[[1]]
coefs <- smooth_obj$first.para : smooth_obj$last.para
beta <- coef(m2)

# Loop over quartiles and compute partial effects
plot_dat_list <- lapply(names(vpd_quartiles), function(q) {
  vpd_fixed <- vpd_quartiles[q]
  newdat <- data.frame(vpd = vpd_fixed, L = l_seq, pct_conf = pct_conf_med, 
                       mean_evi = mean_evi_med, fwi_90th = fwi_90th_med, 
                       ecoregion_lvl = ecoregion_lvl)
  newXp <- predict(m2, type = 'lpmatrix', newdata = newdat)
  newXp_adj <- matrix(0, nrow = nrow(newXp), ncol = ncol(newXp))
  newXp_adj[, coefs] <- newXp[, coefs]
  partial_effect <- newXp_adj[, coefs] %*% beta[coefs]
  data.frame(LL = ll_seq, partial_effect = as.vector(partial_effect), quartile = q)
})

plot_dat <- do.call(rbind, plot_dat_list)

plot_dat$quartile <- factor(plot_dat$quartile, levels = c("5%", "25%", "50%", "75%"))

# Plot (unchanged, but legend will now be in this order)
ggplot(plot_dat, aes(x = LL, y = partial_effect, color = quartile)) +
  geom_line(linewidth = 1) +
  labs(x = "Lag (months)", y = "Partial Effect of CMI on Median Severity", 
       title = "Partial Effects at CMI Quartiles", color = "CMI Quartile") +
  scale_x_reverse(breaks = 1:5, labels = paste0("-", 1:5)) +
  theme_minimal()

#how many cmi_wy values below 0
m3 <- gam(pct_95~ 
       te(vpd, L, k = c(30,10))+# cmi with 5 year lag
        s(pct_conf, k = 10) +  
         s(mean_evi, k = 10) + 
         s(fwi_90th, k = 10) +
  s(ecoregion_lvl, bs = "re") ,
      data = moisture_data,
      family = scat(link = "identity"),
      method = "REML")

summary(m3)

