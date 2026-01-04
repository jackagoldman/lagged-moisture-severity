#' check for adequate knots
#' 
#' 
#' 
#' 

check_res_edf <- function(model,  k1, k2) {
  res <- residuals(model)
  data <- as.data.frame(model.frame(model))
  sm_terms <- gratia::smooths(model)
  #weights
  weights_1 <- sub(".*:", "",  sm_terms[[1]])
  weights_2 <- sub(".*:", "",  sm_terms[[2]])
  weights_3 <- sub(".*:", "",  sm_terms[[3]])
  # argument from tensor smooth
  te1 <- sub(":.*$", "",sm_terms[[1]] )
  te_arg <- sub("te\\((.*)\\)", "\\1", te1)

  # if k1 and k2 are not provided, extract from model
  if (missing(k1) | missing(k2)) {
    k1 <- model$smooth[[1]]$margin[[1]]$bs.dim
    k2 <- model$smooth[[1]]$margin[[2]]$bs.dim
      } else {
        k1 <- k1
        k2 <- k2
      }
  # construct formula using the actual column names
  form_str1 <- paste0("res ~ te(", te_arg,  ", by = ", weights_1, ", k = c(", k1,",", k2, "), bs = c('tp','cr'))")
  form_str2 <- paste0("res ~ te(", te_arg,  ", by = ", weights_2, ", k = c(", k1,",", k2, "), bs = c('tp','cr'))")
  form_str3 <- paste0("res ~ te(", te_arg,  ", by = ", weights_3, ", k = c(", k1,",", k2, "), bs = c('tp','cr'))")

  # fit models
  gmod.1 <- mgcv::gam(as.formula(form_str1), data = data, gamma = 1.4)
  gmod.2 <- mgcv::gam(as.formula(form_str2), data = data, gamma = 1.4)
  gmod.3 <- mgcv::gam(as.formula(form_str3), data = data, gamma = 1.4)
    
  
  cat("Check edf for weights:", weights_1, "\n")
  print(gmod.1 %>%
    gratia::edf())
  cat("\nCheck edf for weights:", weights_2, "\n")
  print(gmod.2 %>%
    gratia::edf())
  cat("\nCheck edf for weights:", weights_3, "\n")
  print(gmod.3 %>%
    gratia::edf())

}

