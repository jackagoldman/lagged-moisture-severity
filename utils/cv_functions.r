compare_models_cv <- function(model1, model2, data, response_var, cv_type = "repeated_holdout", k = 5, n_reps = 10, seed = 1) {
  if (cv_type == "repeated_holdout") {
    rmse1_vec <- numeric(n_reps)
    rmse2_vec <- numeric(n_reps)
    
    for (i in 1:n_reps) {
      set.seed(seed + i - 1)
      id <- sample.int(nrow(data), floor(0.2 * nrow(data)))
      train <- data[-id, ]
      test <- data[id, ]
      
      m1_cv <- gam(formula(model1), data = train, family = family(model1), method = "REML")
      m2_cv <- gam(formula(model2), data = train, family = family(model2), method = "REML")
      
      p1 <- predict(m1_cv, newdata = test)
      p2 <- predict(m2_cv, newdata = test)
      
      rmse1_vec[i] <- sqrt(mean((test[[response_var]] - p1)^2, na.rm = TRUE))
      rmse2_vec[i] <- sqrt(mean((test[[response_var]] - p2)^2, na.rm = TRUE))
    }
    
    c(RMSE_model1 = mean(rmse1_vec), RMSE_model2 = mean(rmse2_vec))
    
  } else if (cv_type == "kfold") {
    set.seed(seed)
    n <- nrow(data)
    indices <- sample(1:n)
    fold_size <- floor(n / k)
    rmse1_vec <- numeric(k)
    rmse2_vec <- numeric(k)
    
    for (i in 1:k) {
      start_idx <- (i - 1) * fold_size + 1
      end_idx <- if (i < k) i * fold_size else n
      test_idx <- indices[start_idx:end_idx]
      train <- data[-test_idx, ]
      test <- data[test_idx, ]
      
      m1_cv <- gam(formula(model1), data = train, family = family(model1), method = "REML")
      m2_cv <- gam(formula(model2), data = train, family = family(model2), method = "REML")
      
      p1 <- predict(m1_cv, newdata = test)
      p2 <- predict(m2_cv, newdata = test)
      
      rmse1_vec[i] <- sqrt(mean((test[[response_var]] - p1)^2, na.rm = TRUE))
      rmse2_vec[i] <- sqrt(mean((test[[response_var]] - p2)^2, na.rm = TRUE))
    }
    
    c(RMSE_model1 = mean(rmse1_vec), RMSE_model2 = mean(rmse2_vec))
    
  } else if (cv_type == "loo") {
    n <- nrow(data)
    rmse1_vec <- numeric(n)
    rmse2_vec <- numeric(n)
    
    for (i in 1:n) {
      train <- data[-i, ]
      test <- data[i, , drop = FALSE]
      
      m1_cv <- gam(formula(model1), data = train, family = family(model1), method = "REML")
      m2_cv <- gam(formula(model2), data = train, family = family(model2), method = "REML")
      
      p1 <- predict(m1_cv, newdata = test)
      p2 <- predict(m2_cv, newdata = test)
      
      rmse1_vec[i] <- sqrt((test[[response_var]] - p1)^2)
      rmse2_vec[i] <- sqrt((test[[response_var]] - p2)^2)
    }
    
    c(RMSE_model1 = mean(rmse1_vec), RMSE_model2 = mean(rmse2_vec))
  } else {
    stop("Invalid cv_type. Choose 'repeated_holdout', 'kfold', or 'loo'.")
  }
}


