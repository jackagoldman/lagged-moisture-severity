make_gam_summary_table <- function(model) {
  require(gt)
  require(dplyr)
  s_table <- summary(model)$s.table
  dev_expl <- summary(model)$dev.expl

  # get model family 
  model_family <- family(model)$family
  
  
  
  s_df <- as.data.frame(s_table)
  s_df<- rownames_to_column(s_df, var = "Term")

  
  #turn off scientific notation
  options(scipen = 999)

  # round all columns but if p-value is <0.001 make <0.001
 

  if (data.frame(model_family)== "gaussian"){
  s_df <- s_df %>%
    mutate(`p-value` = ifelse(`p-value` < 0.001, "<0.001", round(`p-value`, 3))) %>%
    mutate(across(c(edf, Ref.df, F), ~ round(., 3))) } else {
    s_df <- s_df %>%
      mutate(`p-value` = ifelse(`p-value` < 0.001, "<0.001", round(`p-value`, 3))) %>%
      mutate(across(c(edf, Ref.df, Chi.sq), ~ round(., 3)))
  }

  

  # remove ref.df column
  s_df <- s_df |> dplyr::select(-Ref.df)

 

  # rename Chi.sq to Test statistic
  if (data.frame(model_family) == "gaussian"){
  colnames(s_df)[colnames(s_df) == "F"] <- "Test statistic"
  } else {
  colnames(s_df)[colnames(s_df) == "Chi.sq"] <- "x²"
  }
  
  response <- as.character(model$formula)[2]
  combined <- s_df

  # add response column before ecoregion and fill first row with response and rest with -
  combined <- combined %>%
    mutate(Response = response) %>%
    dplyr::select(Response, everything()) %>%
    mutate(Response = ifelse(row_number() == 1, response, "-"))

  # if response = pct_95 change to extreme
  combined$Response <- recode(combined$Response, 
                              "pct_95" = "extreme")

  # In term column if value is pct_conf change to percent conifer, change fwi_90th to FWI , change mean_evi to EVI and s(ecoregion_lvl) to s(ecoregion)
  combined$Term <- recode(combined$Term, 
                          "s(pct_conf)" = "s(percent conifer)", 
                          "s(dc_90th)" = "s(Drought Code)",
                          "s(fyear)" = "s(fire year)",
                          "s(ecoregion_lvl)" = "s(ecoregion)")
  
  # get deviance explained from model summary and return as column after ecoregion with the value in the first row only
  combined <- combined |> 
    mutate(`Deviance Explained (%)` = ifelse(row_number() == 1, round(dev_expl * 100, 1), '-')) %>%
    dplyr::select(Response, `Deviance Explained (%)`, everything())
  # rename edf column (e)df
  colnames(combined)[colnames(combined) == "edf"] <- "(e)df"

  
  return(combined)
  }



