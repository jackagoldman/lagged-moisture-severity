
#' Title: Brief one-line title
#'
#' @description
#' Short description of what the function does. Provide enough detail for a reader
#' to understand the function's purpose and behavior.
#'
#' @param x <type> Description of parameter x: expected values, length, units, or structure.
#' @param y <type> Description of parameter y: allowed options and how it affects results.
#' @param ... Additional arguments passed to methods or lower-level functions.
#'
#' @return <type> Description of the returned value. If returning a list or data.frame,
#'   document important fields/columns and their types.
#'
#' @details
#' More extensive notes on algorithm, side effects (e.g. modifies global state, files),
#' NA/NULL handling, performance characteristics, and any invariants or guarantees.
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' res <- your_function(x = 1, y = "option")
#'
#' # Example showing edge-case behavior
#' }
#'
#' @seealso
#' \code{\link[base]{other_function}}, other_related_function
#'
#' @author Your Name
#' @export
#' @keywords internal
#' 
#' 
make_gam_summary_table <- function(model) {
  require(gt)
  require(dplyr)
  p_table <- summary(model)$p.table
  s_table <- summary(model)$s.table
  dev_expl <- summary(model)$dev.expl

  # get model family 
  model_family <- family(model)$family
  
  # Convert to data frames
  p_df <- as.data.frame(p_table)
  p_df <- rownames_to_column(p_df, var = "Term")
  
  s_df <- as.data.frame(s_table)
  s_df<- rownames_to_column(s_df, var = "Term")

  # add ecoregion column
  s_df$ecoregion <- NA
  s_df$ecoregion[1:12] <- sub(".*:", "", s_df$Term[1:12])

  # if model name contains cmi, recode ecoregion values based on weights_btl2, weights_ln2, weights_lsj2 else keep it the same
  if (grepl("cmi", deparse(substitute(model)))) {
    s_df$ecoregion <- recode(s_df$ecoregion, 
                            "weights_btl2" = "BTL", 
                            "weights_ln2" = "LN", 
                            "weights_lsj2" = "LSJ")
  } else {
    s_df$ecoregion <- recode(s_df$ecoregion, 
                            "weights_btl" = "BTL", 
                            "weights_ln" = "LN", 
                            "weights_lsj" = "LSJ")} 
  
  
  s_df$ecoregion <- recode(s_df$ecoregion, 
                            "ecoregion_lvlBTL" = "BTL", 
                            "ecoregion_lvlLN" = "LN", 
                            "ecoregion_lvlLSJ" = "LSJ")
  

  
  # move to before term
  s_df <- s_df |> dplyr::select(ecoregion, everything())

  # make <NA> to -
  s_df$ecoregion[is.na(s_df$ecoregion)] <- "-"

  #turn off scientific notation
  options(scipen = 999)

  # round all columns but if p-value is <0.001 make <0.001
  if (data.frame(model_family) ==  "gaussian"){
  p_df <- p_df %>%
    mutate(`Pr(>|t|)` = ifelse(`Pr(>|t|)` < 0.001, "<0.001", round(`Pr(>|t|)`, 3))) %>%
    mutate(across(c(Estimate, `Std. Error`, `t value`), ~ round(., 3)))
    # change p_df z value to Test statistic and Pr(>|z|) to p-value
  colnames(p_df)[colnames(p_df) == "t value"] <- "Test statistic"
  colnames(p_df)[colnames(p_df) == "Pr(>|t|)"] <- "p-value"
  } else {
    p_df <- p_df %>%
      mutate(`Pr(>|z|)` = ifelse(`Pr(>|z|)` < 0.001, "<0.001", round(`Pr(>|z|)`, 3))) %>%
      mutate(across(c(Estimate, `Std. Error`, `z value`), ~ round(., 3)))
    # change p_df z value to Test statistic and Pr(>|z|) to p-value
    colnames(p_df)[colnames(p_df) == "z value"] <- "Test statistic"
    colnames(p_df)[colnames(p_df) == "Pr(>|z|)"] <- "p-value"
  }

  if (data.frame(model_family)== "gaussian"){
  s_df <- s_df %>%
    mutate(`p-value` = ifelse(`p-value` < 0.001, "<0.001", round(`p-value`, 3))) %>%
    mutate(across(c(edf, Ref.df, F), ~ round(., 3))) } else {
    s_df <- s_df %>%
      mutate(`p-value` = ifelse(`p-value` < 0.001, "<0.001", round(`p-value`, 3))) %>%
      mutate(across(c(edf, Ref.df, Chi.sq), ~ round(., 3)))
  }

  # add edf column with values of one before Test statistic and only keep edf, Test statistic, p-value
  p_df <- p_df %>%
    mutate(edf = 1) %>%
    dplyr::select(Term,  edf, `Test statistic`, `p-value`)

  # remove (Intercept) row in Term
  p_df <- p_df %>% filter(Term != "(Intercept)")
  # take everything from rowname before the , and add a )
  s_df$Term <- paste0(sub("(.*),.*", "\\1", s_df$Term), ")")

  # for the last two row names and remove the last )
  s_df$Term[(nrow(s_df)-1):nrow(s_df)] <- sub("\\)$", "", s_df$Term[(nrow(s_df)-1):nrow(s_df)])
  
  # remove ref.df column
  s_df <- s_df |> dplyr::select(-Ref.df)

  # remove everything after :from term names on rows 4 -12
  s_df$Term[4:12] <- sub(":.*", "", s_df$Term[4:12])

  # rename Chi.sq to Test statistic
  if (data.frame(model_family) == "gaussian"){
  colnames(s_df)[colnames(s_df) == "F"] <- "Test statistic"
  } else {
  colnames(s_df)[colnames(s_df) == "Chi.sq"] <- "Test statistic"
  }
  # to p_df add ecoregion column with values of -
  if (nrow(p_df) == 0) {
    combined <- s_df
  } else {
    p_df$ecoregion <- "-"
    p_df <- p_df |> dplyr::select(ecoregion, everything())
    combined <- dplyr::bind_rows(p_df, s_df)
  }
  response <- as.character(model$formula)[2]

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
                          "s(fwi_90th)" = "s(FWI)", 
                          "s(mean_evi)" = "s(EVI)", 
                          "s(fyear)" = "s(fire year)",
                          "s(ecoregion_lvl)" = "s(ecoregion)")
  
  # get deviance explained from model summary and return as column after ecoregion with the value in the first row only
  combined <- combined |> 
    mutate(`Deviance Explained (%)` = ifelse(row_number() == 1, round(dev_expl * 100, 1), '-')) %>%
    dplyr::select(Response, ecoregion, `Deviance Explained (%)`, everything())
  # rename edf column (e)df
  colnames(combined)[colnames(combined) == "edf"] <- "(e)df"

  
  return(combined)
  }



