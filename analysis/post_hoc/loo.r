# Script for calculating loo estimates
library(dplyr)
library(kableExtra)
here::i_am("analysis/post_hoc/loo.r")
setwd(here::here())

outcomes <- c("heavy", "alc", "stim", "thc")
covar_models <- c("rifactor","ri", "ind", "un", "cs", "ar", "ad",  "indcv", "uncv", "cscv", "arcv", "adcv")

loo_df <- data.frame(outcome = character(), covar_model = character(),
                     looic = numeric(), elpd_loo = numeric(), p_loo = numeric())
#colnames(loo_df) <- c("outcome", "covar_model", "looic", "elpd_loo", "p_loo")

for (outcome in outcomes){
  for (covmod in covar_models){
    # Load the fit results
    loo_result <- readRDS(paste0("analysis/fit_results/loo/", outcome, "_", covmod, "_loo.rds"))
    # Calculate the loo estimate
    loo_df <- rbind(loo_df, data.frame(outcome = outcome, covar_model = covmod,
                                       looic = as.numeric(loo_result$estimates[3,1]),
                                       elpd_loo = as.numeric(loo_result$estimates[1,1]),
                                       p_loo = as.numeric(loo_result$estimates[2,1])))
  }
}


# Now turn into a table
loo_table <- loo_df |> select(outcome, covar_model, looic) |>
    tidyr::pivot_wider(names_from = outcome, values_from = looic) |>
    mutate(across(heavy:thc, ~round(.))) |>
    arrange(heavy) |> kableExtra::kable(format = 'latex', booktabs = T, format.args = list(big.mark = ","),
    caption = "LOOIC estimates for each outcome and covariate codel, sorted by fit for heavy drinking outcome.")

loo_table
