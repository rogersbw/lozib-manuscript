# Script to produce posterior summary data set

library(dplyr)

### Directory organization

#here::i_am("analysis/post_hoc/post_predict.r")

#proj_path <- here::here()

source("post_functions.r")

model_list <- c("ri", "ind", "indcv", "cs", "cscv", "ar", "arcv", "ad", "adcv", "un", "uncv")
outcomes <- c("heavy", "alc", "stim", "thc")
visits <- c("visit1", "visit2", "visit3", "visit4", "visit2trt", "visit3trt", "visit4trt")

post_summary_df <- structure(list(model = character(), outcome = character(), est = character(), avg = numeric(), lower = numeric(), upper = numeric()), class = "data.frame")
  
for (outcome in outcomes) {
  for (model in model_list) {
    post <- post_means(outcome, model)
    tryCatch({
    post_summary <- post_summary_means(post$mu_est)
    },  error = function(e) {
      message("Error in model: ", model, ", outcome: ", outcome)
      message("Error message: ", e$message)
      next  # Skip to the next iteration
    })
    post_summary_df <- rbind(post_summary_df,
                             data.frame(model = model,
                                        outcome = outcome,
                                        visit = visits,
                                        est = "avg",
                                        avg = post_summary$avg,
                                        lower = post_summary$lower,
                                        upper = post_summary$upper))
    tryCatch({post_summary <- post_summary_means(post$theta_est)},  error = function(e) {
      message("Error in model: ", model, ", outcome: ", outcome)
      message("Error message: ", e$message)
      next  # Skip to the next iteration
    })
    post_summary_df <- rbind(post_summary_df,
                             data.frame(model = model,
                                        outcome = outcome,
                                        visit = visits,
                                        est = "theta",
                                        avg = post_summary$avg,
                                        lower = post_summary$lower,
                                        upper = post_summary$upper))

    tryCatch({
    post_summary <- post_summary_means(post$pi_est)},  error = function(e) {
      message("Error in model: ", model, ", outcome: ", outcome)
      message("Error message: ", e$message)
      next  # Skip to the next iteration
    })
    post_summary_df <- rbind(post_summary_df,
                             data.frame(model = model,
                                        outcome = outcome,
                                        visit = visits,
                                        est = "pi",
                                        avg = post_summary$avg,
                                        lower = post_summary$lower,
                                        upper = post_summary$upper))

    post_summary <- calc_DoD(post$mu_est)
    post_summary_df <- rbind(post_summary_df,
                             data.frame(model = model,
                                        outcome = outcome,
                                        visit = visits[2:4],
                                        est = "DoD",
                                        avg = post_summary$avg,
                                        lower = post_summary$lower,
                                        upper = post_summary$upper))
    
    rm(post)
  }
}


saveRDS(post_summary_df, file.path("..", "data", "processed", "post_summary_means.rds"))
