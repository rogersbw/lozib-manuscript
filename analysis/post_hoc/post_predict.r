# To simulate from posterior predictive distribution
library(dplyr)
library(ggplot2)

### Directory organization

here::i_am("analysis/post_hoc/post_predict.r")

proj_path <- here::here()

source(file.path(proj_path, "analysis", "post_hoc", "post_functions.r"))



# For in-sample predictions, we will need the data

sbirt <- readRDS(file.path(proj_path, "analysis/data/sbirt_clean.rds"))
features_ctrl <- sbirt |> filter(group == 0) |> select(id, starts_with("visit"))
features_trt <- sbirt |> filter(group == 1) |> select(id, starts_with("visit"))


##############################################################################
######### Posterior Draws
##############################################################################

##################
#### INDcv model
##################

set.seed(9731)


### First let's look at posterior mean estimates from RI and INDcv models

indcv_post <- post_means("heavy", "indcv")
indcv_mean <- post_summary_means(indcv_post$mu_est)
indcv_mean

paste0(round(indcv_mean$avg, 2), " (", round(indcv_mean$lower, 2), ", ", round(indcv_mean$upper, 2), ")")

rm(indcv_post)

ri_post <- post_means("heavy", "ri")
ri_mean <- post_summary_means(ri_post$mu_est)

paste0(round(ri_mean$avg, 2), " (", round(ri_mean$lower, 2), ", ", round(ri_mean$upper, 2), ")")


rm(ri_post)

ar_post <- post_means("heavy", "ar")
ar_mean <- post_summary_means(ar_post$mu_est)

# Observed Data means:

# Baseline mean
baseline <- sbirt |> filter(visit == 1)
mean(baseline$heavy, na.rm = TRUE)

# Follow-up means
sbirt |> filter(visit > 1) |> group_by(visit, group) |> summarise(mean_heavy = mean(heavy, na.rm = TRUE))


### Function to get plot of all models
## Step 1: Create data frame with model, outcome, mean, lb, and ub

model_list <- c("ri", "rifactor", "ind", "indcv", "cs", "cscv", "ar", "arcv", "ad", "adcv", "un", "uncv")
outcomes <- c("heavy", "alc", "stim", "thc")
visits <- c("visit1", "visit2", "visit3", "visit4", "visit2trt", "visit3trt", "visit4trt")
post_summary_df <- structure(list(model = character(), outcome = character(), est = character(), avg = numeric(), lower = numeric(), upper = numeric()), class = "data.frame")

for (outcome in outcomes) {
  for (model in model_list){
    post <- post_means(outcome, model)
    post_summary <- post_summary_means(post$mu_est)
    post_summary_df <- rbind(post_summary_df,
                             data.frame(model = model,
                                        outcome = outcome,
                                        est = visits,
                                        avg = post_summary$avg,
                                        lower = post_summary$lower,
                                        upper = post_summary$upper))
  
  }
}


## Let's get a data frame of mean, lower, and upper bounds for each model and each substance

indcv_DoD <- calc_DoD(indcv_post$mu_est)


#indcv_post_predict <- data.frame(post_predict(indcv_post$theta_est, indcv_post$pi_est))
colnames(indcv_post_predict) <- c("visit1", "visit2", "visit3", "visit4", "visit2trt", "visit3trt", "visit4trt")

ggplot(data = indcv_post_predict, aes(x = visit1)) + 
  geom_density(color = "darkblue", fill = "lightblue")



# First load in INDcv and RI model results
indcv_post <- extract_draws("heavy", "indcv")

#Zero model fixed effects
visit1_ids <- sbirt |> filter(visit == 1) |> select(id) |> unique() |> pull()
visit2_ids <- sbirt |> filter(visit == 2) |> select(id) |> unique() |> pull()
id_group_map <- sbirt |> select(id, group) |> unique()


# Select index of posterior samples
draw_iter <- seq(1, dim(indcv_post_df)[1], by = 500)

y_draws_baseline_indcv_df = data.frame(sample = factor(), y = integer())
iter <- 1

for(i in draw_iter){

    theta_visit1 <- expit(beta1[i,1] + sigma1[i]*as.numeric(gamma1[i,visit1_ids]))
    pie_visit1 <- expit(beta2[i,1] + psi[i]*sigma1[i] + sigma2[i]*as.numeric(gamma2_1[i,visit1_ids]))

    y_draws_temp <- rbinom(length(theta_visit1), 1, theta_visit1) * rbinom(length(pie_visit1), 90, pie_visit1)
    sample <- rep(iter, length(y_draws_temp))
    y_draws_baseline_indcv_df <- rbind(y_draws_baseline_indcv_df, data.frame(sample = sample, y = y_draws_temp))
    iter <- iter + 1
}

sample_plots_baseline_indcv <- y_draws_baseline_indcv_df |> 
    filter(sample < 11) |> 
    mutate(sample = as.factor(sample)) |>
    ggplot() + 
    geom_histogram(aes(x = y), binwidth = 1) +
    facet_wrap(vars(sample)) +


sample_plots_baseline_indcv

##################
#### RI model
##################


ri_post <- readRDS("analysis/fit_results/parameter_draws/draws_heavy_ri.rds")

chains <- dim(ri_post)[2]

ri_post_df <- data.frame(ri_post[,1,])
for (i in 2:chains){
  ri_post_df <- rbind(ri_post_df, data.frame(ri_post[,i,]))
}

beta1 <- ri_post_df[,1:7]
beta2 <- ri_post_df[,8:14]
sigma1 <- ri_post_df[,15]
sigma2 <- ri_post_df[,16]
psi <- ri_post_df[,17]
gamma1 <- ri_post_df[,18:735 ]
gamma2 <- ri_post_df[,736:1453]



# Select index of posterior samples
draw_iter <- seq(1, dim(ri_post_df)[1], by = 500)

y_draws_baseline_ri_df = data.frame(sample = factor(), y = integer())
iter <- 1

for(i in draw_iter){

    theta_visit1 <- expit(beta1[i,1] + sigma1[i]*as.numeric(gamma1[i,visit1_ids]))
    pie_visit1 <- expit(beta2[i,1] + psi[i]*sigma1[i] + sigma2[i]*as.numeric(gamma2[i,visit1_ids]))

    y_draws_temp <- rbinom(length(theta_visit1), 1, theta_visit1) * rbinom(length(pie_visit1), 90, pie_visit1)
    sample <- rep(iter, length(y_draws_temp))
    y_draws_baseline_ri_df <- rbind(y_draws_baseline_ri_df, data.frame(sample = sample, y = y_draws_temp))
    iter <- iter + 1
}



sample_plots_baseline_ri <- y_draws_baseline_ri_df |> 
    filter(sample < 11) |> 
    mutate(sample = as.factor(sample)) |>
    ggplot() + geom_histogram(aes(x = y), binwidth = 1) + facet_wrap(vars(sample))

sample_plots_baseline_ri

# Then means at each time point for in sample ()

# Need to touch up these plots

# Also want posterior distribution of the ranom intercepts


# Mean of treatment and control group at baseline and each time point
# What is the posterior distribution of group means


##########
# Density plots of Random effects distributions
# Predictive mean and variance for RI and INDcv
# 