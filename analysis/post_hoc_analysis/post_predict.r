# To simulate from posterior predictive distribution
library(dplyr)
library(ggplot2)

### Directory organization

here::i_am("analysis/post_hoc_analysis/post_predict.r")

proj_path <- here::here()

source(file.path(proj_path, "analysis", "post_hoc_analysis", "post_functions.r"))



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