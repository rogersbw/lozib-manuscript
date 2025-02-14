# To simulate from posterior predictive distribution
library(dplyr)
library(ggplot2)
# We will do posterior predictive distribution from IND model and RI model

# Two auxilary functions
logit<-function(x){log(x/(1-x))}
expit<-function(x){exp(x)/(1+exp(x))}

# For in-sample predictions, we will need the data

sbirt <- readRDS("analysis/data/sbirt_clean.rds")
features_ctrl <- sbirt |> filter(group == 0) |> select(id, starts_with("visit"))
features_trt <- sbirt |> filter(group == 1) |> select(id, starts_with("visit"))


##############################################################################
######### Posterior Draws
##############################################################################

# First load in INDcv and RI model results

indcv_post <- readRDS("analysis/fit_results/parameter_draws/draws_heavy_indcv.rds")
ri_post <- readRDS("analysis/fit_results/parameter_draws/draws_heavy_ri.rds")


chains <- dim(indcv_post)[2]

indcv_post_df <- data.frame(indcv_post[,1,])
for (i in 2:chains){
  indcv_post_df <- rbind(indcv_post_df, data.frame(indcv_post[,i,]))
}


beta1 <- indcv_post_df[,1:7]
beta2 <- indcv_post_df[,8:14]
sigma1 <- indcv_post_df[,15]
sigma2 <- indcv_post_df[,16]
psi <- indcv_post_df[,17]
gamma1 <- indcv_post_df[,18:735 ]
gamma2_1 <- indcv_post_df[,736:1453]
gamma2_2 <- indcv_post_df[,1454:2171]
gamma2_3 <- indcv_post_df[,2172:2889]
gamma2_4 <- indcv_post_df[,2890:3607]

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

sample_plots_baseline <- y_draws_baseline_indcv_df |> 
    filter(sample < 11) |> 
    mutate(sample = as.factor(sample)) |>
    ggplot() + geom_histogram(aes(x = y), binwidth = 1) + facet_wrap(vars(sample))

sample_plots_baseline



# Fix up these plots and save them for appendix


# Now do the same for the RI model


# Then means at each time point for in sample ()

