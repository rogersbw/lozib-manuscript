require(dplyr)
# Post processing functions

here::i_am("analysis/post_hoc_analysis/post_functions.r")

proj_path <- here::here()



# Logit function
logit <- function(x) { 
  log(x / (1 - x))
}

# Inverse logit function
expit <- function(x) {
  exp(x) / (1 + exp(x))
}


## Function takes the stan draws and extracts the relevant parameters for the model
extract_draws <- function(outcome, model_name){
  if(!is.character(model_name) || !is.character(outcome)){
      stop("model_name and outcome name must be character strings")
  }
  draws_path <- file.path(proj_path,paste0("analysis/fit_results/parameter_draws/draws_", outcome, "_", model_name, ".rds"))
  post <- readRDS(draws_path)
  chains <- dim(post)[2]
  post_df <- data.frame(post[,1,])
  for (i in 2:chains){
    post_df <- rbind(post_df, data.frame(post[,i,]))
  }
  beta1 <- select(post_df, starts_with("beta1"))
  beta2 <- select(post_df, starts_with("beta2"))
  sigma1 <- select(post_df, starts_with("sigma1"))
  sigma2 <- select(post_df, starts_with("sigma2"))
  psi <- select(post_df, starts_with("psi"))
  gamma1 <- select(post_df, starts_with("gamma1"))
  gamma2 <- select(post_df, starts_with("gamma2"))

  params <- list(beta1 = beta1, beta2 = beta2, sigma1 = sigma1, sigma2 = sigma2, psi = psi, gamma1 = gamma1)

  if (model_name %in% c("ri", "rifactor")){
    params <- append(params, list(gamma2 = gamma2))
  }
  if (model_name %in% c("ind", "indcv", "cs", "cscv", "ar", "arcv", "ad", "adcv", "un", "uncv")){
    gamma2_1 <- select(gamma2, ends_with(".1."))
    gamma2_2 <- select(gamma2, ends_with(".2."))
    gamma2_3 <- select(gamma2, ends_with(".3."))
    gamma2_4 <- select(gamma2, ends_with(".4."))

    params <- append(params, list(gamma2_1 = gamma2_1, gamma2_2 = gamma2_2, gamma2_3 = gamma2_3, gamma2_4 = gamma2_4))

    if (model_name %in% c("ar", "arcv", "ad", "adcv")){
        rho <- select(post_df, starts_with("rho"))

        params <- append(params, list(rho = rho))
    }
    if (model_name %in% c("un", "uncv")){
        L_omega <- select(post_df, starts_with("L_omega"))

        params <- append(params, list(L_omega = L_omega))
    }
  }
  return(params)
}

# Function to transform gamma draws for cs model
transform_cs <- function(gam2, sigma2, sigma2_cs) {
  gam2_scaled <- list()
  gam2_cs <- rnorm(dim(gam2)[1])
  gam2_cs <- outer(sigma2_cs, gam2_cs)
  for (t in 1:dim(gam2)[2]){
    gam2_scaled[[t]] <- outer(sigma2[,t], gam2[,t]) + gam2_cs 
}
  return(gam2_scaled)
}

transform_ind <- function(gam2, sigma2) {
  gam2_scaled <- list()
  for (t in 1:dim(gam2)[2]){
    gam2_scaled[[t]] <- outer(sigma2[,t], gam2[,t])
  }
  return(gam2_scaled)
}

transform_ar <- function(gam2, sigma2, rho, cv = FALSE) {
  gam2_scaled <- list()
  if(cv) {
    sigma2[, 1] <- sigma2[, 1]/sqrt(1-rho^2)
  }
  gam2_scaled[[1]] <- outer(sigma2[,1], gam2[,1]) 
  for (t in 2:dim(gam2)[2]){
    gam2_scaled[[t]] <- rho * gam2_scaled[[t-1]] + outer(sigma2 [, t], gam2[,t])
  }
  return(gam2_scaled)
}

transform_ad <- function(gam2, sigma2, rho, cv = FALSE) {
  gam2_scaled <- list()
  if(cv) {
    sigma2[, 1] <- sigma2[, 1]/sqrt(1-rho[,1]^2)
  }
  gam2_scaled[[1]] <- outer(sigma2[,1], gam2[,1]) 
  for (t in 2:dim(gam2)[2]){
    gam2_scaled[[t]] <- rho[, t-1] * gam2_scaled[[t-1]] + outer(sigma2 [, t], gam2[,t])
  }
  return(gam2_scaled)
}

#transform_un <- function(gam2, sigma2, L_omega, cv=FALSE) {
# In this we probably want to got through the gam samples one by one and chol multiply each one,
# We can build out an array, n_post, n_gamma_samples, n_time_points
# Should switch the others to also return this array.
#}


# To generate draws of zero-inflated model parameters
# Returns posterior samples for theta_est, pi_est, and mu_est
# Need to construct AR, AD and UN gammas
post_means <- function(outcome, model, gam_draws = 500, out_of_sample = TRUE) {

  post_draws <- extract_draws(outcome, model)

  iter <- dim(post_draws[[1]])[1]
  pi_est <- array(data = NA, c(iter, gam_draws, 7))
  theta_est <- array(data = NA, c(iter, gam_draws, 7))
  mu_est <- array(data = NA, c(iter, gam_draws, 7))
  gam1 <- rnorm(gam_draws)

  if (model %in% c("ri", "rifactor")) {
    gam2 <- matrix(rep(rnorm(gam_draws), 4), ncol = 4)
  }else {
    gam2 <- matrix(rnorm(4 * gam_draws), ncol = 4)
  }
 
  beta1 <- post_draws$beta1
  beta2 <- post_draws$beta2
  sigma1 <- post_draws$sigma1
  sigma2 <- post_draws$sigma2
  psi <- post_draws$psi

  if (grepl("cv$", model) || grepl("ri", model)) {
    sigma2 <- matrix(rep(sigma2[, 1], 4), ncol = 4) 
  }

  gam1_scaled <- outer(sigma1[, 1], gam1)
  psi_scaled <- outer(psi[, 1], gam1)

  #For control group
  for (t in 1:4){
    theta_est[, , t] <- expit(beta1[,t] + gam1_scaled)
    pi_est[, , t] <- expit(beta2[,t] + psi_scaled + outer(sigma2[, t], gam2[, t]))
    mu_est[, , t] <- theta_est[, , t] * pi_est[, , t] * 90
  }

  #For treatment group
  for (t in 1:3){
    theta_est[, , t + 4] <- expit(beta1[, t + 1] + beta1[, t + 4] + gam1_scaled)
    pi_est[, , t + 4] <- expit(beta2[,t + 1] + beta2[, t + 4] + psi_scaled + outer(sigma2[, t + 1], gam2[, t + 1]))
    mu_est[, , t + 4] <- theta_est[, , t + 4] * pi_est[, , t + 4] * 90
  }
  
  return(list(theta_est = theta_est, pi_est = pi_est, mu_est = mu_est))
}


## Take an array of posterior samples (integrating over random effects) and return the mean and 95% credible interval for the mean
post_summary_means <- function(est_array) {
    means <- apply(est_array, c(1, 3), mean)
    avg <- apply(means, 2, mean)
    lower <- apply(means, 2, quantile, probs = c(.025))
    upper <- apply(means, 2, quantile, probs = c(.975))
    return(list(avg = avg, lower = lower, upper = upper))
}


calc_DoD <- function(est) {
    DoD <- array(data = NA, c(dim(est)[1], dim(est)[2], 3))
    DoD[, , 1] = est[, , 5] - est[, , 2]
    DoD[, , 2] = est[, , 6] - est[, , 3]
    DoD[, , 3] = est[, , 7] - est[, , 4]

    DoD_summary <- post_summary_means(DoD)

    return(DoD_summary)
}

# posterior predictive draws

post_predict <- function(theta_draws, pi_draws) {
    iter <- dim(theta_draws)[1]
    y_draws <- array(data = NA, c(iter, 7))
    for(i in 1:iter){
        y_draws[i,] <- rbinom(7, 1, theta_draws[i,]) * rbinom(7, 90, pi_draws[i,])
    }
    return(y_draws)
}
