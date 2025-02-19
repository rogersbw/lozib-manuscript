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

  if(model_name %in% c("ri", "rifactor")){
    return(list(beta1 = beta1, beta2 = beta2, sigma1 = sigma1, sigma2 = sigma2, psi = psi, gamma1 = gamma1, gamma2 = gamma2))
  }
  if(model_name %in% c("ind", "indcv", "cs", "cscv", "ar", "arcv", "ad", "adcv", "un", "uncv")){
    gamma2_1 <- select(gamma2, ends_with(".1."))
    gamma2_2 <- select(gamma2, ends_with(".2."))
    gamma2_3 <- select(gamma2, ends_with(".3."))
    gamma2_4 <- select(gamma2, ends_with(".4."))
    return(list(beta1 = beta1, beta2 = beta2, sigma1 = sigma1, sigma2 = sigma2, psi = psi, gamma1 = gamma1, gamma2_1 = gamma2_1, gamma2_2 = gamma2_2, gamma2_3 = gamma2_3, gamma2_4 = gamma2_4))
  }
  else{
    stop("model not found")
  }
}


# To generate draws of zero-inflated model parameters
# Returns posterior samples for theta_est, pi_est, and mu_est
post_means <- function(outcome, model_name, out_of_sample = TRUE) {

  post_draws <- extract_draws(outcome, model_name)

  iter = dim(post_draws[[1]])[1]
  pi_est <- array(data = NA, c(iter, 7))
  theta_est <- array(data = NA, c(iter, 7))
  mu_est <- array(data = NA, c(iter, 7))
  
  gam1 <- rnorm(iter)

  if(model_name %in% c("ri", "rifactor")) {
    gam2 <- matrix(rep(rnorm(iter), 4), ncol = 4)
  }
  else {
    gam2 <- matrix(rnorm(4 * iter), ncol = 4)
  }
 
  beta1 <- post_draws$beta1
  beta2 <- post_draws$beta2
  sigma1 <- post_draws$sigma1
  sigma2 <- post_draws$sigma2
  psi <- post_draws$psi


  #For control group
  for(t in 1:4){
    theta_est[, t] <- pull(expit(beta1[,t] + sigma1 * gam1))
    pi_est[, t] <- pull(expit(beta2[,t] + psi * sigma1 * gam1 + sigma2 * gam2[,t]))
    mu_est[, t] <- theta_est[, t] * pi_est[, t] * 90
  }

  #For treatment group
  for(t in 1:3){
    theta_est[, t + 4] <- pull(expit(beta1[,t] + beta1[, t + 4] + sigma1 * gam1))
    pi_est[, t + 4] <- pull(expit(beta2[,t] + beta2[, t + 4] + psi * sigma1 * gam1 + sigma2 * gam2[,t]))
    mu_est[, t + 4] <- theta_est[, t + 4] * pi_est[, t + 4] * 90
  }

  return(list(theta_est = theta_est, pi_est = pi_est, mu_est = mu_est))
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

test <- post_means("heavy", "indcv")
test_post_pred <- post_predict(test$theta_est, test$pi_est)