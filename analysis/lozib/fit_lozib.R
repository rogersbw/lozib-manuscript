# Description: Function to fit the lozib model to data
here::i_am("analysis/lozib/fit_lozib.R")

############
# Load data
############

#I'd like args to consists of: loo or est and the name of the model


fit_lozib <- function(){
    #This function will run everything
}


fit_lozib <- function(X, ll, y, visit, V, model, iter, warmup, chains,  seed) {
  # Set up stan data list
  require("cmdstanr")

  N <- length(y)
  P <- ncol(X)
  J <- length(unique(ll))
  zib_dat <- list(P = P, J = J, N = N, ll = ll,
                  X = X, y = y, visit = visit, V = V)
  # Set up stan model
  stanfile <- file.path(paste0("ZIB", model, ".stan"))
  zib_stan_model <- cmdstan_model(stanfile)
  # Select which parameters to keep draws for
  # For all models we want these 4 parameters
  params <- c("beta1", "beta2", "sigma1", "sigma2")
  # Then depending on the mode, we want different additional parameters
  if (model %in% c("ri", "ind", "indcv")) {
    params <- c(params, "psi")
  }else if (model %in% c("un", "uncv")) {
    params <- c(params, "psi", "Omega2")
  }else if (model %in% c("cs", "cscv")) {
    params <- c(params, "psi", "sigma2_cs")
  }else if (model %in% c("ar", "arcv", "ad", "adcv")) {
    params <- c(params, "psi", "rho")
  }

  if (model == "rifactor"){
    params <- c(params, "gamma1")
  }else {
    params <- c(params, "gamma1", "gamma2")
  }
  # If we want to calculate loo, we need to keep the log_lik
  params <- c(params, "log_lik")


  # Fit the model
  zib_fit <- zib_stan_model$sample(data = zib_dat,
                                   chains = chains,
                                   pars = params,
                                   parallel_chains = chains,
                                   iter_sampling = iter,
                                   iter_warmup = warmup,
                                   thin = 2,
                                   seed = seed,
                                   max_treedepth = 15)
  return(zib_fit)
}

