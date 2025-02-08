# Description: Function to fit the lozib model to data
here::i_am("analysis/lozib/fit_lozib.R")

fit_lozib <- function(X, ll, y, visit, V, model, iter, warmup, chains,  seed) {
  # Set up stan data list
  require("cmdstanr")

  N <- length(y)
  P <- ncol(X)
  J <- length(unique(ll))
  zib_dat <- list(P = P, J = J, N = N, ll = ll,
                  X = X, y = y, visit = visit, V = V)
  # Set up stan model
  stanfile <- file.path(paste0("analysis/lozib/stan_models/ZIB", model, ".stan"))
  zib_stan_model <- cmdstan_model(stanfile)
  # Fit the model
  zib_fit <- zib_stan_model$sample(data = zib_dat,
                                   chains = chains,
                                   parallel_chains = chains,
                                   iter_sampling = iter,
                                   iter_warmup = warmup,
                                   thin = 2,
                                   seed = seed,
                                   max_treedepth = 15)
  return(zib_fit)
}