# Function to run estimated values stan models


library(cmdstanr)

#Loading and tidying the data







#Make a function, that takes a model, and runs through the appropriate stan code

est_run <- function(X, ll, y, visit, V, model, iter, warmup, chains,  seed){
  N = length(y)
  P = ncol(X)
  J = length(unique(ll))
  
  zib_dat <- list(P = P, J = J, N = N, ll= ll, X = X, y = y, visit = visit, V=V)
  
  stanfile <- file.path(paste0("ZIB",model,"_est.stan"))
  
  zib_stan_model <- cmdstan_model(stanfile)
  
  zib_fit <- zib_stan_model$sample(data = zib_dat,
                                   chains=chains, 
                                   parallel_chains=chains,
                                   iter_sampling=iter,
                                   iter_warmup=warmup,
                                   seed=seed,
                                   max_treedepth=15)  
  
  params <- c("beta1", "beta2", "sigma1", "sigma2")
  
  if(model %in% c("ri", "ind", "indcv")){
    model_draws <- zib_fit$draws(c(params, "psi"))
  }else if(model %in% c("un", "uncv")){
    model_draws <- zib_fit$draws(c(params, "psi", "Omega2"))
  }else if(model %in% c("cs", "cscv")){
    model_draws <- zib_fit$draws(c(params, "psi", "sigma2_cs"))
  }else if(model %in% c("ar", "arcv", "ad", "adcv")){
    model_draws <- zib_fit$draws(c(params, "psi", "rho"))
  }else{
    model_draws <- zib_fit$draws(c(params))
  }
  
  return(model_draws)
}
