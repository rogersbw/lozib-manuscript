# R script
library(cmdstanr)

#Loading and tidying the data




#Make a function, that takes a model, and runs through the appropriate stan code

para_run <- function(X, ll, y, visit, V, model, iter, warmup, chains,  seed){
  N = length(y)
  P = ncol(X)
  J = length(unique(ll))
  
  zib_dat <- list(P = P, J = J, N = N, ll= ll, X = X, y = y, visit = visit, V=V)
  
  stanfile <- file.path(paste0("ZIB",model,".stan"))
  
  zib_stan_model <- cmdstan_model(stanfile)
  
  zib_fit <- zib_stan_model$sample(data = zib_dat,
                                   chains=chains, 
                                   parallel_chains=chains,
                                   iter_sampling=iter,
                                   iter_warmup=warmup,
                                   seed=seed,
                                   max_treedepth=15)  
  
  return(zib_fit)
}


