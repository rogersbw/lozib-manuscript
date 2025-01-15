#Switching this to parallelize in the sbatch and command line rather than through R
#The hope is this allows for each cmdstan run to multicore
library(loo)

source("ZIBModels.R")

args <- commandArgs(trailingOnly = TRUE)

covmod <- args[1]

iter <- 20000
warmup <- 5000
chains <- 4

sbirt <- readRDS("sbirt_clean.rds")

# Fitting the Model

idx <- !is.na(sbirt$heavy)

P = 7 # Number of predictor columns (visits*treatment groups)
J = length(unique(sbirt$id)) # Unique subjects
N = nrow(sbirt[idx,]) # Number of total observations
ll = sbirt$id[idx] # List of subject ids
X = as.matrix(sbirt[idx,8:14]) # Predictor variables
y = sbirt$heavy[idx] # Outcome
visit = sbirt$visit[idx] # Visit number
V = 4 # Max number of visits


fit_result <- para_run(X, ll, y, visit, V=V, model=covmod, iter, warmup, chains,  seed=1234)

## This code I'm just editing to run locally, should be removed for cluster runs
# 
# zib_stan_model <- cmdstan_model("ZIBrifactor.stan")
# covmod = "rifactor"
# 
# zib_dat <- list(P = P, J = J, N = N, ll= ll, X = X, y = y, visit = visit, V=V)
# 
# zib_fit <- zib_stan_model$sample(data = zib_dat,
#                                  chains=chains, 
#                                  parallel_chains=chains,
#                                  iter_sampling=iter,
#                                  iter_warmup=warmup,
#                                  seed=1234,
#                                  max_treedepth=15)  
# 
# out.file <- file.path(paste0("output_",covmod,".rds"))
# 
# zib_fit$loo()

####

out.file <- file.path(paste0("output_",covmod,".rds"))

saveRDS(fit_result, out.file)

fit_result$loo(cores = 3)

fit_result$summary(variables = c("beta1","beta2", "sigma1", "sigma2", "lp__"))

print(covmod)




#### Old version
# 
# library(doParallel) #To run all models in parallel
# 
# source("ZIBModels.R")
# 
# iter <- 10000
# warmup <- 3000
# chains <- 4
# 
# sbirt <- readRDS("sbirt_clean.rds")
# 
# # Fitting the Model
# 
# idx <- !is.na(sbirt$heavy)
# 
# P = 7 # Number of predictor columns (visits*treatment groups)
# J = length(unique(sbirt$id)) # Unique subjects
# N = nrow(sbirt[idx,]) # Number of total observations
# ll = sbirt$id[idx] # List of subject ids
# X = as.matrix(sbirt[idx,8:14]) # Predictor variables
# y = sbirt$heavy[idx] # Outcome
# visit = sbirt$visit[idx] # Visit number
# V = 4 # Max number of visits
# 
# 
# 
# #Test for just ri model
# #model_outs <- para_run(X, ll, y, visit, V=V, model="ar", iter, warmup, chains,  seed=1234)
# 
# # 
# models <- c("ri", "arcv", "ar", "adcv", "ad", "uncv", "un")
# 
# # 
# # 
# registerDoParallel(length(models))
# 
# # 
# model_outs <- foreach(i=1:length(models)) %dopar% para_run(X, ll, y, visit, V=V, model=models[i], iter, warmup, chains,  seed=(1234*i))
# 
# saveRDS(model_outs, "model_output.rds")



