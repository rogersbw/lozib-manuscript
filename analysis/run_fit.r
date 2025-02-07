# To actually run the lozib model code
library(loo)

here::i_am("analysis/lozib/run_fit.R")

source("fit_lozib.R")

args <- commandArgs(trailingOnly = TRUE)

covmod <- args[1]

iter <- 20000
warmup <- 5000
chains <- 4

sbirt <- readRDS("analysis/data/sbirt_clean.rds")

# Fitting the Model

idx <- !is.na(sbirt$heavy)

P <- 7 # Number of predictor columns (visits*treatment groups)
J <- length(unique(sbirt$id)) # Unique subjects
N <- nrow(sbirt[idx,]) # Number of total observations
ll <- sbirt$id[idx] # List of subject ids
X <- as.matrix(sbirt[idx,8:14]) # Predictor variables
y <- sbirt$heavy[idx] # Outcome
visit <- sbirt$visit[idx] # Visit number
V <- 4 # Max number of visits

#Fit the stan model
fit_result <- fit_lozib(X, ll, y, visit, V=V, model=covmod, iter, warmup, chains,  seed=1234)

# Parse results into files of interest

loo_result <- fit_result$loo(cores = 3)
saveRDS(loo_result, file.path(paste0(covmod,"_loo.rds")))

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


out.file <- file.path(paste0("output_",covmod,".rds"))

saveRDS(fit_result, out.file)

fit_result$loo(cores = 3)

fit_result$summary(variables = c("beta1","beta2", "sigma1", "sigma2", "lp__"))

print(covmod)
