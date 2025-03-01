# To actually run the lozib model code
options(mc.cores = parallel::detectCores())

require("loo")
source("lozib/fit_lozib.R")

args <- commandArgs(trailingOnly = TRUE)

covmod <- args[1]
outcome <- args[2]

iter <- 40000
warmup <- 10000
chains <- 4

sbirt <- readRDS("data/sbirt_clean.rds")

# Fitting the Model

idx <- !is.na(sbirt[[outcome]])

P <- 7 # Number of predictor columns (visits*treatment groups)
J <- length(unique(sbirt$id)) # Unique subjects
N <- nrow(sbirt[idx,]) # Number of total observations
ll <- sbirt$id[idx] # List of subject ids
X <- as.matrix(sbirt[idx, 8:14]) # Predictor variables
y <- sbirt[[outcome]][idx] # Outcome
visit <- sbirt$visit[idx] # Visit number
V <- 4 # Max number of visits

#Fit the stan model
fit_result <- fit_lozib(X, ll, y, visit, V = V, model = covmod,
                        iter, warmup, chains,  seed = 9731)

# Parse results into files of interest

loo_result <- fit_result$loo(cores = 3)
saveRDS(loo_result,
        file.path(paste0("fit_results/loo/", outcome, "_", covmod, "_loo.rds")))

# Select which parameters to keep draws for
# We want these four parameters for all models
params <- c("beta1", "beta2", "sigma1", "sigma2")
# Then depending on the mode, we want different additional parameters
if (covmod %in% c("ri", "ind", "indcv")) {
  params <- c(params, "psi")
}else if (covmod %in% c("un", "uncv")) {
  params <- c(params, "psi", "L_Omega")
}else if (covmod %in% c("cs", "cscv")) {
  params <- c(params, "psi", "sigma2_cs")
}else if (covmod %in% c("ar", "arcv", "ad", "adcv")) {
  params <- c(params, "psi", "rho")
}
# Same as above, but for the random effects
if (covmod == "rifactor"){
  params <- c(params, "gamma1")
}else {
  params <- c(params, "gamma1", "gamma2")
}

# Extract relevant draws
model_draws <- fit_result$draws(c(params))
#Save to file
saveRDS(model_draws,
        paste0("fit_results/parameter_draws/draws_",
               outcome, "_", covmod, ".rds"))
