# Update this to be able to run all models
# Will need to update all stan models

# if on own computer:
# setwd("~/Documents/Research/ZIBpaper/loZIB/ModelEstimates")

source("EstimateFunctions.R")

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


fit_result <- est_run(X, ll, y, visit, V=V, model=covmod, iter, warmup, chains,  seed=1234)

out.file <- file.path(paste0("draws_",covmod,".rds"))

print(covmod)

saveRDS(fit_result, out.file)



