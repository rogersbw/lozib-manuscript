# Exploring correlation posteriors

# Steps:
# Load posterior samples from AD model
# Extract correlation coefficients

# a1 = (sigma2[2] / sigma2[1]) *rho, so rho1 = a1 * sigma2[1]/sigma2[2]

# We want posterior distributions of rho[1], rho[2] and rho[3]
# Also sigma2[1], sigma2[2], ...
# AD Covariance Posterior summaries:


library(dplyr)

### Directory organization

#here::i_am("analysis/post_hoc/post_predict.r")

#proj_path <- here::here()

here::i_am("analysis/post_hoc/post_corr.R")

source("analysis/post_hoc/post_functions.r")

setwd("analysis/post_hoc")

post_draws <- extract_draws("heavy", "ad")


sigma2 <- post_draws$sigma2
rho <- post_draws$rho

sigma2_1 <- sigma2[, 1]
sigma2_2 <- sqrt(sigma2[,1]^2 + sigma2[,2]^2)
sigma2_3 <- sqrt(sigma2[,1]^2 + sigma2[,2]^2 + sigma2[,3]^2)
sigma2_4 <- sqrt(sigma2[,1]^2 + sigma2[,2]^2 + sigma2[,3]^2 + sigma2[,4]^2)
rho1 <-  rho[, 1] * sigma2_1 / sigma2_2
rho2 <-  rho[, 2] * sigma2_2 / sigma2_3
rho3 <-  rho[, 3] * sigma2_3 / sigma2_4


post_summary_df <- data.frame(
  sigma2_1 = sigma2[, 1],
  sigma2_2 = sqrt((rho[, 1] * sigma2[, 1])^2 + sigma2[, 2]^2),
  sigma2_3 = sqrt((rho[, 2] * sigma2_2)^2 + sigma2[, 3]^2),
  sigma2_4 = sqrt((rho[, 3] * sigma2_3)^2 + sigma2[, 4]^2),
  rho1 = rho[, 1] * sigma2_1 / sigma2_2,
  rho2 = rho[, 2] * sigma2_2 / sigma2_3,
  rho3 = rho[, 3] * sigma2_3 / sigma2_4
)

# To get the covariance posterior summaries
apply(post_summary_df, 2, function(x) {
  c(mean(x), quantile(x, c(0.025, 0.975)))
})

# Now for the UN model

