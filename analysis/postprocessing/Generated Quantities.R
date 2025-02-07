# Use stan output to create generated quantities
library(bayesplot)
library(RColorBrewer)
library(tidyverse)
library(forcats)

setwd(here::here())



###############################################################################
############## Utility Functions ##############################################
###############################################################################

logit<-function(x){log(x/(1-x))}
expit<-function(x){exp(x)/(1+exp(x))}


#Constructs an AR correlation matrix given a rho and a dimension
AR_mat <- function(rho, J){
  
  Omega <- matrix(0, nrow=J, ncol = J)
  for(i in 1:(J-1)){
    for(j in (i+1):J){
      Omega[i, j] <- rho^(j-i)
    }
  }
  Omega <- Omega + t(Omega) + diag(J)
  
  #diag(sigma2) %*% Omega %*% diag(sigma2)
  Omega
}


#Creates an AD matrix given a vector of correlations
AD_mat <- function(rho){
  J <- length(rho)+1
  
  Omega <- matrix(0, nrow=J, ncol = J)
  for(i in 1:(J-1)){
    for(j in (i+1):J){
      Omega[i, j] <- prod(rho[i:(j-1)])
    }
  }
  #Omega[1,2] <- rho[1]
  #Omega[1,3] <- rho[1]*rho[2]
  #Omega[1,4] <- rho[1]*rho[2]*rho[3]
  
  Omega <- Omega + t(Omega) + diag(J)
  
  #diag(sigma2) %*% Omega %*% diag(sigma2)
  Omega
}


## Fill Mean Estimate Arrays
fill_DoD <- function(i, j, idx){
  
  mu_est[i,j,,idx] <<- theta_est[i,j,,idx] * pi_est[i,j,,idx]
  
  DoD_pi[i,j,1,idx] <<- pi_est[i,j,5,idx] - pi_est[i,j,2,idx]
  DoD_pi[i,j,2,idx] <<- pi_est[i,j,6,idx] - pi_est[i,j,3,idx]
  DoD_pi[i,j,3,idx] <<- pi_est[i,j,7,idx] - pi_est[i,j,4,idx]
  DoD_pi[i,j,4,idx] <<- (DoD_pi[i,j,1,idx] + DoD_pi[i,j,2,idx] + DoD_pi[i,j,3,idx])/3
  
  DoD_theta[i,j,1,idx] <<- theta_est[i,j,5,idx] - theta_est[i,j,2,idx]
  DoD_theta[i,j,2,idx] <<- theta_est[i,j,6,idx] - theta_est[i,j,3,idx]
  DoD_theta[i,j,3,idx] <<- theta_est[i,j,7,idx] - theta_est[i,j,4,idx]
  DoD_theta[i,j,4,idx] <<- (DoD_theta[i,j,1,idx] + DoD_theta[i,j,2,idx] + DoD_theta[i,j,3,idx])/3
  
  DoD_mu[i,j,1,idx] <<- mu_est[i,j,5,idx] - mu_est[i,j,2,idx]
  DoD_mu[i,j,2,idx] <<- mu_est[i,j,6,idx] - mu_est[i,j,3,idx]
  DoD_mu[i,j,3,idx] <<- mu_est[i,j,7,idx] - mu_est[i,j,4,idx]
  DoD_mu[i,j,4,idx] <<- (DoD_mu[i,j,1,idx] + DoD_mu[i,j,2,idx] + DoD_mu[i,j,3,idx])/3
  
}


## Raw data values:

sbirt <- readRDS("sbirt_clean.rds")


mean_baseline <- sbirt |> filter(visit==1) |> pull(heavy) |> mean() 

mean_followup <- sbirt |> filter(visit>1) |> group_by(group, visit) |> summarize(mean = mean(heavy, na.rm=TRUE)) 

set.seed(1989)



#####Simulate standard normal for gamma1 and gamma2

gam1 <- rnorm(1000)
gam2 <- rnorm(1000)




###############################################################################
############.  1 shared RI Model ###############
###############################################################################
# Do the same as above for 1RI model

ri1_draws <- readRDS("Draws/draws_rifactor.rds")
str(ri1_draws)

dim(ri1_draws)


# Initialize arrays of model estimates 
pi_est <- array(data=NA, c(dim(ri1_draws)[1],dim(ri1_draws)[2],7, 12))
theta_est <- array(data=NA, c(dim(ri1_draws)[1],dim(ri1_draws)[2],7, 12))
mu_est <- array(data=NA, c(dim(ri1_draws)[1],dim(ri1_draws)[2],7, 12))

DoD_pi <- array(data=NA, c(dim(ri1_draws)[1],dim(ri1_draws)[2],4, 12))
DoD_theta <- array(data=NA, c(dim(ri1_draws)[1],dim(ri1_draws)[2],4, 12))
DoD_mu <- array(data=NA, c(dim(ri1_draws)[1],dim(ri1_draws)[2],4, 12))



#Sigma1
mean(ri1_draws[,,15], na.rm =T)
quantile(ri1_draws[,,15], probs=c(.025), na.rm=T)
quantile(ri1_draws[,,15], probs=c(.975), na.rm=T)


#Sigma2
mean(ri1_draws[,,16], na.rm =T)
quantile(ri1_draws[,,16], probs=c(.025), na.rm=T)
quantile(ri1_draws[,,16], probs=c(.975), na.rm=T)


for(i in 1:dim(ri1_draws)[1]){
  for(j in 1:dim(ri1_draws)[2]){
    
    # Zero model mean estimates for visit x treatment group  
    pi_est[i,j,1,1] <- mean(expit(as.vector(ri1_draws[i,j,1]) + as.numeric(ri1_draws[i,j,15])*gam1))
    pi_est[i,j,2,1] <- mean(expit(as.vector(ri1_draws[i,j,2]) + as.numeric(ri1_draws[i,j,15])*gam1))
    pi_est[i,j,3,1] <- mean(expit(as.vector(ri1_draws[i,j,3]) + as.numeric(ri1_draws[i,j,15])*gam1))
    pi_est[i,j,4,1] <- mean(expit(as.vector(ri1_draws[i,j,4]) + as.numeric(ri1_draws[i,j,15])*gam1))
    pi_est[i,j,5,1] <- mean(expit(as.vector(ri1_draws[i,j,2]) + as.vector(ri1_draws[i,j,5]) + as.numeric(ri1_draws[i,j,15])*gam1))
    pi_est[i,j,6,1] <- mean(expit(as.vector(ri1_draws[i,j,3]) + as.vector(ri1_draws[i,j,6]) + as.numeric(ri1_draws[i,j,15])*gam1))
    pi_est[i,j,7,1] <- mean(expit(as.vector(ri1_draws[i,j,4]) + as.vector(ri1_draws[i,j,7]) + as.numeric(ri1_draws[i,j,15])*gam1))
    
    
    #Count model mean estimates for visit x treatment group  
    
    constant <- as.numeric(ri1_draws[i,j,16])*gam1
    
    theta_est[i,j,1,1] <- 90*mean(expit(as.vector(ri1_draws[i,j,8]) +  constant))
    theta_est[i,j,2,1] <- 90*mean(expit(as.vector(ri1_draws[i,j,9]) +  constant))
    theta_est[i,j,3,1] <- 90*mean(expit(as.vector(ri1_draws[i,j,10]) +  constant))
    theta_est[i,j,4,1] <- 90*mean(expit(as.vector(ri1_draws[i,j,11]) + constant))
    theta_est[i,j,5,1] <- 90*mean(expit(as.vector(ri1_draws[i,j,9]) + as.vector(ri1_draws[i,j,12])  + constant))
    theta_est[i,j,6,1] <- 90*mean(expit(as.vector(ri1_draws[i,j,10]) + as.vector(ri1_draws[i,j,13])  + constant))
    theta_est[i,j,7,1] <- 90*mean(expit(as.vector(ri1_draws[i,j,11]) + as.vector(ri1_draws[i,j,14]) + constant))
    
    fill_DoD(i, j, 1)
  }
}



###############################################################################
############.  2 RI Model ###############
###############################################################################

ri_draws <- readRDS("Draws/draws_ri.rds")

str(ri_draws)

dimnames(ri_draws)[3]

#Psi posterior
mean(ri_draws[,,15], na.rm =T)
quantile(ri_draws[,,15], probs=c(.025), na.rm=T)
quantile(ri_draws[,,15], probs=c(.975), na.rm=T)

#Sigma1
mean(ri_draws[,,16], na.rm =T)
quantile(ri_draws[,,16], probs=c(.025), na.rm=T)
quantile(ri_draws[,,16], probs=c(.975), na.rm=T)


#Sigma2
mean(ri_draws[,,17], na.rm =T)
quantile(ri_draws[,,17], probs=c(.025), na.rm=T)
quantile(ri_draws[,,17], probs=c(.975), na.rm=T)


for(i in 1:dim(ri_draws)[1]){
  for(j in 1:dim(ri_draws)[2]){
     
     # Zero model mean estimates for visit x treatment group  
     pi_est[i,j,1,2] <- mean(expit(as.vector(ri_draws[i,j,1]) + as.numeric(ri_draws[i,j,16])*gam1))
     pi_est[i,j,2,2] <- mean(expit(as.vector(ri_draws[i,j,2]) + as.numeric(ri_draws[i,j,16])*gam1))
     pi_est[i,j,3,2] <- mean(expit(as.vector(ri_draws[i,j,3]) + as.numeric(ri_draws[i,j,16])*gam1))
     pi_est[i,j,4,2] <- mean(expit(as.vector(ri_draws[i,j,4]) + as.numeric(ri_draws[i,j,16])*gam1))
     pi_est[i,j,5,2] <- mean(expit(as.vector(ri_draws[i,j,2]) + as.vector(ri_draws[i,j,5]) + as.numeric(ri_draws[i,j,16])*gam1))
     pi_est[i,j,6,2] <- mean(expit(as.vector(ri_draws[i,j,3]) + as.vector(ri_draws[i,j,6]) + as.numeric(ri_draws[i,j,16])*gam1))
     pi_est[i,j,7,2] <- mean(expit(as.vector(ri_draws[i,j,4]) + as.vector(ri_draws[i,j,7]) + as.numeric(ri_draws[i,j,16])*gam1))
     
     
     #Count model mean estimates for visit x treatment group  
     
     constant <- as.numeric(ri_draws[i,j,15]) *  gam1 + as.numeric(ri_draws[i,j,17])*gam2
     
     theta_est[i,j,1,2] <- 90*mean(expit(as.vector(ri_draws[i,j,8]) +  constant))
     theta_est[i,j,2,2] <- 90*mean(expit(as.vector(ri_draws[i,j,9]) +  constant))
     theta_est[i,j,3,2] <- 90*mean(expit(as.vector(ri_draws[i,j,10]) +  constant))
     theta_est[i,j,4,2] <- 90*mean(expit(as.vector(ri_draws[i,j,11]) + constant))
     theta_est[i,j,5,2] <- 90*mean(expit(as.vector(ri_draws[i,j,9]) + as.vector(ri_draws[i,j,12])  + constant))
     theta_est[i,j,6,2] <- 90*mean(expit(as.vector(ri_draws[i,j,10]) + as.vector(ri_draws[i,j,13])  + constant))
     theta_est[i,j,7,2] <- 90*mean(expit(as.vector(ri_draws[i,j,11]) + as.vector(ri_draws[i,j,14]) + constant))
     
     fill_DoD(i, j, 2)
  }
}


###############################################################################
############                OLRE Models                         ###############
###############################################################################

# Sample gamma2 matrix for Observation Level Random Effects (OLRE)
gam2 <- matrix(rnorm(4000), ncol=4)


###############################################################################
############                  INDcv Model                       ###############
###############################################################################

indcv_draws <- readRDS("Draws/draws_indcv.rds")

dimnames(indcv_draws)[3]

#Psi posterior
mean(indcv_draws[,,15], na.rm =T)
quantile(indcv_draws[,,15], probs=c(.025), na.rm=T)
quantile(indcv_draws[,,15], probs=c(.975), na.rm=T)

#Sigma1
mean(indcv_draws[,,16], na.rm =T)
quantile(indcv_draws[,,16], probs=c(.025), na.rm=T)
quantile(indcv_draws[,,16], probs=c(.975), na.rm=T)


#Sigma2
mean(indcv_draws[,,17], na.rm =T)
quantile(indcv_draws[,,17], probs=c(.025), na.rm=T)
quantile(indcv_draws[,,17], probs=c(.975), na.rm=T)


for(i in 1:dim(indcv_draws)[1]){
  for(j in 1:dim(indcv_draws)[2]){
    
    constant <- as.numeric(indcv_draws[i,j,16])*gam1
    
    # Zero model mean estimates for visit x treatment group  
    pi_est[i,j,1,3] <- mean(expit(as.vector(indcv_draws[i,j,1]) + constant))
    pi_est[i,j,2,3] <- mean(expit(as.vector(indcv_draws[i,j,2]) + constant))
    pi_est[i,j,3,3] <- mean(expit(as.vector(indcv_draws[i,j,3]) + constant))
    pi_est[i,j,4,3] <- mean(expit(as.vector(indcv_draws[i,j,4]) + constant))
    pi_est[i,j,5,3] <- mean(expit(as.vector(indcv_draws[i,j,2]) + as.vector(indcv_draws[i,j,5]) + constant))
    pi_est[i,j,6,3] <- mean(expit(as.vector(indcv_draws[i,j,3]) + as.vector(indcv_draws[i,j,6]) + constant))
    pi_est[i,j,7,3] <- mean(expit(as.vector(indcv_draws[i,j,4]) + as.vector(indcv_draws[i,j,7]) + constant))
    
    
    #Count model mean estimates for visit x treatment group  
    
    constant <- as.numeric(indcv_draws[i,j,15]) *  gam1
    gam2_scaled <- as.numeric(indcv_draws[i,j,17]) * gam2
    
    theta_est[i,j,1,3] <- 90*mean(expit(as.vector(indcv_draws[i,j,8]) +  constant + gam2_scaled[,1]))
    theta_est[i,j,2,3] <- 90*mean(expit(as.vector(indcv_draws[i,j,9]) +  constant  + gam2_scaled[,2]))
    theta_est[i,j,3,3] <- 90*mean(expit(as.vector(indcv_draws[i,j,10]) +  constant  + gam2_scaled[,3]))
    theta_est[i,j,4,3] <- 90*mean(expit(as.vector(indcv_draws[i,j,11]) + constant +  + gam2_scaled[,4]))
    theta_est[i,j,5,3] <- 90*mean(expit(as.vector(indcv_draws[i,j,9]) + as.vector(indcv_draws[i,j,12])  + constant + gam2_scaled[,2]))
    theta_est[i,j,6,3] <- 90*mean(expit(as.vector(indcv_draws[i,j,10]) + as.vector(indcv_draws[i,j,13])  + constant + gam2_scaled[,3]))
    theta_est[i,j,7,3] <- 90*mean(expit(as.vector(indcv_draws[i,j,11]) + as.vector(indcv_draws[i,j,14]) + constant + gam2_scaled[,4]))
    
    fill_DoD(i, j, 3)
  }
}




###############################################################################
############                  IND Model                       ###############
###############################################################################

ind_draws <- readRDS("Draws/draws_ind.rds")

dimnames(ind_draws)[3]

#Psi posterior
mean(ind_draws[,,15], na.rm =T)
quantile(ind_draws[,,15], probs=c(.025), na.rm=T)
quantile(ind_draws[,,15], probs=c(.975), na.rm=T)

#Sigma1
mean(ind_draws[,,16], na.rm =T)
quantile(ind_draws[,,16], probs=c(.025), na.rm=T)
quantile(ind_draws[,,16], probs=c(.975), na.rm=T)


for(i in 1:dim(ind_draws)[1]){
  for(j in 1:dim(ind_draws)[2]){
    
    constant <- as.numeric(ind_draws[i,j,16])*gam1
    
    # Zero model mean estimates for visit x treatment group  
    pi_est[i,j,1,4] <- mean(expit(as.vector(ind_draws[i,j,1]) + constant))
    pi_est[i,j,2,4] <- mean(expit(as.vector(ind_draws[i,j,2]) + constant))
    pi_est[i,j,3,4] <- mean(expit(as.vector(ind_draws[i,j,3]) + constant))
    pi_est[i,j,4,4] <- mean(expit(as.vector(ind_draws[i,j,4]) + constant))
    pi_est[i,j,5,4] <- mean(expit(as.vector(ind_draws[i,j,2]) + as.vector(ind_draws[i,j,5]) + constant))
    pi_est[i,j,6,4] <- mean(expit(as.vector(ind_draws[i,j,3]) + as.vector(ind_draws[i,j,6]) + constant))
    pi_est[i,j,7,4] <- mean(expit(as.vector(ind_draws[i,j,4]) + as.vector(ind_draws[i,j,7]) + constant))
    
    
    #Count model mean estimates for visit x treatment group  
    
    constant <- as.numeric(ind_draws[i,j,15]) *  gam1
    gam2_scaled <- as.numeric(ind_draws[i,j,17]) * gam2
    
    theta_est[i,j,1,4] <- 90*mean(expit(as.vector(ind_draws[i,j,8]) +  constant + gam2_scaled[,1]))
    theta_est[i,j,2,4] <- 90*mean(expit(as.vector(ind_draws[i,j,9]) +  constant  + gam2_scaled[,2]))
    theta_est[i,j,3,4] <- 90*mean(expit(as.vector(ind_draws[i,j,10]) +  constant  + gam2_scaled[,3]))
    theta_est[i,j,4,4] <- 90*mean(expit(as.vector(ind_draws[i,j,11]) + constant +  + gam2_scaled[,4]))
    theta_est[i,j,5,4] <- 90*mean(expit(as.vector(ind_draws[i,j,9]) + as.vector(ind_draws[i,j,12])  + constant + gam2_scaled[,2]))
    theta_est[i,j,6,4] <- 90*mean(expit(as.vector(ind_draws[i,j,10]) + as.vector(ind_draws[i,j,13])  + constant + gam2_scaled[,4]))
    theta_est[i,j,7,4] <- 90*mean(expit(as.vector(ind_draws[i,j,11]) + as.vector(ind_draws[i,j,14]) + constant + gam2_scaled[,4]))
    
    fill_DoD(i, j, 4)
  }
}





###############################################################################
############                  CScv Model                          ###############
###############################################################################

cscv_draws <- readRDS("Draws/draws_cscv.rds")

dimnames(cscv_draws)[3]

#sigma1 posterior
mean(cscv_draws[,,15], na.rm =T)
quantile(cscv_draws[,,15], probs=c(.025), na.rm=T)
quantile(cscv_draws[,,15], probs=c(.975), na.rm=T)

#Sigma2
mean(cscv_draws[,,16], na.rm =T)
quantile(cscv_draws[,,16], probs=c(.025), na.rm=T)
quantile(cscv_draws[,,16], probs=c(.975), na.rm=T)


#Psi
mean(cscv_draws[,,17], na.rm =T)
quantile(cscv_draws[,,17], probs=c(.025), na.rm=T)
quantile(cscv_draws[,,17], probs=c(.975), na.rm=T)

## Need an RI element to construct CS model
gam2_ri <- rnorm(1000)

for(i in 1:dim(cscv_draws)[1]){
  for(j in 1:dim(cscv_draws)[2]){
    
    constant <- as.numeric(cscv_draws[i,j,16])*gam1
    
    # Zero model mean estimates for visit x treatment group  
    pi_est[i,j,1,5] <- mean(expit(as.vector(cscv_draws[i,j,1]) + constant))
    pi_est[i,j,2,5] <- mean(expit(as.vector(cscv_draws[i,j,2]) + constant))
    pi_est[i,j,3,5] <- mean(expit(as.vector(cscv_draws[i,j,3]) + constant))
    pi_est[i,j,4,5] <- mean(expit(as.vector(cscv_draws[i,j,4]) + constant))
    pi_est[i,j,5,5] <- mean(expit(as.vector(cscv_draws[i,j,2]) + as.vector(cscv_draws[i,j,5]) + constant))
    pi_est[i,j,6,5] <- mean(expit(as.vector(cscv_draws[i,j,3]) + as.vector(cscv_draws[i,j,6]) + constant))
    pi_est[i,j,7,5] <- mean(expit(as.vector(cscv_draws[i,j,4]) + as.vector(cscv_draws[i,j,7]) + constant))
    
    
    #Count model mean estimates for visit x treatment group  
    
    constant <- as.numeric(cscv_draws[i,j,17]) *  gam1
    gam2_scaled <-  as.numeric(cscv_draws[i,j,18]) * gam2_ri + as.numeric(cscv_draws[i,j,16]) * gam2
    
    theta_est[i,j,1,5] <- 90*mean(expit(as.vector(cscv_draws[i,j,8]) +  constant + gam2_scaled[,1]))
    theta_est[i,j,2,5] <- 90*mean(expit(as.vector(cscv_draws[i,j,9]) +  constant  + gam2_scaled[,2]))
    theta_est[i,j,3,5] <- 90*mean(expit(as.vector(cscv_draws[i,j,10]) +  constant  + gam2_scaled[,3]))
    theta_est[i,j,4,5] <- 90*mean(expit(as.vector(cscv_draws[i,j,11]) + constant + gam2_scaled[,4]))
    theta_est[i,j,5,5] <- 90*mean(expit(as.vector(cscv_draws[i,j,9]) + as.vector(cscv_draws[i,j,12])  + constant + gam2_scaled[,2]))
    theta_est[i,j,6,5] <- 90*mean(expit(as.vector(cscv_draws[i,j,10]) + as.vector(cscv_draws[i,j,13])  + constant + gam2_scaled[,4]))
    theta_est[i,j,7,5] <- 90*mean(expit(as.vector(cscv_draws[i,j,11]) + as.vector(cscv_draws[i,j,14]) + constant + gam2_scaled[,4]))
    
    fill_DoD(i, j, 5)
  }
}



###############################################################################
############                  CS Model                          ###############
###############################################################################

cs_draws <- readRDS("Draws/draws_cs.rds")

dimnames(cs_draws)[3]

#Sigma1 posterior
mean(cs_draws[,,15], na.rm =T)
quantile(cs_draws[,,15], probs=c(.025), na.rm=T)
quantile(cs_draws[,,15], probs=c(.975), na.rm=T)

#psi
mean(cs_draws[,,20], na.rm =T)
quantile(cs_draws[,,20], probs=c(.025), na.rm=T)
quantile(cs_draws[,,20], probs=c(.975), na.rm=T)


## Need an RI element to construct CS model
gam2_ri <- rnorm(1000)

for(i in 1:dim(cs_draws)[1]){
  for(j in 1:dim(cs_draws)[2]){
    
    constant <- as.numeric(cs_draws[i,j,16])*gam1
    
    # Zero model mean estimates for visit x treatment group  
    pi_est[i,j,1,6] <- mean(expit(as.vector(cs_draws[i,j,1]) + constant))
    pi_est[i,j,2,6] <- mean(expit(as.vector(cs_draws[i,j,2]) + constant))
    pi_est[i,j,3,6] <- mean(expit(as.vector(cs_draws[i,j,3]) + constant))
    pi_est[i,j,4,6] <- mean(expit(as.vector(cs_draws[i,j,4]) + constant))
    pi_est[i,j,5,6] <- mean(expit(as.vector(cs_draws[i,j,2]) + as.vector(cs_draws[i,j,5]) + constant))
    pi_est[i,j,6,6] <- mean(expit(as.vector(cs_draws[i,j,3]) + as.vector(cs_draws[i,j,6]) + constant))
    pi_est[i,j,7,6] <- mean(expit(as.vector(cs_draws[i,j,4]) + as.vector(cs_draws[i,j,7]) + constant))
    
    
    #Count model mean estimates for visit x treatment group  
    
    constant <- as.numeric(cs_draws[i,j,20]) *  gam1
    gam2_scaled <-  as.numeric(cs_draws[i,j,21]) * gam2_ri + gam2 %*% diag(as.numeric(cs_draws[i,j,16:19]))
    
    theta_est[i,j,1,6] <- 90*mean(expit(as.vector(cs_draws[i,j,8]) +  constant + gam2_scaled[,1]))
    theta_est[i,j,2,6] <- 90*mean(expit(as.vector(cs_draws[i,j,9]) +  constant  + gam2_scaled[,2]))
    theta_est[i,j,3,6] <- 90*mean(expit(as.vector(cs_draws[i,j,10]) +  constant  + gam2_scaled[,3]))
    theta_est[i,j,4,6] <- 90*mean(expit(as.vector(cs_draws[i,j,11]) + constant + gam2_scaled[,4]))
    theta_est[i,j,5,6] <- 90*mean(expit(as.vector(cs_draws[i,j,9]) + as.vector(cs_draws[i,j,12])  + constant + gam2_scaled[,2]))
    theta_est[i,j,6,6] <- 90*mean(expit(as.vector(cs_draws[i,j,10]) + as.vector(cs_draws[i,j,13])  + constant + gam2_scaled[,4]))
    theta_est[i,j,7,6] <- 90*mean(expit(as.vector(cs_draws[i,j,11]) + as.vector(cs_draws[i,j,14]) + constant + gam2_scaled[,4]))
    
    fill_DoD(i, j, 6)
  }
}


###############################################################################
############                  ARcv Model                          ###############
###############################################################################


arcv_draws <- readRDS("Draws/draws_arcv.rds")

dimnames(arcv_draws)[3]

#Rho posterior
mean(arcv_draws[,,18], na.rm =T)
quantile(arcv_draws[,,18], probs=c(.025), na.rm=T)
quantile(arcv_draws[,,18], probs=c(.975), na.rm=T)

#Psi posterior
mean(arcv_draws[,,15], na.rm =T)
quantile(arcv_draws[,,15], probs=c(.025), na.rm=T)
quantile(arcv_draws[,,15], probs=c(.975), na.rm=T)

#Sigma1
mean(arcv_draws[,,16], na.rm =T)
quantile(arcv_draws[,,16], probs=c(.025), na.rm=T)
quantile(arcv_draws[,,16], probs=c(.975), na.rm=T)


#Sigma2
mean(arcv_draws[,,17], na.rm =T)
quantile(arcv_draws[,,17], probs=c(.025), na.rm=T)
quantile(arcv_draws[,,17], probs=c(.975), na.rm=T)





#apply(arcv_draws[,,18], 3, mean, na.rm=T)
#apply(arcv_draws[,,18], 3, function(x) quantile(x, probs=c(.025), na.rm=T))
#apply(arcv_draws[,,18], 3, function(x) quantile(x, probs=c(.975), na.rm=T))

 
for(i in 1:dim(arcv_draws)[1]){
  for(j in 1:dim(arcv_draws)[2]){
    
    chol_Sig <- chol(AR_mat(as.numeric(arcv_draws[i,j,18]), 4))
    
    gam2_scaled <- as.numeric(arcv_draws[i,j,17]) * gam2 %*% chol_Sig 
    
    # Zero model mean estimates for visit x treatment group  
    pi_est[i,j,1,7] <- mean(expit(as.vector(arcv_draws[i,j,1]) + as.numeric(arcv_draws[i,j,16])*gam1))
    pi_est[i,j,2,7] <- mean(expit(as.vector(arcv_draws[i,j,2]) + as.numeric(arcv_draws[i,j,16])*gam1))
    pi_est[i,j,3,7] <- mean(expit(as.vector(arcv_draws[i,j,3]) + as.numeric(arcv_draws[i,j,16])*gam1))
    pi_est[i,j,4,7] <- mean(expit(as.vector(arcv_draws[i,j,4]) + as.numeric(arcv_draws[i,j,16])*gam1))
    pi_est[i,j,5,7] <- mean(expit(as.vector(arcv_draws[i,j,2]) + as.vector(arcv_draws[i,j,5]) + as.numeric(arcv_draws[i,j,16])*gam1))
    pi_est[i,j,6,7] <- mean(expit(as.vector(arcv_draws[i,j,3]) + as.vector(arcv_draws[i,j,6]) + as.numeric(arcv_draws[i,j,16])*gam1))
    pi_est[i,j,7,7] <- mean(expit(as.vector(arcv_draws[i,j,4]) + as.vector(arcv_draws[i,j,7]) + as.numeric(arcv_draws[i,j,16])*gam1))
    
    
    #Count model mean estimates for visit x treatment group  
    
    constant <- as.numeric(arcv_draws[i,j,15]) *  gam1 
    
    theta_est[i,j,1,7] <- 90*mean(expit(as.vector(arcv_draws[i,j,8]) +  constant + gam2_scaled[,1]))
    theta_est[i,j,2,7] <- 90*mean(expit(as.vector(arcv_draws[i,j,9]) +  constant + gam2_scaled[,2]))
    theta_est[i,j,3,7] <- 90*mean(expit(as.vector(arcv_draws[i,j,10]) +  constant + gam2_scaled[,3]))
    theta_est[i,j,4,7] <- 90*mean(expit(as.vector(arcv_draws[i,j,11]) + constant + gam2_scaled[,4]))
    theta_est[i,j,5,7] <- 90*mean(expit(as.vector(arcv_draws[i,j,9]) + as.vector(arcv_draws[i,j,12] )  + constant + gam2_scaled[,2]))
    theta_est[i,j,6,7] <- 90*mean(expit(as.vector(arcv_draws[i,j,10]) + as.vector(arcv_draws[i,j,13])  + constant + gam2_scaled[,3]))
    theta_est[i,j,7,7] <- 90*mean(expit(as.vector(arcv_draws[i,j,11]) + as.vector(arcv_draws[i,j,14]) + constant + gam2_scaled[,4]))
    
    
    fill_DoD(i, j, 7)
  }
}




###############################################################################
############                  AR Model                          ###############
###############################################################################

ar_draws <- readRDS("Draws/draws_ar.rds")

for(i in 1:dim(ar_draws)[1]){
  for(j in 1:dim(ar_draws)[2]){
    
    chol_Sig <- chol(AR_mat(as.numeric(ar_draws[i,j,21]), 4))
    
    gam2_scaled <-  gam2 %*% chol_Sig %*% diag(as.numeric(ar_draws[i,j,17:20]))
    
    # Zero model mean estimates for visit x treatment group  
    pi_est[i,j,1,8] <- mean(expit(as.vector(ar_draws[i,j,1]) + as.numeric(ar_draws[i,j,16])*gam1))
    pi_est[i,j,2,8] <- mean(expit(as.vector(ar_draws[i,j,2]) + as.numeric(ar_draws[i,j,16])*gam1))
    pi_est[i,j,3,8] <- mean(expit(as.vector(ar_draws[i,j,3]) + as.numeric(ar_draws[i,j,16])*gam1))
    pi_est[i,j,4,8] <- mean(expit(as.vector(ar_draws[i,j,4]) + as.numeric(ar_draws[i,j,16])*gam1))
    pi_est[i,j,5,8] <- mean(expit(as.vector(ar_draws[i,j,2]) + as.vector(ar_draws[i,j,5]) + as.numeric(ar_draws[i,j,16])*gam1))
    pi_est[i,j,6,8] <- mean(expit(as.vector(ar_draws[i,j,3]) + as.vector(ar_draws[i,j,6]) + as.numeric(ar_draws[i,j,16])*gam1))
    pi_est[i,j,7,8] <- mean(expit(as.vector(ar_draws[i,j,4]) + as.vector(ar_draws[i,j,7]) + as.numeric(ar_draws[i,j,16])*gam1))
    
    
    #Count model mean estimates for visit x treatment group  
    
    constant <- as.numeric(ar_draws[i,j,15]) *  gam1 
    
    theta_est[i,j,1,8] <- 90*mean(expit(as.vector(ar_draws[i,j,8]) +  constant + gam2_scaled[,1]))
    theta_est[i,j,2,8] <- 90*mean(expit(as.vector(ar_draws[i,j,9]) +  constant + gam2_scaled[,2]))
    theta_est[i,j,3,8] <- 90*mean(expit(as.vector(ar_draws[i,j,10]) +  constant + gam2_scaled[,3]))
    theta_est[i,j,4,8] <- 90*mean(expit(as.vector(ar_draws[i,j,11]) + constant + gam2_scaled[,4]))
    theta_est[i,j,5,8] <- 90*mean(expit(as.vector(ar_draws[i,j,9]) + as.vector(ar_draws[i,j,12] )  + constant + gam2_scaled[,2]))
    theta_est[i,j,6,8] <- 90*mean(expit(as.vector(ar_draws[i,j,10]) + as.vector(ar_draws[i,j,13])  + constant + gam2_scaled[,3]))
    theta_est[i,j,7,8] <- 90*mean(expit(as.vector(ar_draws[i,j,11]) + as.vector(ar_draws[i,j,14]) + constant + gam2_scaled[,4]))
    
    
    fill_DoD(i, j, 8)
  }
}





###############################################################################
############                  ADcv Model                          ###############
###############################################################################

adcv_draws <- readRDS("Draws/draws_adcv.rds")

for(i in 1:dim(adcv_draws)[1]){
  for(j in 1:dim(adcv_draws)[2]){
    
    chol_Sig <- chol(AD_mat(as.numeric(adcv_draws[i,j,18:20])))
    
    gam2_scaled <- as.numeric(adcv_draws[i,j,17]) * gam2 %*% chol_Sig 
    
    # Zero model mean estimates for visit x treatment group  
    pi_est[i,j,1,9] <- mean(expit(as.vector(adcv_draws[i,j,1]) + as.numeric(adcv_draws[i,j,16])*gam1))
    pi_est[i,j,2,9] <- mean(expit(as.vector(adcv_draws[i,j,2]) + as.numeric(adcv_draws[i,j,16])*gam1))
    pi_est[i,j,3,9] <- mean(expit(as.vector(adcv_draws[i,j,3]) + as.numeric(adcv_draws[i,j,16])*gam1))
    pi_est[i,j,4,9] <- mean(expit(as.vector(adcv_draws[i,j,4]) + as.numeric(adcv_draws[i,j,16])*gam1))
    pi_est[i,j,5,9] <- mean(expit(as.vector(adcv_draws[i,j,2]) + as.vector(adcv_draws[i,j,5]) + as.numeric(adcv_draws[i,j,16])*gam1))
    pi_est[i,j,6,9] <- mean(expit(as.vector(adcv_draws[i,j,3]) + as.vector(adcv_draws[i,j,6]) + as.numeric(adcv_draws[i,j,16])*gam1))
    pi_est[i,j,7,9] <- mean(expit(as.vector(adcv_draws[i,j,4]) + as.vector(adcv_draws[i,j,7]) + as.numeric(adcv_draws[i,j,16])*gam1))
    
    
    #Count model mean estimates for visit x treatment group  
    
    constant <- as.numeric(adcv_draws[i,j,15]) *  gam1 
    
    theta_est[i,j,1,9] <- 90*mean(expit(as.vector(adcv_draws[i,j,8]) +  constant + gam2_scaled[,1]))
    theta_est[i,j,2,9] <- 90*mean(expit(as.vector(adcv_draws[i,j,9]) +  constant + gam2_scaled[,2]))
    theta_est[i,j,3,9] <- 90*mean(expit(as.vector(adcv_draws[i,j,10]) +  constant + gam2_scaled[,3]))
    theta_est[i,j,4,9] <- 90*mean(expit(as.vector(adcv_draws[i,j,11]) + constant + gam2_scaled[,4]))
    theta_est[i,j,5,9] <- 90*mean(expit(as.vector(adcv_draws[i,j,9]) + as.vector(adcv_draws[i,j,12] )  + constant + gam2_scaled[,2]))
    theta_est[i,j,6,9] <- 90*mean(expit(as.vector(adcv_draws[i,j,10]) + as.vector(adcv_draws[i,j,13])  + constant + gam2_scaled[,3]))
    theta_est[i,j,7,9] <- 90*mean(expit(as.vector(adcv_draws[i,j,11]) + as.vector(adcv_draws[i,j,14]) + constant + gam2_scaled[,4]))
    
    
    fill_DoD(i, j, 9)
    
  }
  
  
}



###############################################################################
############                  AD Model                          ###############
###############################################################################


ad_draws <- readRDS("Draws/draws_ad.rds")

for(i in 1:dim(ad_draws)[1]){
  for(j in 1:dim(ad_draws)[2]){
    
    chol_Sig <- chol(AD_mat(as.numeric(ad_draws[i,j,21:23])))
    
    gam2_scaled <-   gam2 %*% chol_Sig %*% diag(as.numeric(ad_draws[i,j,17:20]))
    
    # Zero model mean estimates for visit x treatment group  
    pi_est[i,j,1,10] <- mean(expit(as.vector(ad_draws[i,j,1]) + as.numeric(ad_draws[i,j,16])*gam1))
    pi_est[i,j,2,10] <- mean(expit(as.vector(ad_draws[i,j,2]) + as.numeric(ad_draws[i,j,16])*gam1))
    pi_est[i,j,3,10] <- mean(expit(as.vector(ad_draws[i,j,3]) + as.numeric(ad_draws[i,j,16])*gam1))
    pi_est[i,j,4,10] <- mean(expit(as.vector(ad_draws[i,j,4]) + as.numeric(ad_draws[i,j,16])*gam1))
    pi_est[i,j,5,10] <- mean(expit(as.vector(ad_draws[i,j,2]) + as.vector(ad_draws[i,j,5]) + as.numeric(ad_draws[i,j,16])*gam1))
    pi_est[i,j,6,10] <- mean(expit(as.vector(ad_draws[i,j,3]) + as.vector(ad_draws[i,j,6]) + as.numeric(ad_draws[i,j,16])*gam1))
    pi_est[i,j,7,10] <- mean(expit(as.vector(ad_draws[i,j,4]) + as.vector(ad_draws[i,j,7]) + as.numeric(ad_draws[i,j,16])*gam1))
    
    
    #Count model mean estimates for visit x treatment group  
    
    constant <- as.numeric(ad_draws[i,j,15]) *  gam1 
    
    theta_est[i,j,1,10] <- 90*mean(expit(as.vector(ad_draws[i,j,8]) +  constant + gam2_scaled[,1]))
    theta_est[i,j,2,10] <- 90*mean(expit(as.vector(ad_draws[i,j,9]) +  constant + gam2_scaled[,2]))
    theta_est[i,j,3,10] <- 90*mean(expit(as.vector(ad_draws[i,j,10]) +  constant + gam2_scaled[,3]))
    theta_est[i,j,4,10] <- 90*mean(expit(as.vector(ad_draws[i,j,11]) + constant + gam2_scaled[,4]))
    theta_est[i,j,5,10] <- 90*mean(expit(as.vector(ad_draws[i,j,9]) + as.vector(ad_draws[i,j,12] )  + constant + gam2_scaled[,2]))
    theta_est[i,j,6,10] <- 90*mean(expit(as.vector(ad_draws[i,j,10]) + as.vector(ad_draws[i,j,13])  + constant + gam2_scaled[,3]))
    theta_est[i,j,7,10] <- 90*mean(expit(as.vector(ad_draws[i,j,11]) + as.vector(ad_draws[i,j,14]) + constant + gam2_scaled[,4]))
    
    
    fill_DoD(i, j, 10)
  }
}




###############################################################################
############                UNcv Model                          ###############
###############################################################################


uncv_draws <- readRDS("Draws/draws_uncv.rds")

for(i in 1:dim(uncv_draws)[1]){
  for(j in 1:dim(uncv_draws)[2]){
    
    chol_Sig <- chol(matrix(as.numeric(uncv_draws[i,j,18:33]), nrow=4, byrow=F))
    
    gam2_scaled <-  as.numeric(uncv_draws[i,j,17]) * gam2 %*% chol_Sig #%*% diag(as.numeric(uncv_draws[i,j,17:20]))
    
    # Zero model mean estimates for visit x treatment group  
    pi_est[i,j,1,11] <- mean(expit(as.vector(uncv_draws[i,j,1]) + as.numeric(uncv_draws[i,j,16])*gam1))
    pi_est[i,j,2,11] <- mean(expit(as.vector(uncv_draws[i,j,2]) + as.numeric(uncv_draws[i,j,16])*gam1))
    pi_est[i,j,3,11] <- mean(expit(as.vector(uncv_draws[i,j,3]) + as.numeric(uncv_draws[i,j,16])*gam1))
    pi_est[i,j,4,11] <- mean(expit(as.vector(uncv_draws[i,j,4]) + as.numeric(uncv_draws[i,j,16])*gam1))
    pi_est[i,j,5,11] <- mean(expit(as.vector(uncv_draws[i,j,2]) + as.vector(uncv_draws[i,j,5]) + as.numeric(uncv_draws[i,j,16])*gam1))
    pi_est[i,j,6,11] <- mean(expit(as.vector(uncv_draws[i,j,3]) + as.vector(uncv_draws[i,j,6]) + as.numeric(uncv_draws[i,j,16])*gam1))
    pi_est[i,j,7,11] <- mean(expit(as.vector(uncv_draws[i,j,4]) + as.vector(uncv_draws[i,j,7]) + as.numeric(uncv_draws[i,j,16])*gam1))
    
    
    #Count model mean estimates for visit x treatment group  
    
    constant <- as.numeric(uncv_draws[i,j,15]) *  gam1 
    
    theta_est[i,j,1,11] <- 90*mean(expit(as.vector(uncv_draws[i,j,8]) +  constant + gam2_scaled[,1]))
    theta_est[i,j,2,11] <- 90*mean(expit(as.vector(uncv_draws[i,j,9]) +  constant + gam2_scaled[,2]))
    theta_est[i,j,3,11] <- 90*mean(expit(as.vector(uncv_draws[i,j,10]) +  constant + gam2_scaled[,3]))
    theta_est[i,j,4,11] <- 90*mean(expit(as.vector(uncv_draws[i,j,11]) + constant + gam2_scaled[,4]))
    theta_est[i,j,5,11] <- 90*mean(expit(as.vector(uncv_draws[i,j,9]) + as.vector(uncv_draws[i,j,12] )  + constant + gam2_scaled[,2]))
    theta_est[i,j,6,11] <- 90*mean(expit(as.vector(uncv_draws[i,j,10]) + as.vector(uncv_draws[i,j,13])  + constant + gam2_scaled[,3]))
    theta_est[i,j,7,11] <- 90*mean(expit(as.vector(uncv_draws[i,j,11]) + as.vector(uncv_draws[i,j,14]) + constant + gam2_scaled[,4]))
     
    
    fill_DoD(i, j, 11)
    
  }
}



###############################################################################
############                 UN Model                          ###############
###############################################################################

un_draws <- readRDS("Draws/draws_un.rds")

for(i in 1:dim(un_draws)[1]){
  for(j in 1:dim(un_draws)[2]){
    
    chol_Sig <- chol(matrix(as.numeric(un_draws[i,j,21:36]), nrow=4, byrow=F))
    
    gam2_scaled <- gam2 %*% chol_Sig %*% diag(as.numeric(un_draws[i,j,17:20]))
    
    # Zero model mean estimates for visit x treatment group  
    pi_est[i,j,1,12] <- mean(expit(as.vector(un_draws[i,j,1]) + as.numeric(un_draws[i,j,16])*gam1))
    pi_est[i,j,2,12] <- mean(expit(as.vector(un_draws[i,j,2]) + as.numeric(un_draws[i,j,16])*gam1))
    pi_est[i,j,3,12] <- mean(expit(as.vector(un_draws[i,j,3]) + as.numeric(un_draws[i,j,16])*gam1))
    pi_est[i,j,4,12] <- mean(expit(as.vector(un_draws[i,j,4]) + as.numeric(un_draws[i,j,16])*gam1))
    pi_est[i,j,5,12] <- mean(expit(as.vector(un_draws[i,j,2]) + as.vector(un_draws[i,j,5]) + as.numeric(un_draws[i,j,16])*gam1))
    pi_est[i,j,6,12] <- mean(expit(as.vector(un_draws[i,j,3]) + as.vector(un_draws[i,j,6]) + as.numeric(un_draws[i,j,16])*gam1))
    pi_est[i,j,7,12] <- mean(expit(as.vector(un_draws[i,j,4]) + as.vector(un_draws[i,j,7]) + as.numeric(un_draws[i,j,16])*gam1))
    
    
    #Count model mean estimates for visit x treatment group  
    
    constant <- as.numeric(un_draws[i,j,15]) *  gam1 
    
    theta_est[i,j,1,12] <- 90*mean(expit(as.vector(un_draws[i,j,8]) +  constant + gam2_scaled[,1]))
    theta_est[i,j,2,12] <- 90*mean(expit(as.vector(un_draws[i,j,9]) +  constant + gam2_scaled[,2]))
    theta_est[i,j,3,12] <- 90*mean(expit(as.vector(un_draws[i,j,10]) +  constant + gam2_scaled[,3]))
    theta_est[i,j,4,12] <- 90*mean(expit(as.vector(un_draws[i,j,11]) + constant + gam2_scaled[,4]))
    theta_est[i,j,5,12] <- 90*mean(expit(as.vector(un_draws[i,j,9]) + as.vector(un_draws[i,j,12] )  + constant + gam2_scaled[,2]))
    theta_est[i,j,6,12] <- 90*mean(expit(as.vector(un_draws[i,j,10]) + as.vector(un_draws[i,j,13])  + constant + gam2_scaled[,3]))
    theta_est[i,j,7,12] <- 90*mean(expit(as.vector(un_draws[i,j,11]) + as.vector(un_draws[i,j,14]) + constant + gam2_scaled[,4]))
    
    fill_DoD(i, j, 12)    

  }
}


#Summary of Omega2

Omega_est_post <- array(dim = c(4,16))
dimnames(Omega_est_post)[[1]] <-  c("mean", "sd", "LB", "UB")

Omega_est_post[1,] <- apply(un_draws[,,21:36], 3, mean, na.rm=T)
Omega_est_post[2,]  <- apply(un_draws[,,21:36], 3, sd, na.rm=T)
Omega_est_post[3,]  <- apply(un_draws[,,21:36], 3, function(x) quantile(x, probs=c(.025), na.rm=T))
Omega_est_post[4,]  <- apply(un_draws[,,21:36], 3, function(x) quantile(x, probs=c(.975), na.rm=T))

matrix(Omega_est_post[1,], nrow=4, byrow=F)
matrix(Omega_est_post[3,], nrow=4, byrow=F)
matrix(Omega_est_post[4,], nrow=4, byrow=F)


#### Plot of Omega parameters posterior distributions
  as_tibble(reshape2::melt(un_draws[,,21:36])) %>% mutate(chain=as.factor(chain)) %>% 
    ggplot(aes(x=value, fill=chain, color=chain)) +
    geom_density() +
      facet_wrap(~variable)


#Summary of Omega2cv

Omega_est_postcv <- array(dim = c(4,16))
dimnames(Omega_est_postcv)[[1]] <-  c("mean", "sd", "LB", "UB")

Omega_est_postcv[1,] <- apply(uncv_draws[,,18:33], 3, mean, na.rm=T)
Omega_est_postcv[2,]  <- apply(uncv_draws[,,18:33], 3, sd, na.rm=T)
Omega_est_postcv[3,]  <- apply(uncv_draws[,,18:33], 3, function(x) quantile(x, probs=c(.025), na.rm=T))
Omega_est_postcv[4,]  <- apply(uncv_draws[,,18:33], 3, function(x) quantile(x, probs=c(.975), na.rm=T))

matrix(Omega_est_postcv[1,], nrow=4, byrow=F)
matrix(Omega_est_postcv[3,], nrow=4, byrow=F)
matrix(Omega_est_postcv[4,], nrow=4, byrow=F)



### Now generate summaries of the generated quantities draws

# Fr each quantity, I want mean, CrI LB, CrI UB, SD

#Code for it:

#Now create a df for each quantity I'd like to plot (6 in total), 4 columns, one for each quantity
# theta_est, pi_est, mu_est, DoD_theta, DoD_pi, DoD_mu

timepoints <-  c("Baseline", "3 Month", "6 Month", "12 Month", "3 Month x trt", "6 Month x trt", "12 Month x trt")
modelnames <- c("1RI", "2RI", "INDcv", "IND", "CScv", "CS", "ARcv", "AR", "ADcv", "AD", "UNcv", "UN")
post_stat <- c("mean", "sd", "LB", "UB")

# theta_est_post
theta_est_post <- array(dim = c(4,7,12))
dimnames(theta_est_post)[[1]] <-  post_stat
dimnames(theta_est_post)[[2]] <- timepoints
dimnames(theta_est_post)[[3]] <- modelnames

theta_est_post[1,,] <- apply(theta_est, c(3,4), mean, na.rm=T)
theta_est_post[2,,]  <- apply(theta_est, c(3,4), sd, na.rm=T)
theta_est_post[3,,]  <- apply(theta_est, c(3,4), function(x) quantile(x, probs=c(.025), na.rm=T))
theta_est_post[4,,]  <- apply(theta_est, c(3,4), function(x) quantile(x, probs=c(.975), na.rm=T))


# pi_est_post
pi_est_post <- array(dim = c(4,7,12))
dimnames(pi_est_post)[[1]] <- post_stat
dimnames(pi_est_post)[[2]] <- timepoints
dimnames(pi_est_post)[[3]] <- modelnames

pi_est_post[1,,] <- apply(pi_est, c(3,4), mean, na.rm=T)
pi_est_post[2,,] <- apply(pi_est, c(3,4), sd, na.rm=T)
pi_est_post[3,,] <- apply(pi_est, c(3,4), function(x) quantile(x, probs=c(.025), na.rm=T))
pi_est_post[4,,] <- apply(pi_est, c(3,4), function(x) quantile(x, probs=c(.975), na.rm=T))



# mu_est_post
mu_est_post <- array(dim = c(4,7,12))
dimnames(mu_est_post)[[1]] <-  post_stat
dimnames(mu_est_post)[[2]] <- timepoints
dimnames(mu_est_post)[[3]] <- modelnames

mu_est_post[1,,] <- apply(mu_est, c(3,4), mean)
mu_est_post[2,,] <- apply(mu_est, c(3,4), sd)
mu_est_post[3,,] <- apply(mu_est, c(3,4), function(x) quantile(x, probs=c(.025)))
mu_est_post[4,,] <- apply(mu_est, c(3,4), function(x) quantile(x, probs=c(.975)))


####### DoDs

# DoD_pi_post
DoD_pi_post <- array(dim = c(5,4,12))
dimnames(DoD_pi_post)[[1]] <-  c("mean", "sd", "LB", "UB", "P<0")
dimnames(DoD_pi_post)[[2]] <- c("3 Month", "6 Month", "12 Month", "Average")
dimnames(DoD_pi_post)[[3]] <- modelnames

DoD_pi_post[1,,] <- apply(DoD_pi, c(3,4), mean)
DoD_pi_post[2,,] <- apply(DoD_pi, c(3,4), sd)
DoD_pi_post[3,,] <- apply(DoD_pi, c(3,4), function(x) quantile(x, probs=c(.025)))
DoD_pi_post[4,,] <- apply(DoD_pi, c(3,4), function(x) quantile(x, probs=c(.975)))
DoD_pi_post[5,,] <- apply(DoD_pi, c(3,4), function(x) {ecdf(x)(0)})

# DoD_theta_post
DoD_theta_post <- array(dim = c(5,4,12))
dimnames(DoD_theta_post)[[1]] <-  c("mean", "sd", "LB", "UB", "P<0")
dimnames(DoD_theta_post)[[2]] <- c("3 Month", "6 Month", "12 Month", "Average")
dimnames(DoD_theta_post)[[3]] <- modelnames

DoD_theta_post[1,,] <- apply(DoD_theta, c(3,4), mean)
DoD_theta_post[2,,] <- apply(DoD_theta, c(3,4), sd)
DoD_theta_post[3,,] <- apply(DoD_theta, c(3,4), function(x) quantile(x, probs=c(.025)))
DoD_theta_post[4,,] <- apply(DoD_theta, c(3,4), function(x) quantile(x, probs=c(.975)))
DoD_theta_post[5,,] <- apply(DoD_theta, c(3,4), function(x) {ecdf(x)(0)})



# DoD_mu_post
DoD_mu_post <- array(dim = c(5,4,12))
dimnames(DoD_mu_post)[[1]] <-  c("mean", "sd", "LB", "UB", "P<0")
dimnames(DoD_mu_post)[[2]] <- c("3 Month", "6 Month", "12 Month", "Average")
dimnames(DoD_mu_post)[[3]] <- modelnames

DoD_mu_post[1,,] <- apply(DoD_mu, c(3,4), mean)
DoD_mu_post[2,,] <- apply(DoD_mu, c(3,4), sd)
DoD_mu_post[3,,] <- apply(DoD_mu, c(3,4), function(x) quantile(x, probs=c(.025)))
DoD_mu_post[4,,] <- apply(DoD_mu, c(3,4), function(x) quantile(x, probs=c(.975)))
DoD_mu_post[5,,] <- apply(DoD_mu, c(3,4), function(x) {ecdf(x)(0)})



####### Now, we plot!!!!

#probability of reduction (what proportion of DoD's were less than 0)

#First Plot: Overall difference of differences, referenced against 0

#These are given by DoD_mu_post

dim(DoD_mu_post)

#First need to melt into a 2 dimensional tibble:


DoD_mu_plotdat <- as_tibble(reshape2::melt(DoD_mu_post)) %>% 
  rename(Estimate = Var1, Visit = Var2, Model = Var3, DoD=value) %>% 
  filter(Visit != "Average") %>% 
  separate(Visit, into =c("Month", NA), sep="^\\S*\\K\\s+") %>% 
  mutate(Month = as.numeric(Month)) %>% 
  pivot_wider(names_from = Estimate, values_from = DoD)
  
# sz <- ifelse(DoD_mu_plotdat$Model %in% c("RI", "UN"),1.5,1)
# alp <- ifelse(DoD_mu_plotdat$Model %in% c("RI", "UN"),1,.4)
# 
# pd <- position_dodge(.5) #stagger bars for visualization 
# DoD_mean_plot <- ggplot(DoD_mu_plotdat, aes(x=Month, y=mean, colour=Model, group=Model)) + 
#   geom_errorbar(aes(ymin=LB, ymax=UB), alpha=alp, width=.1, position=pd) +
#   geom_line(aes(),position=pd, alpha=alp) +
#   geom_point(position=pd, size=3, shape=21, fill="white", alpha=alp) + # 21 is filled circle
#   geom_hline(yintercept=0, linetype="dashed", 
#              color = "black") +
#   scale_color_brewer(palette="Set1") + 
#   xlab("Month") +
#   ylab("Difference of Differences") +
#   labs(title =  "Difference of Differences in Heavy Drinking",
#          subtitle = "SBIRT - Control Group") +
#   theme_bw() +
#   theme(legend.justification=c(.975,.025),
#         legend.position=c(.975,.025),
#         panel.background = element_blank(),
#         axis.title.x = element_text(size=20),
#         axis.title.y = element_blank(),
#         axis.text.x = element_text(size = 20),
#         panel.grid.major=element_line(colour="grey90"),
#         title = element_text(size = 30))               # Position legend in bottom right
# 
# 
# ggsave("DoDMeanPlot.pdf", plot=DoD_mean_plot )
# 
# 





####
#Plots:
#Zero mean estimate

###########################################################################
#############################     Pi Plot      ##########################
###########################################################################


pi_plotdat_ctrl <- as_tibble(reshape2::melt(pi_est_post[,1:4,])) %>%
  mutate(Group="Ctrl")

pi_plotdat_trt <- as_tibble(reshape2::melt(pi_est_post[,c(1,5:7),])) %>%
  mutate(Group="Trt")

#Make visits match between the ctrl and treatment, used to be "visit 3 x trt", etc.
pi_plotdat_trt$Var2 <- pi_plotdat_ctrl$Var2

pi_plotdat <- rbind(pi_plotdat_ctrl, pi_plotdat_trt) %>% 
  rename(Estimate = Var1, Visit = Var2, Model = Var3) %>% 
  pivot_wider(names_from = Estimate, values_from = value) 

sz <- ifelse(pi_plotdat$Model %in% c("2RI", "UN"),1.5,1)
alp <- ifelse(pi_plotdat$Model %in% c("2RI", "UN"),1,.3)

# Colorblind-friendly palette:
#cbp2 <- c( "#0072B2", "orange", "#009E73", "#CC79A7",  "#56B4E9", "#FF0000","#F0E442" )


pd <- position_dodge(.5) #stagger bars for visualization 
pi_plot <- ggplot(pi_plotdat, aes(x=Visit, y=mean, colour=Model, group=interaction(Model, Group), linetype=Group)) + 
  geom_errorbar(aes(ymin=LB, ymax=UB), alpha=alp, position=pd, linewidth=sz, width=.5) +
  geom_line(aes(),position=pd, alpha=alp, linewidth=sz) +
  geom_point(position=pd, size=3, shape=21, fill="white", alpha=alp) + # 21 is filled circle
  #scale_color_manual(values=cbp2) + 
  xlab("Visit") +
  labs(title =  "Zero Model") +
  theme_bw() +
  theme(#legend.justification=c(.975,.975),
        legend.position="none",
        panel.background = element_blank(),
        axis.title.x = element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 20),
        panel.grid.major=element_line(colour="grey90"),
        plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5),
        title = element_text(size = 25))               # Position legend in bottom right

pi_plot


ggsave("PiPlot.pdf", plot=pi_plot )

###########################################################################
#############################     Theta Plot      ##########################
###########################################################################

theta_plotdat_ctrl <- as_tibble(reshape2::melt(theta_est_post[,1:4,])) %>%
  mutate(Group="Ctrl")

theta_plotdat_trt <- as_tibble(reshape2::melt(theta_est_post[,c(1,5:7),])) %>%
  mutate(Group="Trt")

#Make visits match between the ctrl and treatment, used to be "visit 3 x trt", etc.
theta_plotdat_trt$Var2 <- theta_plotdat_ctrl$Var2

theta_plotdat <- rbind(theta_plotdat_ctrl, theta_plotdat_trt) %>% 
  rename(Estimate = Var1, Visit = Var2, Model = Var3) %>% 
  pivot_wider(names_from = Estimate, values_from = value) 

sz <- ifelse(theta_plotdat$Model %in% c("2RI", "UN"),1.5,1)
alp <- ifelse(theta_plotdat$Model %in% c("2RI", "UN"),1,.3)

# Colorblind-friendly palette:
cbp2 <- c( "#0072B2", "orange", "#009E73"
           , "#CC79A7",  "#56B4E9", "#FF0000","#F0E442" )
                   
pd <- position_dodge(.5) #stagger bars for visualization 
theta_plot <- ggplot(theta_plotdat, aes(x=Visit, y=mean, colour=Model, group=interaction(Model, Group), linetype=Group)) + 
  geom_errorbar(aes(ymin=LB, ymax=UB), alpha=alp, position=pd, linewidth=sz, width=.5) +
  geom_line(aes(),position=pd, alpha=alp, linewidth=sz) +
  geom_point(position=pd, size=3, shape=21, fill="white", alpha=alp) + # 21 is filled circle
  #scale_color_manual(values=cbp2) + 
  xlab("Visit") +
  labs(title =  "Count Model") +
  theme_bw() +
  theme(#legend.justification=c(.025,.025),
        legend.position="none",
        panel.background = element_blank(),
        axis.title.x = element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 20),
        panel.grid.major=element_line(colour="grey90"),
        plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5),
        title = element_text(size = 25))               # Position legend in bottom right

theta_plot


ggsave("thetaPlot.pdf", plot=theta_plot )



###########################################################################
#############################     Mean Plot      ##########################
###########################################################################

#Add raw data means to this
sbirt <- readRDS("~/Documents/Research/loZIBpaper/R files/sbirt_clean.rds")

sbirt <- sbirt |> mutate(group_plot = ifelse(visit == 1, 0, group))

ctrl_mean <- sbirt |> group_by(group_plot, visit) |>
  summarize(mean = mean(heavy, na.rm=T)) |>
  filter(group_plot==0) |>
  select(visit, mean) |>
  mutate(Visit = factor(visit, labels = c("Baseline", "3 Month", "6 Month", "12 Month"))
  )

trt_mean <- sbirt |> group_by(group_plot, visit) |>
  summarize(mean = mean(heavy, na.rm=T)) |>
  filter(group_plot==1) |>
  select(visit, mean) |>
  mutate(Visit = factor(visit, labels = c("3 Month", "6 Month", "12 Month"))
  )

trt_mean = rbind(ctrl_mean[1,], trt_mean)

mu_plotdat_ctrl <- as_tibble(reshape2::melt(mu_est_post[,1:4,])) %>%
  mutate(Group="Ctrl")

mu_plotdat_trt <- as_tibble(reshape2::melt(mu_est_post[,c(1,5:7),])) %>%
  mutate(Group="Trt")

#Make visits match between the ctrl and treatment, used to be "visit 3 x trt", etc.
mu_plotdat_trt$Var2 <- mu_plotdat_ctrl$Var2

mu_plotdat <- rbind(mu_plotdat_ctrl, mu_plotdat_trt) %>% 
  rename(Estimate = Var1, Visit = Var2, Model = Var3) %>% 
  pivot_wider(names_from = Estimate, values_from = value) 

sz <- ifelse(mu_plotdat$Model %in% c("2RI", "UN"),1.5,1)
alp <- ifelse(mu_plotdat$Model %in% c("2RI", "UN"),1,.3)

# Colorblind-friendly palette:
#cbp2 <- c( "#0072B2", "orange", "#009E73"
#           , "#CC79A7",  "#56B4E9", "#FF0000","#F0E442" )
# I think I wqnt to split this into two separate plots, one for each Control and Trt group

#color_vec = scales::hue_pal()(12)
color_vec <- c("#F8766D", "#DE8C00", "#B79F00", "#7CAE00", "#00BA38", "#00C08B", "#00BFC4", "#00B4F0", "#C77CFF", "#F564E3", "#FF64B0", "#619CFF")

pd <- position_dodge(.5) #stagger bars for visualization 
mu_plot <- ggplot() + 
  geom_errorbar(data = mu_plotdat, aes(ymin=LB, ymax=UB, x=Visit, y=mean, colour=Model, group=interaction(Model, Group), linetype=Group), alpha=alp, position=pd, linewidth=sz, width=.5) +
  geom_line(data = mu_plotdat, aes(x=Visit, y=mean, colour=Model, group=interaction(Model, Group), linetype=Group),position=pd, alpha=alp, linewidth=sz) +
  geom_point(data = mu_plotdat, position=pd, size=3, shape=21, fill="white", alpha=alp, aes(x=Visit, y=mean, colour=Model, group=interaction(Model, Group))) + # 21 is filled circle
  geom_line(data = ctrl_mean, mapping = aes(x = Visit, y = mean, group=1), color = "black", linewidth = 1, alpha = .6) +
  geom_point(data = ctrl_mean, mapping = aes(x= Visit, y = mean), color = "black", fill = "black", alpha = .6) + 
  geom_line(data = trt_mean, mapping = aes(x = Visit, y = mean, group=2, linetype="Trt"), color = "black", linewidth = 1, alpha = .6) +
  geom_point(data = trt_mean, mapping = aes(x= Visit, y = mean), color = "black", fill = "black", alpha = .6) + 
  scale_color_manual(values=color_vec) + 
  xlab("Visit") +
  labs(title =  "Mean Days of Heavy Drinking") +
  theme_bw() +
  theme(#legend.justification=c(.975,.975),
        legend.position="right",
        panel.background = element_blank(),
        axis.title.x = element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 20),
        panel.grid.major=element_line(colour="grey90"),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15),
        plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5),
        title = element_text(size = 25))               # Position legend in bottom right

mu_plot




ggsave("muPlot.pdf", plot=mu_plot )


#A version of this but with only UN, 2RI and raw data:

mu_plotdat_UN = mu_plotdat[mu_plotdat$Model %in% c("2RI", "UN"),]

color_vec = c("blue", "darkorange")

pd = position_dodge(.2)

mu_plot <- ggplot() + 
  geom_errorbar(data = mu_plotdat_UN, aes(ymin=LB, ymax=UB, x=Visit, y=mean, colour=Model, group=interaction(Group, Model), linetype=Group), position=pd, width=.5, linewidth = 1) +
  geom_line(data = mu_plotdat_UN, aes(x=Visit, y=mean, colour=Model, group=interaction(Group, Model), linetype=Group),position=pd, linewidth = 1) +
  geom_point(data = mu_plotdat_UN, position=pd, size=3, shape=21, aes(x=Visit, y=mean, colour=Model, fill = Model, group=interaction(Group, Model))) + # 21 is filled circle
  geom_line(data = ctrl_mean, mapping = aes(x = Visit, y = mean, group=1), color = "black", linewidth = 1) +
  geom_point(data = ctrl_mean, mapping = aes(x= Visit, y = mean), color = "black", fill = "black", size = 4, shape = 17) + 
  geom_line(data = trt_mean, mapping = aes(x = Visit, y = mean, group=2, linetype="Trt"), color = "black", linewidth = 1) +
  geom_point(data = trt_mean, mapping = aes(x= Visit, y = mean), color = "black", fill = "black", size = 4, shape = 17) + 
  scale_color_manual(values=color_vec) + 
  scale_fill_manual(values = color_vec) +
  xlab("Visit") +
  labs(title =  "Mean Days of Heavy Drinking") +
  theme_bw() +
  theme(#legend.justification=c(.975,.975),
    legend.position="right",
    panel.background = element_blank(),
    axis.title.x = element_text(size=20),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 20),
    panel.grid.major=element_line(colour="grey90"),
    legend.title = element_text(size=15),
    legend.text = element_text(size=15),
    plot.title = element_text(hjust = .5),
    plot.subtitle = element_text(hjust = .5),
    title = element_text(size = 25))               # Position legend in bottom right

mu_plot

ggsave("muPlot.pdf", plot=mu_plot )
# I like this one better

#Quesetion: Can I do half violin plots of the posterior distribution as error bars? I think this would make a good additional plot - posterior distributions of group means at each time point with the observed group means included for reference.


###########################################################################
#############################     DoD Pi Plot      ##########################
###########################################################################
#DoD Pi Plot



DoD_pi_plotdat <- as_tibble(reshape2::melt(DoD_pi_post)) %>% 
  rename(Estimate = Var1, Visit = Var2, Model = Var3, DoD=value) %>% 
  filter(Visit != "Average") %>% 
  #separate(Visit, into =c("Month", NA), sep="^\\S*\\K\\s+") %>% 
  #mutate(Month = as.numeric(Month)) %>% 
  pivot_wider(names_from = Estimate, values_from = DoD)

sz <- ifelse(DoD_pi_plotdat$Model %in% c("2RI", "UN"),1.5,1)
alp <- ifelse(DoD_pi_plotdat$Model %in% c("2RI", "UN"),1,.3)

# Colorblind-friendly palette:
cbp2 <- c( "#0072B2", "orange", "#009E73"
                    , "#CC79A7",  "#56B4E9", "#FF0000","#F0E442" )

pd <- position_dodge(.5) #stagger bars for visualization 
DoD_pi_plot <- ggplot(DoD_pi_plotdat, aes(x=Visit, y=mean, colour=Model, group=Model)) + 
  geom_errorbar(aes(ymin=LB, ymax=UB), alpha=alp, position=pd, linewidth=1.1, width=.5) +
  geom_line(aes(),position=pd, alpha=alp, linewidth=sz) +
  geom_point(position=pd, size=3, shape=21, fill="white", alpha=alp) + # 21 is filled circle
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black") +
  #scale_color_manual(values=cbp2) + 
  #scale_x_continuous(breaks=c(0, 3, 6, 12)) +
  xlab("Month") +
  ylab("Difference of Differences") +
  labs(title =  "Zero Model DoD") +
  theme_bw() +
  theme(legend.justification=element_blank(),
        legend.position= "none",
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        panel.grid.major=element_line(colour="grey90"),
        title = element_text(size = 25),
        plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5))               # Position legend in bottom right

DoD_pi_plot

setwd("/Users/bwrogers/Documents/Research/loZIBpaper/Latex Files")

ggsave("DoDPiPlot.pdf", plot=DoD_pi_plot )





###########################################################################
#############################     DoD Theta Plot      ##########################
###########################################################################
# DoD Theta Plot


DoD_theta_plotdat <- as_tibble(reshape2::melt(DoD_theta_post)) %>% 
  rename(Estimate = Var1, Visit = Var2, Model = Var3, DoD=value) %>% 
  filter(Visit != "Average") %>% 
  #separate(Visit, into =c("Month", NA), sep="^\\S*\\K\\s+") %>% 
  #mutate(Month = as.numeric(Month)) %>% 
  pivot_wider(names_from = Estimate, values_from = DoD)

sz <- ifelse(DoD_pi_plotdat$Model %in% c("2RI", "UN"),1.5,1)
alp <- ifelse(DoD_pi_plotdat$Model %in% c("2RI", "UN"),1,.3)

# Colorblind-friendly palette:
cbp2 <- c( "#0072B2", "orange", "#009E73"
           , "#CC79A7",  "#56B4E9", "#FF0000","#F0E442" )
                    
pd <- position_dodge(.5) #stagger bars for visualization 
DoD_theta_plot <- ggplot(DoD_theta_plotdat, aes(x=Visit, y=mean, colour=Model, group=Model)) + 
  geom_errorbar(aes(ymin=LB, ymax=UB), alpha=alp, position=pd, linewidth=sz, width=.5) +
  geom_line(aes(),position=pd, alpha=alp, linewidth=1.1) +
  geom_point(position=pd, size=3, shape=21, fill="white", alpha=alp) + # 21 is filled circle
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black") +
  #scale_color_manual(values=cbp2) + 
  #scale_x_continuous(breaks=c(0, 3, 6, 12)) +
  xlab("Month") +
  ylab("Difference of Differences") +
  labs(title =  "Count Model DoD") +
  theme_bw() +
  theme(legend.justification=element_blank(),
        legend.position= "none",
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        panel.grid.major=element_line(colour="grey90"),
        title = element_text(size = 25),
        plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5))               # Position legend in bottom right

DoD_theta_plot

setwd("/Users/bwrogers/Documents/Research/loZIBpaper/Latex Files")

ggsave("DoDThetaPlot.pdf", plot=DoD_theta_plot )



###########################################################################
#############################     DoD Mu Plot      ##########################
###########################################################################
#DoD Mu Plot



DoD_mu_plotdat <- as_tibble(reshape2::melt(DoD_mu_post)) %>% 
  rename(Estimate = Var1, Visit = Var2, Model = Var3, DoD=value) %>% 
  filter(Visit != "Average") %>% 
  #separate(Visit, into =c("Month", NA), sep="^\\S*\\K\\s+") %>% 
  #mutate(Month = as.numeric(Month)) %>% 
  pivot_wider(names_from = Estimate, values_from = DoD)

sz <- ifelse(DoD_pi_plotdat$Model %in% c("2RI", "UN"),1.5,1)
alp <- ifelse(DoD_pi_plotdat$Model %in% c("2RI", "UN"),1,.3)

# Colorblind-friendly palette:
cbp2 <- c( "#0072B2", "orange", "#009E73"
           , "#CC79A7",  "#56B4E9", "#FF0000","#F0E442" )
                    
pd <- position_dodge(.5) #stagger bars for visualization 
DoD_mu_plot <- ggplot(DoD_mu_plotdat, aes(x=Visit, y=mean, colour=Model, group=Model)) + 
  geom_errorbar(aes(ymin=LB, ymax=UB), alpha=alp, position=pd, linewidth=sz, width=.5) +
  geom_line(aes(),position=pd, alpha=alp, linewidth=1.1) +
  geom_point(position=pd, size=3, shape=21, fill="white", alpha=alp) + # 21 is filled circle
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black") +
  #scale_color_manual(values=cbp2) + 
  #scale_x_continuous(breaks=c(0, 3, 6, 12)) +
  xlab("Visit") +
  ylab("Difference of Differences") +
  labs(title =  "DoD in Days of Heavy Drinking",
       subtitle = "SBIRT - Control Group") +
  theme_bw() +
  theme(#legend.justification=c(.975,.975),
    legend.position="right",
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        panel.grid.major=element_line(colour="grey90"),
        title = element_text(size = 25),
        plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5),
        legend.text=element_text(size=15),
        legend.title=element_blank())               # Position legend in bottom right

DoD_mu_plot

setwd("/Users/bwrogers/Documents/Research/loZIBpaper/Latex Files")

ggsave("DoDMeanPlot.pdf", plot=DoD_mean_plot )


#Now I'd like to arrange DoD plots on a single grid

library(ggpubr)
theme_set(theme_pubr())

# #With labels
# DoD_all_plot <- ggarrange(
#   DoD_mu_plot,                # First row plot
#   # Second row with box and dot plots
#   ggarrange(DoD_theta_plot, DoD_pi_plot, ncol = 2, labels = c("B.", "C."), font.label = list(size = 30, font=4)), 
#   nrow = 2, 
#   labels = "A.",       # Label of the top plot
#   font.label = list(size = 30)
# ) 

DoD_all_plot <- ggarrange(
  DoD_mu_plot,                # First row plot
  # Second row with box and dot plots
  ggarrange(DoD_theta_plot, DoD_pi_plot, ncol = 2), 
  nrow = 2
) 

DoD_all_plot

ggsave("DoDAllPlots.pdf", plot=DoD_all_plot, width=10, height=7 )



mean_all_plot <- ggarrange(
  mu_plot,                # First row plot
  # Second row with box and dot plots
  ggarrange(theta_plot, pi_plot, ncol = 2), 
  nrow = 2
) 

mean_all_plot

ggsave("MeanAllPlots.pdf", plot=mean_all_plot,  width=10, height=7  )





############################################################################################
####################### Correlation Posterior Density Plots ################################
############################################################################################

#### Plot of Omega parameters posterior distributions
omega_plotdat <- as_tibble(reshape2::melt(un_draws[,,21:36])) %>% 
  mutate(chain=as.factor(chain)) %>% 
  filter(variable %in% c("Omega2[1,2]", "Omega2[1,3]","Omega2[1,4]","Omega2[2,3]", "Omega2[2,4]", "Omega2[3,4]"))
head(omega_plotdat)


  ggplot(omega_plotdat, aes(x=value, fill=chain, color=chain)) +
  geom_density(alpha=.5) +
    theme(        axis.title.x = element_blank(),
                  legend.position = "none") +
    geom_vline(xintercept = 0, linetype="dashed", 
                color = "black", size=.5) +
  facet_wrap(~variable)
  
  #To fix this plot: 
  #     1) Change title of each facet, 
  #     2) lay the facets out to mirror a correlation matrix
  #     3) Add title
  #     
  
  
  #Code to create a grid of plots:

  library(gridExtra)
  
  row.labs <- c("3 Month", "6 Month", "12 Month", "", "", "")
  col.labs <- c("","","Baseline", "", "3 Month", "6 Month")
  
  lowerCrI <- quantile()
  
  plots <- lapply(1:6, function(id) 
    filter(omega_plotdat, variable == unique(omega_plotdat$variable)[id]) %>% 
      ggplot(aes(x=value, fill=chain, color=chain)) +
      geom_density(alpha=.5) +
      scale_x_continuous(limits=c(-1,1),labels = seq(-1, 1, by = .5)) +
      theme(legend.position = "none",
            axis.title.y=element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.title.x = element_blank(),
            plot.title = element_text(hjust = 0.5)
            ) +
      geom_vline(xintercept = 0, linetype="dashed", 
                 color = "black", size=.5) +
      labs(title=unique(omega_plotdat$variable)[id], hjust="center")
      
    )
  
  grid.titles = c("Baseline", "3 Month", "6 Month", "12 Month")

  
  
  m <- matrix(NA, 3, 3)
  m[upper.tri(m, diag = T)] <- 1:6
  grid.arrange(grobs = plots, layout_matrix = m)
  g <- arrangeGrob(grobs = plots, layout_matrix = m)
  
  
  ggsave(file="CorrPost.pdf", g)
  
  
  
  #Next Coding task:
  
  #Calculate Difference of Differences from data
  
  diffs_data <- sbirt |> select(id, group, visit, heavy) |> 
    filter(heavy < 80) |>
    pivot_wider(values_from = heavy, names_from = visit) |> 
    rename(
     visit1 = 3,
      visit2 = 4,
      visit3 = 5,
      visit4 = 6
    ) |>
    mutate(DoD1 = visit2 - visit1, DoD2 = visit3 - visit1 , DoD3 = visit4 - visit1) |>
    group_by(group) |>
    summarize(mean_DoD1 = mean(DoD1, na.rm=T), mean_DoD2 = mean(DoD2, na.rm=T), mean_DoD3 = mean(DoD3, na.rm=T))
    
  DoDs_data = diffs_data[2,2:4] - diffs_data[1, 2:4]
  DoDs_data
  

  outlier_subjects <- sbirt |> select(id, group, visit, heavy) |>
    filter(heavy >80) |> pull(id) |> unique()

  outlier_subjects  

  outliers <-   sbirt |> select(id, group, visit, heavy) |> 
    filter(id %in% outlier_subjects) |>
    pivot_wider(values_from = heavy, names_from = visit) 

  apply(outliers[outliers$group==1,3:6], 2, mean, na.rm = TRUE)  
  apply(outliers[outliers$group==0,3:6], 2, mean, na.rm = TRUE) 
  
  
  
  
  
  #Posterior samples plot
  