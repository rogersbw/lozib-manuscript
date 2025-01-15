data{
int<lower=0> P; //P is the number of predictor var
int<lower=0> J; //Number of subjects
int<lower=0> N; //Number of obs
int<lower=1> V; //number of visits
array[N] int<lower=0, upper=J> ll; //vector of each subject id associated with each y
array[N] int<lower=1, upper=V> visit; //visit number associated with each observation

array[N] row_vector[P] X; // NxP design matrix
array[N] int<lower=0> y; // vector of outcomes

}

parameters{
vector[P] beta1;
vector[P] beta2;

array[J] real  gamma1;
array[J] vector[V] gamma2z;

real<lower=0> sigma1;
vector<lower=0>[V] sigma2; 
real<lower=0> psi;

array[J] real ri_cs; //RI using CS = RI + sigma2^2 I parameterization 
real<lower=0> sigma2_cs; //RI std dev

}

transformed parameters {
  vector[N] theta;
  vector[N] pie;
  array[J] vector[V] gamma2;
  
  vector[N] log_lik;
  
  for(j in 1:J){
    gamma2[j] = sigma2_cs*ri_cs[j] + sigma2 .* gamma2z[j];
  }

  for(n in 1:N){
    theta[n] = inv_logit(X[n]*beta1 + sigma1*gamma1[ll[n]]);
    pie[n] = inv_logit(X[n]*beta2 + psi*gamma1[ll[n]] + gamma2[ll[n],visit[n]]);
    
  if (y[n] == 0) {
      log_lik[n] = log_sum_exp(bernoulli_lpmf(0 | theta[n]),
                            bernoulli_lpmf(1 | theta[n])
                              + binomial_lpmf(y[n]| 90, pie[n]));
    } else {
      log_lik[n] = bernoulli_lpmf(1 | theta[n])
                  + binomial_lpmf(y[n] | 90, pie[n]);
    }
    
  }
}

model{
  //Priors
  
  //Betas are iid N(0,5)
  beta1 ~ multi_normal(rep_vector(0,P),diag_matrix(rep_vector(5,P)));
  beta2 ~ multi_normal(rep_vector(0,P),diag_matrix(rep_vector(5,P)));
  
  //Since Sigmas are restricted to positive, these are half normal
  sigma1 ~ normal(0, .5);

  for(v in 1:V){
    sigma2[v] ~ normal(0, .5);
  }

  sigma2_cs ~ normal(0, .5);
  
  psi ~ normal(0,1);
  
  for(j in 1:J){
    gamma1[j] ~ normal(0, 1);
    gamma2z[j] ~ normal(0, 1);
    ri_cs[j] ~ normal(0,1);
  }
  
  //likelihood
  target += sum(log_lik);
}


