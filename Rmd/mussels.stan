data {
  int<lower=0> n_obs;      // number of observations 
  int y_i[n_obs];          // counts
  real<lower=0> km[n_obs]; // covariates 
  real<lower=0> my_km;     // which km to predict 
}
parameters {
  real beta_0;
  real beta_1; 
}
transformed parameters {
  vector[n_obs] lambda;
  for(i in 1:n_obs){
    lambda[i] = exp(beta_0 + beta_1*km[i]); 
  }
}
model {
  // specify the priors 
  beta_0 ~ normal(0,3); 
  beta_1 ~ normal(0,3);
  
  // specify the likelihood
  y_i ~ poisson(lambda); 
  
  // NOTE: Stan is doing ln(post) = ln(priors) + ln(like) in background
}
generated quantities {
 int<lower=0> y_pred = poisson_rng(exp(beta_0 + beta_1*my_km));
}
