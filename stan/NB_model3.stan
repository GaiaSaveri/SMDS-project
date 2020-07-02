functions {
  int neg_binomial_2_log_safe_rng(real eta, real phi) {
    real gamma_rate = gamma_rng(phi, phi / exp(eta));
    if (gamma_rate >= exp(20.79))
      return -9;
      
    return poisson_rng(gamma_rate);
  }
}
data {
  int<lower=1> N;                     
  int<lower=0> cases[N];              
  vector<lower=0>[N] time;  
  vector<lower=1>[N] age;              
  
  //department-level data
  int<lower=1> K; //number of attributes for the dep
  int<lower=1> J; //number of departments
  int<lower=1, upper=J> dep_id[N];
  matrix[J,K] dep_data;
}
parameters {
  real<lower=0> inv_phi;   
  real beta_t; //time coefficient  
  real beta_a; //age coefficient              
  vector[J] alpha_d; //intercept for each department            
  real<lower=0> sigma_alpha; //slope for each department
  real mu; //intercept for the mean of alpha_d          
  vector[K] zeta; //coefficients for the dep-specific attributes  
  //to fix divergences
  //vector[J] alpha; //<-- new! 
}
transformed parameters {
  real phi = inv(inv_phi);
  //reparametrization
  //vector[J] alpha_d = mu + dep_data * zeta + sigma_alpha * alpha;
}
model {
  //alpha ~ normal(0,1);
  alpha_d ~ normal(mu + dep_data * zeta, sigma_alpha);
  sigma_alpha ~ normal(0, 1);
  mu ~ normal(log(4), 1);
  zeta ~ normal(0, 1); 
  
  
  beta_t ~ normal(0.5, 1);
  beta_a ~ normal(0,1);
  //inv_phi ~ cauchy(0., 5);
  inv_phi ~ normal(0, 1);
  cases ~ neg_binomial_2_log(alpha_d[dep_id] + beta_t*time+beta_a*age, phi);
} 
generated quantities {
  int y_rep[N];
  vector[N] log_lik;
  for (n in 1:N) {
    real eta_n = alpha_d[dep_id[n]] + beta_t*time[n]+beta_a*age[n];
    y_rep[n] = neg_binomial_2_log_safe_rng(eta_n, phi);
    log_lik[n] = neg_binomial_2_log_lpmf(cases[n]|eta_n, phi); 
  }
}

