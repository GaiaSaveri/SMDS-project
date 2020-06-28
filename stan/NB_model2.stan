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
  int<lower=1> D; //number of departments
  int<lower=1, upper=D> id[N];
}
parameters {
  real<lower=0> inv_phi;   
  real beta;               
  real beta_d; 
  vector[J] alpha;            
  real<lower=0> sigma; 
  real mu;
}
transformed parameters {
  real phi = inv(inv_phi);
}
model {
  alpha ~ normal(mu + id*beta_d, sigma);
  sigma ~ normal(0, 1);
  mu ~ normal(log(4), 1); 
  beta ~ normal(0.5, 1);
  beta_d ~ normal(0,1);
  inv_phi ~ cauchy(0., 5);
  
  cases ~ neg_binomial_2_log(alpha[id], phi);
} 
generated quantities {
  int y_rep[N];
  vector[N] log_lik;
  for (n in 1:N) {
    real eta_n = alpha[id[n]] + beta * time[n];
    y_rep[n] = neg_binomial_2_log_safe_rng(eta_n, phi);
    log_lik[n] = neg_binomial_2_log_lpmf(cases[n]|eta_n, phi);  
  }
}
