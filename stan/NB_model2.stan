functions {
<<<<<<< HEAD
=======
  /*
  * Alternative to neg_binomial_2_log_rng() that 
  * avoids potential numerical problems during warmup
  */
>>>>>>> d1161e2cf231d60c7b1be97bc1eca5dbee8877f1
  int neg_binomial_2_log_safe_rng(real eta, real phi) {
    real gamma_rate = gamma_rng(phi, phi / exp(eta));
    if (gamma_rate >= exp(20.79))
      return -9;
      
    return poisson_rng(gamma_rate);
  }
}
data {
<<<<<<< HEAD
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
=======
  int<lower=0> D; //number of departments
  
  int<lower=1> N;
  int<lower=0> cases[N];
  int<lower=1, upper=D> dep[N];
  vector<lower=0>[N] time;
}
parameters {
  vector[D] alpha;
  real<lower=0> sigma;
  real beta;
  real mu;
  real<lower=0> inv_phi;
>>>>>>> d1161e2cf231d60c7b1be97bc1eca5dbee8877f1
}
transformed parameters {
  real phi = inv(inv_phi);
}
model {
<<<<<<< HEAD
  alpha ~ normal(mu + id*beta_d, sigma);
  sigma ~ normal(0, 1);
  mu ~ normal(log(4), 1); 
  beta ~ normal(0.5, 1);
  beta_d ~ normal(0,1);
  inv_phi ~ cauchy(0., 5);
  
  cases ~ neg_binomial_2_log(alpha[id], phi);
=======
  alpha ~ normal(mu, sigma);
  sigma ~ normal(0,1);
  mu ~ normal(log(4), 1);
  //vector[N] eta = alpha + beta * time;
  //alpha ~ normal(0, 10);
  beta ~ normal(0, 2.5);
  //inv_phi ~ normal(0, 1);
  inv_phi ~ cauchy(0., 5);
  
  cases ~ neg_binomial_2_log(alpha[dep]+beta*time, phi);
>>>>>>> d1161e2cf231d60c7b1be97bc1eca5dbee8877f1
} 
generated quantities {
  int y_rep[N];
  vector[N] log_lik;
  for (n in 1:N) {
<<<<<<< HEAD
    real eta_n = alpha[id[n]] + beta * time[n];
=======
    real eta_n = alpha[dep[n]] + beta * time[n];
>>>>>>> d1161e2cf231d60c7b1be97bc1eca5dbee8877f1
    y_rep[n] = neg_binomial_2_log_safe_rng(eta_n, phi);
    log_lik[n] = neg_binomial_2_log_lpmf(cases[n]|eta_n, phi);  
  }
}
