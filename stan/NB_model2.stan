functions {
  /*
  * Alternative to neg_binomial_2_log_rng() that 
  * avoids potential numerical problems during warmup
  */
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
  vector<lower=0>[N] log_people;
  //int<lower=1> D; //number of departments
  //int<lower=1, upper=D> dep[N];
  //int<lower=1, upper=D> id[N];
}
parameters {
  real<lower=0> inv_phi;   
  real beta_t; //time coefficient               
  //real beta_d; 
  real beta_a; //age coefficient
  real alpha;            
  real<lower=0> sigma; //sd of alpha
  real mu; //mean of alpha
  
}
transformed parameters {
  real phi = inv(inv_phi);
}
model {
  vector[N] eta = alpha + beta_t * time + beta_a * age + log_people;
  
  alpha ~ normal(mu , sigma);
  sigma ~ normal(0, 1);
  mu ~ normal(log(4), 1); 
  beta_t ~ normal(0.5, 1);
  beta_a ~ normal(0,1);
  inv_phi ~ cauchy(0., 5);
  
  cases ~ neg_binomial_2_log(eta, phi);
  //alpha ~ normal(mu, sigma);
  //sigma ~ normal(0,1);
  //mu ~ normal(log(4), 1);
  
  cases ~ neg_binomial_2_log(alpha+beta_t*time+beta_a*age, phi);
} 
generated quantities {
  int y_rep[N];
  vector[N] log_lik;
  for (n in 1:N) {
    real eta_n = alpha + beta_t * time[n] + beta_a * age[n] + log_people[n];
    y_rep[n] = neg_binomial_2_log_safe_rng(eta_n, phi);
    log_lik[n] = neg_binomial_2_log_lpmf(cases[n]|eta_n, phi);  
  }
}
