data {
  int<lower=1> N;
  int<lower=0> cases[N];
  vector<lower=0>[N] time;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> inv_phi;
}
transformed parameters {
  real phi = inv(inv_phi);
}
model {
  vector[N] eta = alpha + beta * time;
  alpha ~ normal(0, 10);
  beta ~ normal(0, 2.5);
  //inv_phi ~ normal(0, 1);
  inv_phi ~ cauchy(0., 5);
  
  cases ~ neg_binomial_2_log(eta, phi);
} 
generated quantities {
  int y_rep[N];
  vector[N] log_lik;
  for (n in 1:N) {
    real eta_n = alpha + beta * time[n];
    y_rep[n] = neg_binomial_2_log_rng(eta_n, phi);
    log_lik[n] = neg_binomial_2_log_lpmf(cases[n]|eta_n, phi);  
  }
}
