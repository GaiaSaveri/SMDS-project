data {
  int<lower=1> N;
  int<lower=0> cases[N];
  vector<lower=0>[N] time;
}
parameters {
  real alpha;
  real beta;
}
model {
  beta ~ normal(0.25, 1);
  alpha ~ normal(log(4), 1);
  
  // poisson_log(eta) is more efficient and stable alternative to poisson(exp(eta))
  cases ~ poisson_log(alpha + beta * time);
} 
generated quantities {
  // sample predicted values from the model for posterior predictive checks
  int<lower=0> y_rep[N];
  vector[N] log_lik;
  for (n in 1:N) {
    real eta_n = alpha + beta * time[n];
    y_rep[n] = poisson_log_rng(eta_n);
    log_lik[n] = poisson_log_lpmf(cases[n]|eta_n);
  }
}
