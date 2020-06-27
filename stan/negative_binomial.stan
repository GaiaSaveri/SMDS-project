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
  //vector[N] mu;
  //vector[N] log_lik;
  vector[N] y_rep;
  //mu = exp(eta);
  for (n in 1:N) {
    real eta_n = alpha + beta * time[n];
    y_rep[n] = neg_binomial_2_rng(exp(eta_n), phi);
    //log_lik[i] = neg_binomial_2_log_lpmf(cases[i] | eta[i], phi);
    //y_rep[i] = neg_binomial_2_rng(mu[i], phi);
  }
}
