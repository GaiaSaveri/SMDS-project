functions {
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {

      real S = y[1];
      real I = y[2];
      real R = y[3];
      real N = x_i[1];
      
      real beta = theta[1];
      real gamma = theta[2];
      
      real dS_dt = -beta * I * S / N;
      real dI_dt =  beta * I * S / N - gamma * I;
      real dR_dt =  gamma * I;
      
      return {dS_dt, dI_dt, dR_dt};
  }
}
data {
  int<lower=1> n_days;
  real y0[3];
  real t0;
  real ts[n_days];
  int N;
  int cases[n_days];
}
transformed data {
  real x_r[0];
  int x_i[1] = { N };
}
parameters {
  real<lower=0, upper=1> gamma;
  real<lower=0, upper=1> beta;
  real<lower=0> phi_inv;
}
transformed parameters{
  real y[n_days, 3];
  real phi = 1. / phi_inv;
  {
    real theta[2];
    theta[1] = beta;
    theta[2] = gamma;

    y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
  }
}
model {
  //priors
  //beta ~ normal(2, 1);--> no good 
  beta ~ lognormal(log(0.4), 0.5); 
  //beta ~ uniform(0.60, 0.80);
  //gamma ~ uniform(0.35, 0.50);
  //gamma ~ normal(0.4, 0.5); --> no_good
  gamma ~ normal(0.4, 0.2);
  //phi_inv ~ exponential(5);
  phi_inv ~ cauchy(0., 5);
  
  //sampling distribution
  //col(matrix x, int n) - The n-th column of matrix x. Here the number of infected people 
  cases ~ neg_binomial_2(col(to_matrix(y), 2), phi);
}
generated quantities {
  real R0 = beta / gamma;
  real recovery_time = 1 / gamma;
  real pred_cases[n_days];
  vector[n_days] log_lik;
  //pred_cases = neg_binomial_2_rng(col(to_matrix(y), 2), phi);
  for(n in 1:n_days) {
    pred_cases[n] = neg_binomial_2_rng(y[n,2], phi);
    log_lik[n] = neg_binomial_2_log_lpmf(cases[n]| col(to_matrix(y), 2), phi);
  }
}

