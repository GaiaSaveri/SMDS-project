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
      //real dI_dt =  beta * I * S - gamma * I;
      real dI_dt =  beta * I * S / N - gamma * I;
      real dR_dt =  gamma * I;
      
      return {dS_dt, dI_dt, dR_dt};
  }
}
data {
  int<lower=1> n_days; //number of days observed
  real y0[3];
  real t0; //initial time point
  real ts[n_days-1]; //time points observed
  int N; //population size
  int cases[n_days]; //total number of infected individual 
}
transformed data {
  real x_r[0];
  int x_i[1] = { N };
  //int x_i[0];
}
parameters {
  real<lower=0> gamma;
  real<lower=0> beta;
  //real<lower=0> theta[2]; //beta=theta[1], gamma=theta[2]
  real<lower=0> phi_inv;
  //real<lower=0, upper=1> S0; //initial fraction of susceptible individuals
}
transformed parameters{
  //ODE solutions
  real<lower=0, upper=1> y[n_days, 3];
  //initial SIR fractions
  //real<lower=0, upper=1> y_init[3];
  //real<lower=0> lambda[n_days];
  
  //y_init[1] = S0;
  //y_init[2] = 1 - S0;
  //y_init[3] = 0;
  
  //runge-kutta solver
  //y_hat = integrate_ode_rk45(sir, y_init, t0, ts, theta, x_r, x_i);
  
  //for(i in 1:n_days) {
    //lambda[i] = y_hat[i,2] * N;
  //}
  real incidence[n_days]; //this is new
  
  real phi = 1. / phi_inv;
  {
    real theta[2];
    theta[1] = beta;
    theta[2] = gamma;

    y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
  }
  
  for(i in 1:n_days-1) {
    incidence[i] = y[i,1] - y[i+1,1];
  }
  incidence[n_days] = cases[n_days];
}
model {
  //priors
  //beta ~ lognormal(log(0.4), 0.5); 
  beta ~ normal(log(0.4), 0.5); 
  //gamma ~ lognormal(log(0.1), 0.2);
  gamma ~ normal(log(0.1), 0.2);
  phi_inv ~ exponential(5);
  //theta ~ lognormal(0,1);
  //S0 ~ beta(1,1);
  
  //sampling distribution
  cases ~ neg_binomial_2(incidence, phi);
  //cases ~ poisson(lambda);
}
generated quantities {
  //real R0 = theta[1]/theta[2];

  real R0 = beta / gamma;
  real recovery_time = 1 / gamma;
  real pred_cases[n_days];
  //vector[n_days] log_lik;
  pred_cases = neg_binomial_2_rng(incidence, phi);

}
