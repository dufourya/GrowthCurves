functions{
  real[] log_growth_ode(real t, real[] y, real[] theta, real[] x_r, int[] x_i){
    real dydt[1];
    dydt[1] <- theta[1] * y[1] * (1 - y[1] / theta[2]);
    return dydt;
  }
}

data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> idx[K];
  real time[N];
  real OD[N];
  real t0;
}

transformed data {
  real x_r[0];
  int x_i[0];
}

parameters {
  real<lower=0, upper=10> A;
  real<lower=1, upper=3000> mu;
  real<lower=0, upper=max(OD)> L;
  real<lower=0, upper=max(OD)> B;

  vector<lower=0, upper=10>[K] A_rep;
  vector<lower=1, upper=3000>[K] mu_rep;
  vector<lower=0, upper=max(OD)>[K] L_rep;
  vector<lower=0, upper=max(OD)>[K] B_rep;

  real<lower=0> sigma_OD;

  real<lower=0> sigma_A;
  real<lower=0> sigma_mu;
  real<lower=0> sigma_L;
  real<lower=0> sigma_B;
}

model {
  
  real OD_hat[N];
  int pos;
  real theta[2];
  real y0[1];

  pos <- 1;
  
  A_rep ~ normal(A,sigma_A);
  mu_rep ~ normal(mu,sigma_mu);
  L_rep ~ normal(L,sigma_L);
  B_rep ~ normal(B, sigma_B);

  for (k in 1:K) {
    real test[idx[k],2];
    theta[2] <- A_rep[k];
    theta[1] <- mu_rep[k];
    y0[1] <- L_rep[k];
    //print("test=", test, " y=", integrate_ode(log_growth_ode, y0, t0, segment(time, pos, idx[k]), theta, x_r, x_i));
    test <- integrate_ode(log_growth_ode, y0, t0, segment(time, pos, idx[k]), theta, x_r, x_i);
    print(test);
    //segment(OD_hat, pos, idx[k]) <- test;
    pos <- pos + idx[k];
  }
  
  OD ~ normal(OD_hat, sigma_OD);
  
}

generated quantities {

  vector[100] time_pred;
  matrix[100,K] OD_pred;
  vector[100] OD_pred_all;
  
}