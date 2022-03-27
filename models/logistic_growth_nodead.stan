functions{
  real log_growth(real A, real mu, real L, real t){
    return A*L*exp(t/mu)/(A+L*(exp(t/mu)-1));
  }
}

data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> idx[K];
  vector<lower=0>[N] time;
  vector[N] OD;
  real<lower=0> maxr;
}

parameters {
  real<lower=0, upper=5> A;
  real<lower=10, upper=maxr> mu;
  real<lower=0, upper=max(OD)> L;

  vector<lower=0, upper=5>[K] A_rep;
  vector<lower=10, upper=maxr>[K] mu_rep;
  vector<lower=0, upper=max(OD)>[K] L_rep;

  real<lower=0> sigma_OD;

  real<lower=0> sigma_A;
  real<lower=0> sigma_mu;
  real<lower=0> sigma_L;
}

model {
  int pos;
  pos <- 1;

  //A ~ normal(1,1);
  //mu ~ normal(200,10);
  //L ~ normal(0.05,0.1);
  //B ~ normal(0.1,0.1);

  //sigma_A ~ exponential(0.1);
  //sigma_mu ~ exponential(0.1);
  //sigma_L ~ exponential(0.1);
  //sigma_B ~ exponential(0.1);
  
  A_rep ~ normal(A,sigma_A);
  mu_rep ~ normal(mu,sigma_mu);
  L_rep ~ normal(L,sigma_L);

  for (k in 1:K) {
    for (i in pos:(pos+idx[k]-1)) {
        OD[i] ~ normal(log_growth(A_rep[k], mu_rep[k], L_rep[k], time[i]), sigma_OD);
    }
        pos <- pos + idx[k];
  }
}

generated quantities {

  vector[101] time_pred;
  matrix[101,K] OD_pred;
  vector[101] OD_pred_all;

for (i in 1:101){
    time_pred[i] <- max(time)/100 * (i-1);
    for (k in 1:K){
      OD_pred[i,k] <- log_growth(A_rep[k], mu_rep[k], L_rep[k], time_pred[i]);
    }
    OD_pred_all[i] <- log_growth(A, mu, L, time_pred[i]);
}
}