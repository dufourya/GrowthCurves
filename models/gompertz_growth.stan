functions{
  real gompertz(real A, real mu, real L, real t, real B){
    return B + (A*exp(log(L)*exp(-t/mu)));
  }
}

data {
  int N;
  int K;
  int idx[K];
  vector[N] time;
  vector[N] OD;
}

parameters {
  real<lower=0, upper=1> A;
  real<lower=1, upper=1000> mu;
  real<lower=0> L;
  real<lower=0, upper=0.1> B;

  vector<lower=0, upper=1>[K] A_rep;
  vector<lower=1, upper=1000>[K] mu_rep;
  vector<lower=0>[K] L_rep;
  vector<lower=0, upper=0.1>[K] B_rep;

  real<lower=0> sigma_OD;
  real<lower=0> sigma_A;
  real<lower=0> sigma_mu;
  real<lower=0> sigma_L;
  real<lower=0> sigma_B;
}

model {

int pos;

  A ~ normal(0.5,1);
  mu ~ normal(50,10);
  L ~ exponential(0.01);
  B ~ exponential(0.001);


pos <- 1;
  
  A_rep ~ normal(A,sigma_A);
  mu_rep ~ normal(mu,sigma_mu);
  L_rep ~ normal(L,sigma_L);
  B_rep ~ normal(B, sigma_B);

for (k in 1:K) {
    for (i in pos:(pos+idx[k]-1)) {
        OD[i] ~ normal(gompertz(A_rep[k], mu_rep[k], L_rep[k], time[i], B_rep[k]), sigma_OD);
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
      OD_pred[i,k] <- gompertz(A_rep[k], mu_rep[k], L_rep[k], time_pred[i], B_rep[k]);
    }
    OD_pred_all[i] <- gompertz(A, mu, L, time_pred[i], B);
  }
}