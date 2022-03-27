functions{
  real janoschek(real A, real K, real B, real v, real t){
    return K-(K-A)*2^(-(t^v)/B);
  // A: the lower asymptote
  // K: the upper asymptote
  // B: doubling time
  // v: move maximum growth near one of the asymptote
  // t: time in minutes
  }
}

data {
  int<lower=1> N; //number of total data points
  int<lower=1> M; //number of replicates
  int replicate[N]; //replicate number
  real time[N]; //time vector
  real od[N]; //od vector
}

parameters {
  real<lower=0> A;
  real<lower=0> K;
  real<lower=0> B;
  real<lower=0> v;

  vector[M] A_rep;
  vector<lower=0>[M] K_rep;
  vector<lower=0>[M] B_rep;
  vector<lower=0>[M] v_rep;

  real<lower=0> sigma_od;
  real<lower=0> sigma_A;
  real<lower=0> sigma_K;
  real<lower=0> sigma_B;
  real<lower=0> sigma_v;
}

model {

  A ~ cauchy(0,0.1);
  K ~ cauchy(0,2);
  B ~ cauchy(30,100);
  v ~ cauchy(1,1);

  sigma_A ~ cauchy(0,0.5);
  sigma_K ~ cauchy(0,5);
  sigma_B ~ cauchy(0,100);
  sigma_v ~ cauchy(0,1);

  A_rep ~ normal(A,sigma_A);
  K_rep ~ normal(K,sigma_K);
  B_rep ~ normal(B,sigma_B);
  v_rep ~ normal(v,sigma_v);

  for (n in 1:N) {
        od[n] ~ normal(janoschek(A_rep[replicate[n]],K_rep[replicate[n]],B_rep[replicate[n]],v_rep[replicate[n]],time[n]), sigma_od);
  }

}

generated quantities {

  vector[101] time_pred;
  matrix[101,M] OD_pred;
  vector[101] OD_pred_all;

  for (i in 1:101){
    time_pred[i] <- max(time)/100 * (i-1);
    for (k in 1:M){
      OD_pred[i,k] <- janoschek(A_rep[k],K_rep[k],B_rep[k],v_rep[k],time_pred[i]);
    }
    OD_pred_all[i] <- janoschek(A,K,B,v,time_pred[i]);
  }
}
