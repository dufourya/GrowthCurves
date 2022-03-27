functions{

  real generalized_logistic(real A, real K, real C, real Q, real B, real v, real t){
    return A+(K-A)/((C+2^((Q-t)/exp(B)))^exp(v));
  // A: the lower asymptope
  // K: the upper asymptope
  // B: doubling time
  // C: typically 1
  // Q: related to the number of growing cells at t = 0
  // v: move maximum growth near one of the symptopes (v>0, v=1 -> classic logistic)
  // t: time in minutes
  }

}

data {

  int<lower=1> N; //number of total data points
  int<lower=1> M; //number of replicates
  int replicate[N]; //replicate number
  vector[N] time; //time vector
  vector[N] od; //od vector
  real eK;
  real eA;
  real eB;
  real eQ;
  //real ev;

}

transformed data {

  vector[M] C_rep;
  real C;
  C = 1;


  for (i in 1:M){
    C_rep[i] = 1;
  }

}

parameters {

  real<lower=0> A;
  real<lower=0> K;
  real B;
  real Q;
  real v;

  vector<lower=0>[M] A_rep;
  vector<lower=0>[M] K_rep;
  vector[M] B_rep;
  vector[M] Q_rep;
  vector[M] v_rep;

  real<lower=0> sigma_od;
  real<lower=0> sigma_A;
  real<lower=0> sigma_K;
  real<lower=0> sigma_B;
  real<lower=0> sigma_Q;
  real<lower=0> sigma_v;

}

model {

  A ~ normal(eA,0.25);
  K ~ lognormal(log(eK),0.25);
  B ~ normal(log(eB),10);
  Q ~ normal(eQ,10);
  v ~ normal(0,1);

  sigma_A ~ exponential(1);
  sigma_K ~ exponential(1);
  sigma_B ~ exponential(1);
  sigma_Q ~ exponential(1);
  sigma_v ~ exponential(1);
  sigma_od ~ exponential(10);

  A_rep ~ normal(A,sigma_A);
  K_rep ~ normal(K,sigma_K);
  B_rep ~ normal(B,sigma_B);
  Q_rep ~ normal(Q,sigma_Q);
  v_rep ~ normal(v,sigma_v);

  for (n in 1:N) {
        od[n] ~ normal(generalized_logistic(A_rep[replicate[n]],K_rep[replicate[n]],C_rep[replicate[n]],Q_rep[replicate[n]],B_rep[replicate[n]],v_rep[replicate[n]],time[n]),sigma_od);
  }

}

generated quantities {

  vector[101] time_pred;
  matrix[101,M] OD_pred;
  vector[101] OD_pred_all;
  vector[M] min_relative_doubling_rep;
  real min_relative_doubling;

for (k in 1:M){
  min_relative_doubling_rep[k] = exp(B_rep[k])*(exp(v_rep[k])+1)*(K_rep[k]-A_rep[k]+A_rep[k]*((exp(v_rep[k])+1)/exp(v_rep[k]))^(exp(v_rep[k])))/(exp(v_rep[k])*(K_rep[k]-A_rep[k]));
}

min_relative_doubling = exp(B)*(exp(v)+1)*(K-A+A*((exp(v)+1)/exp(v))^(exp(v)))/(exp(v)*(K-A));

  for (i in 1:101){
    time_pred[i] = max(time)/100 * (i-1);
    for (k in 1:M){
      OD_pred[i,k] = generalized_logistic(A_rep[k],K_rep[k],C_rep[k],Q_rep[k],B_rep[k],v_rep[k],time_pred[i]);
    }
    OD_pred_all[i] = generalized_logistic(A,K,C,Q,B,v,time_pred[i]);
  }

}
