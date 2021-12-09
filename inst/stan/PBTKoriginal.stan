//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;// number of data per replicate
  int<lower=0> N_rep; // Number of replicate

  matrix[N,N_rep] Cobs_caecum;
  matrix[N,N_rep] Cobs_cephalon;
  matrix[N,N_rep] Cobs_intestin;
  matrix[N,N_rep] Cobs_reste;

  int<lower=0> tp_position[N];

  real C0_caecum;
  real C0_cephalon;
  real C0_reste;
  real C0_intestin;

  int<lower=0> N_Cw;
  vector[N_Cw] tp_Cw;

  real p;

}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[4] log10ku;
  vector[4] log10ke;

  vector[4] log10k1;
  vector[4] log10k2;
  vector[4] log10k3;
  vector[4] log10k4;

  real<lower=0> sigma[4];

}

transformed parameters{

  vector[4] ku;
  vector[4] ke;

  real k[4,4];

  vector[N_Cw] Cpred_caecum;
  vector[N_Cw] Cpred_cephalon;
  vector[N_Cw] Cpred_reste;
  vector[N_Cw] Cpred_intestin;

  for(i in 1:4){
    ku[i] = 10^log10ku[i];
    ke[i] = 10^log10ke[i];
    k[1,i] = 10^log10k1[i];
    k[2,i] = 10^log10k2[i];
    k[3,i] = 10^log10k3[i];
    k[4,i] = 10^log10k4[i];
  }

  Cpred_caecum[1] = C0_caecum;
  Cpred_cephalon[1] = C0_cephalon;
  Cpred_reste[1] = C0_reste;
  Cpred_intestin[1] = C0_intestin;

  for(i in 2:N_Cw){
    Cpred_intestin[i] = Cpred_intestin[i-1] + (
      ku[1]*tp_Cw[i] - ke[1]*Cpred_intestin[i-1] +
      k[4,1]*Cpred_reste[i-1] + k[2,1]*Cpred_caecum[i-1] + k[3,1]*Cpred_cephalon[i-1] -
      (k[1,2]+k[1,3]+k[1,4])*Cpred_intestin[i-1] )*p ;
    Cpred_caecum[i] = Cpred_caecum[i-1] + (
      ku[2]*tp_Cw[i] - ke[2]*Cpred_caecum[i-1] +
      k[1,2]*Cpred_intestin[i-1] + k[3,2]*Cpred_cephalon[i-1] + k[4,2]*Cpred_reste[i-1] -
      (k[2,1]+k[2,3]+k[2,4])*Cpred_caecum[i-1] )*p ;
    Cpred_cephalon[i] = Cpred_cephalon[i-1] + (
      ku[3]*tp_Cw[i] - ke[3]*Cpred_cephalon[i-1] +
      k[1,3]*Cpred_intestin[i-1] + k[2,3]*Cpred_caecum[i-1]+k[4,3]*Cpred_reste[i-1] -
      (k[3,1]+k[3,2]+k[3,4])*Cpred_cephalon[i-1] )*p ;
    Cpred_reste[i] = Cpred_reste[i-1] + (
      ku[4]*tp_Cw[i] - ke[4]*Cpred_reste[i-1] +
      k[1,4]*Cpred_intestin[i-1] + k[2,4]*Cpred_caecum[i-1] + k[3,4]*Cpred_cephalon[i-1] -
      (k[4,1]+k[4,2]+k[4,3])*Cpred_reste[i-1] )*p ;
  }
}

model {

  target += uniform_lpdf(log10ku | -5, 5);
  target += uniform_lpdf(log10ke | -5, 5);

  target += uniform_lpdf(log10k1 | -5, 5);
  target += uniform_lpdf(log10k2 | -5, 5);
  target += uniform_lpdf(log10k3 | -5, 5);
  target += uniform_lpdf(log10k4 | -5, 5);

  target += gamma_lpdf(sigma | 0.001,0.001);

  for(rep in 1:3){
    for(t in 2:7){
      target += normal_lpdf(Cobs_intestin[t,rep] | Cpred_intestin[tp_position[t]],sigma[1]);
      target += normal_lpdf(Cobs_caecum[t,rep] | Cpred_caecum[tp_position[t]],sigma[2]);
      target += normal_lpdf(Cobs_cephalon[t,rep] | Cpred_cephalon[tp_position[t]],sigma[3]);
      target += normal_lpdf(Cobs_reste[t,rep] | Cpred_reste[tp_position[t]],sigma[4]);
    }
  }
}
generated quantities {
  vector[N] Cgen_caecum;
  vector[N] Cgen_cephalon;
  vector[N] Cgen_reste;
  vector[N] Cgen_intestin;

  for(t in 1:N){
    Cgen_intestin[t] = normal_rng(Cpred_intestin[tp_position[t]], sigma[1]) ;
    Cgen_caecum[t] = normal_rng(Cpred_caecum[tp_position[t]], sigma[2]) ;
    Cgen_cephalon[t] = normal_rng(Cpred_cephalon[tp_position[t]], sigma[3]) ;
    Cgen_reste[t] = normal_rng(Cpred_reste[tp_position[t]], sigma[4]) ;
  }
}


