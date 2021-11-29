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
  int<lower=0> N;

  vector[N] Cobs_caecum;
  vector[N] Cobs_cephalon;
  vector[N] Cobs_intestin;
  vector[N] Cobs_reste;

  real C0_caecum;
  real C0_cephalon;
  real C0_reste;
  real C0_intestin;

  int N_vCw;
  vector[N_vCw] vCw;

  int<lower=0> posTd;
  int<lower=0> Post2;
  int<lower=0> Post5;
  int<lower=0> Post7;
  int<lower=0> Post9;
  int<lower=0> Post15;
  int<lower=0> Post21;
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

  vector[posTd] Cpred_caecum;
  vector[posTd] Cpred_cephalon;
  vector[posTd] Cpred_reste;
  vector[posTd] Cpred_intestin;

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

  for(i in 2:posTd){
    Cpred_intestin[i] = Cpred_intestin[i-1] + (
      ku[1]*vCw[i] - ke[1]*Cpred_intestin[i-1] +
      k[4,1]*Cpred_reste[i-1] + k[2,1]*Cpred_caecum[i-1] + k[3,1]*Cpred_cephalon[i-1] -
      (k[1,2]+k[1,3]+k[1,4])*Cpred_intestin[i-1] )*p ;
    Cpred_caecum[i] = Cpred_caecum[i-1] + (
      ku[2]*vCw[i] - ke[2]*Cpred_caecum[i-1] +
      k[1,2]*Cpred_intestin[i-1] + k[3,2]*Cpred_cephalon[i-1] + k[4,2]*Cpred_reste[i-1] -
      (k[2,1]+k[2,3]+k[2,4])*Cpred_caecum[i-1] )*p ;
    Cpred_cephalon[i] = Cpred_cephalon[i-1] + (
      ku[3]*vCw[i] - ke[3]*Cpred_cephalon[i-1] +
      k[1,3]*Cpred_intestin[i-1] + k[2,3]*Cpred_caecum[i-1]+k[4,3]*Cpred_reste[i-1] -
      (k[3,1]+k[3,2]+k[3,4])*Cpred_cephalon[i-1] )*p ;
    Cpred_reste[i] = Cpred_reste[i-1] + (
      ku[4]*vCw[i] - ke[4]*Cpred_reste[i-1] +
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

  for(i in 1:3){
    Cobs_intestin[i] ~ normal(Cpred_intestin[Post2],sigma[1]);
    Cobs_caecum[i] ~ normal(Cpred_caecum[Post2],sigma[2]);
    Cobs_cephalon[i] ~ normal(Cpred_cephalon[Post2],sigma[3]);
    Cobs_reste[i] ~ normal(Cpred_reste[Post2],sigma[4]);
  }
  for(i in 4:6){
    Cobs_intestin[i] ~ normal(Cpred_intestin[Post5],sigma[1]);
    Cobs_caecum[i] ~ normal(Cpred_caecum[Post5],sigma[2]);
    Cobs_cephalon[i] ~ normal(Cpred_cephalon[Post5],sigma[3]);
    Cobs_reste[i] ~ normal(Cpred_reste[Post5],sigma[4]);
  }
  for(i in 7:9){
    Cobs_intestin[i] ~ normal(Cpred_intestin[Post7],sigma[1]);
    Cobs_caecum[i] ~ normal(Cpred_caecum[Post7],sigma[2]);
    Cobs_cephalon[i] ~ normal(Cpred_cephalon[Post7],sigma[3]);
    Cobs_reste[i] ~ normal(Cpred_reste[Post7],sigma[4]);
  }
  for(i in 10:12){
    Cobs_intestin[i] ~ normal(Cpred_intestin[Post9],sigma[1]);
    Cobs_caecum[i] ~ normal(Cpred_caecum[Post9],sigma[2]);
    Cobs_cephalon[i] ~ normal(Cpred_cephalon[Post9],sigma[3]);
    Cobs_reste[i] ~ normal(Cpred_reste[Post9],sigma[4]);
  }
  for(i in 13:15){
    Cobs_intestin[i] ~ normal(Cpred_intestin[Post15],sigma[1]);
    Cobs_caecum[i] ~ normal(Cpred_caecum[Post15],sigma[2]);
    Cobs_cephalon[i] ~ normal(Cpred_cephalon[Post15],sigma[3]);
    Cobs_reste[i] ~ normal(Cpred_reste[Post15],sigma[4]);
  }
  for(i in 16:18){
    Cobs_intestin[i] ~ normal(Cpred_intestin[Post21],sigma[1]);
    Cobs_caecum[i] ~ normal(Cpred_caecum[Post21],sigma[2]);
    Cobs_cephalon[i] ~ normal(Cpred_cephalon[Post21],sigma[3]);
    Cobs_reste[i] ~ normal(Cpred_reste[Post21],sigma[4]);
  }
}

generated quantities {

  vector[posTd] Cgen_caecum;
  vector[posTd] Cgen_cephalon;
  vector[posTd] Cgen_reste;
  vector[posTd] Cgen_intestin;

  for(t in 1:posTd){
    Cgen_intestin[t] = normal_rng(Cpred_intestin[t], sigma[1]) ;
    Cgen_caecum[t] = normal_rng(Cpred_caecum[t], sigma[2]) ;
    Cgen_cephalon[t] = normal_rng(Cpred_cephalon[t], sigma[3]) ;
    Cgen_reste[t] = normal_rng(Cpred_reste[t], sigma[4]) ;
  }
}


