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
  int<lower=0> N_time;// number of data per replicate
  int<lower=0> N_rep; // Number of replicate
  int<lower=0> N_comp;

  real Cobs_comp[N_time,N_rep,N_comp];

  int<lower=0> tp_position[N_time];

  real C0_comp[N_comp];

  int<lower=0> N_Cw;
  vector[N_Cw] tp_Cw;
  vector[N_Cw] Cw;

  real p;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[N_comp] log10ku;
  vector[N_comp] log10ke;

  matrix[N_comp,N_comp] log10k;

  real<lower=0> sigma[N_comp];

}

transformed parameters{

  vector[N_comp] ku;
  vector[N_comp] ke;

  real k[N_comp,N_comp];

  real Cpred_comp[N_Cw,N_comp];

  for(i in 1:N_comp){
    ku[i] = 10^log10ku[i];
    ke[i] = 10^log10ke[i];
    for(j in 1:N_comp){
       k[i,j] = 10^log10k[i,j];
    }
  }


  for(j in 1:N_comp){
    Cpred_comp[1,j] = C0_comp[j] ;
    for(i in 2:N_Cw){
      Cpred_comp[i,j] = Cpred_comp[i-1,j] + (
        ku[j]*Cw[i] - ke[1]*Cpred_comp[i-1,j] +
        k[j,1]*Cpred_comp[i-1,1] + k[j,2]*Cpred_comp[i-1,2] + k[j,3]*Cpred_comp[i-1,3] + k[j,4]*Cpred_comp[i-1,4] -
        (k[1,j]+k[2,j]+k[3,j]+k[4,j])*Cpred_comp[i-1,j] )*p ;
    }
  }
}

model {

  target += uniform_lpdf(log10ku | -5, 5);
  target += uniform_lpdf(log10ke | -5, 5);

  for(i_comp in 1:N_comp){
     target += uniform_lpdf(log10k[1:N_comp,i_comp] | -5, 5);
  }

  target += gamma_lpdf(sigma | 0.01, 0.01);

  for(rep in 1:N_rep){
    for(t in 2:N_time){
      for(i_comp in 1:N_comp){
        target += normal_lpdf(Cobs_comp[t,rep,i_comp] | Cpred_comp[tp_position[t],i_comp],sigma[i_comp]);
      }
    }
  }
}
generated quantities {
  real Cgen_comp[N_time,N_comp];

  for(t in 1:N_time){
    for(i_comp in 1:N_comp){
      Cgen_comp[t,i_comp] = normal_rng(Cpred_comp[tp_position[t],i_comp], sigma[i_comp]) ;
    }
  }
}


