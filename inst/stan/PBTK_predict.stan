functions {

#include /include/interpolation.stan

}
data {
  int N_predict ;
  real tp_predict[N_predict] ;

  real tacc;

  int<lower=0> N_predict_Cw;
  vector[N_predict_Cw] tp_predict_Cw;
  vector[N_predict_Cw] Cw_predict;

  real C0_caecum;
  real C0_cephalon;
  real C0_reste;
  real C0_intestin;

  int N_samples;
  real theta[N_samples, 24];
  matrix[N_samples, 4] sigma;
}
transformed data{

  real t0;
  real y0[4];

  int<lower=0> x_int[1];

  real x_r[1+N_predict_Cw+N_predict_Cw];

  t0 = 0;
  y0[1] = C0_intestin ;
  y0[2] = C0_caecum ;
  y0[3] = C0_cephalon ;
  y0[4] = C0_reste ;

  x_int[1] = N_predict_Cw ;

  x_r[1] = tacc ;
  for(i in 1:N_predict_Cw){
    x_r[1+i] = tp_predict_Cw[i] ;
  }
  for(i in 1:N_predict_Cw){
    x_r[i+N_predict_Cw+1] = Cw_predict[i] ;
  }
}
parameters {
}
model {
}
generated quantities {

  real<lower=0> y_sim[N_predict, 4] ;

  matrix[N_samples,N_predict] Cgen_caecum;
  matrix[N_samples,N_predict] Cgen_cephalon;
  matrix[N_samples,N_predict] Cgen_reste;
  matrix[N_samples,N_predict] Cgen_intestin;

  for(i in 1:N_samples) {
    y_sim = integrate_ode_rk45(ode_pbtk, y0, t0, tp_predict, theta[i,1:24], x_r, x_int) ;

    // // compartment:
    // Cpred_intestin[t] = y_sim[t,1] ;
    // Cpred_caecum[t] = y_sim[t,2] ;
    // Cpred_cephalon[t] = y_sim[t,3] ;
    // Cpred_reste[t] = y_sim[t,4] ;

    for(t in 1:N_predict){
      Cgen_intestin[i,t] = normal_rng(y_sim[t,1], sigma[i,1]) ;
      Cgen_caecum[i,t] = normal_rng(y_sim[t,2], sigma[i,2]) ;
      Cgen_cephalon[i,t] = normal_rng(y_sim[t,3], sigma[i,3]) ;
      Cgen_reste[i,t] = normal_rng(y_sim[t,4], sigma[i,4]) ;
    }
  }
}


