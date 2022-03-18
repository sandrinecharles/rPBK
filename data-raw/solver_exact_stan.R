library(rstan)
library(dplyr)
library(tidyr)
library(ggplot2)

##################

model_ode <- "functions{

  matrix exact_AD_long(vector time, int n_time, real tacc, matrix E, matrix I, vector U, int k, real Cx){
    matrix[k,k] E_inv = inverse(E) ;
    matrix[k,n_time] E_out ;
    for(i in 1:n_time){
      if(time[i] < tacc){
        E_out[1:k,i] = E_inv * (matrix_exp(time[i] * E) - I) * U * Cx ;
      } else {
        E_out[1:k,i] = E_inv * (matrix_exp(time[i] * E) - matrix_exp((time[i] - tacc) * E)) * U * Cx ;
      }
    }
    return( E_out' ) ;
  }

  matrix matrix_E(vector ke, matrix k, int N_k){

    matrix[N_k, N_k] m = add_diag(k, 0) ;
    vector[N_k] diag_k = - ke - m * rep_vector(1,N_k) ;
    matrix[N_k, N_k] matrix_E = add_diag(k, diag_k) ;

    return(matrix_E) ;
  }

  matrix matrix_I(int k){
    matrix[k,k] I = add_diag(rep_matrix(0, k, k), 1) ;
    return(I);
  }

}"

expose_stan_functions(stanc(model_code = model_ode))

################################################################################

ku = c(1917,1571,91.1,135)
ke = c(.506,.06,.053,.026)
k = matrix(0,nrow=4,ncol=4)
time = c(0,2,5,7,9,15,21)
exposure = rep(0.01108,7)
time_acc = 7

E = matrix_E(ke, k, 4)
I = matrix_I(4)

sol = exact_AD_long(
  time,
  length(time),
  time_acc,
  E,I,
  ku,
  length(ke),
  exposure)

sol_ = t(sol)
df_out_x = as.data.frame(as.matrix(sol_))

df_out_x = df_out_x %>%
  dplyr::mutate(time = par$tp_Cw) %>%
  tidyr::pivot_longer(cols = -"time")

ggplot(data = df_out_x,
       aes(x=time,y=value,color=name)) +
  theme_minimal() +
  geom_line() +
  geom_point()



sol = exact_AD_long(
  par$tp_Cw,
  length(par$tp_Cw),
  par$tacc,
  E,I,
  par$ku,
  length(par$ke),
  par$Cw[1])
