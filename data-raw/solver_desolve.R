library(rstan)
library(deSolve)
library(dplyr)
library(tidyr)
library(ggplot2)

##################
# Build the interpolation function from Stan
model_code <-
  '
  functions {
  int findfirst(real t, vector xt) {
  int i = 0 ;
  if(t == max(xt)){
    i = num_elements(xt) - 1 ;
    return i ;
  } else if(t < min(xt) || t > max(xt)){
    return i ;
  } else {
    while (t >= xt[i+1]){
      i = i+1 ;
    }
    return i ;
  }
}

    real export_interpolate(real x, vector xpt, vector ypt){
  if(x >= min(xpt) && x <= max(xpt)){
    int idx = findfirst(x, xpt) ;
    return ypt[idx] + (x - xpt[idx]) * (ypt[idx+1] - ypt[idx]) / (xpt[idx+1] - xpt[idx]) ;
  } else{
    return 0.0 ;
  }
}
  }
'
expose_stan_functions(stanc(model_code = model_code))


##################


param = list(
  ku = c(1917,1571,91.1,135),
  ke = c(.506,.06,.053,.026),
  k = matrix(0,nrow=4,ncol=4),
  tp_Cw = c(0,2,5,7,9,15,21),
  Cw = rep(0.01108,7),
  tacc = 7
)

model <- function(t, y, param) {
  with(as.list(param), {
    for(i in 1:4){
      if(t < tacc){
        dydt[i] = ku[i]*export_interpolate(t,tp_Cw,Cw) - ke[i]*y[i] +
          as.numeric(t(k[i,1:4]) %*% y[1:4]) - sum(k[1:4,i]) * y[i]
      } else{
        dydt[i] = - ke[i]*y[i] +
          as.numeric(t(k[i,1:4]) %*% y[1:4]) - sum(k[1:4,i]) * y[i]
      }
    }
    return(list(dydt))
  })
}

y = c(.806,.647,0.067,0.048)
times = seq(0, 21, length.out = 100)
dydt = rep(0,4)
out <- ode(y,times, model, parms = param)

df_out = out %>%
  tidyr::as_tibble() %>%
  dplyr::rename("V1"=`1`, "V2"=`2`, "V3"=`3`, "V4"=`4`) %>%
  tidyr::pivot_longer(cols = 2:5)

ggplot() +
  theme_minimal() +
  geom_line(data = df_out,
            aes(x=time,y=value,color=name)) +
  geom_point(data = df_out_A,
             aes(x=time,y=value,color=name))
