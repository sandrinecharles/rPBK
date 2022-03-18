library(dplyr)
library(tidyr)
library(ggplot2)
library(expm)
library(Matrix)

##################
param = list(
  ku = c(1917,1571,91.1,135),
  ke = c(.506,.06,.053,.026),
  k = matrix(0,nrow=4,ncol=4),
  tp_Cw = c(0,2,5,7,9,15,21),
  Cw = rep(0.01108,7),
  tacc = 7
)

E <- param$k
diag(E) <- 0
E_diag <- - param$ke - E %*% rep(1,4)
diag(E) <- E_diag
U <- param$ku
Cx = param$Cw[1]
time = param$tp_Cw
tacc = param$tacc

sol = sapply(1:length(time), function(i){
  if(time[i] < tacc){
    Matrix::solve(E) %*% (Matrix::expm(time[i] * E)-diag(4)) %*% U * Cx
  } else{
    Matrix::solve(E) %*% (Matrix::expm(time[i] * E) - Matrix::expm((time[i] - tacc) * E) ) %*% U * Cx
  }
})
sol_ = t(do.call("cbind", sol))
df_out_x = as.data.frame(as.matrix(sol_))

df_out_x = df_out_x %>%
  dplyr::mutate(time = time) %>%
  tidyr::pivot_longer(cols = -"time")

ggplot() +
  theme_minimal() +
  geom_point(data = df_out_x,
             aes(x=time,y=value,color=name))

