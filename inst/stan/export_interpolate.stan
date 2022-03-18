functions {
#include /include/interpolation.stan
}
data {
  real x;
  int N;
  vector[N] xpt;
  vector[N] ypt;
}
parameters {
}
model {
}
generated quantities {
  real y;
  y = interpolate(x, xpt, ypt);
}
