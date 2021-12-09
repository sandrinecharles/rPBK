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

real interpolate(real x, vector xpt, vector ypt){
  if(x >= min(xpt) && x <= max(xpt)){
    int idx = findfirst(x, xpt) ;
    return ypt[idx] + (x - xpt[idx]) * (ypt[idx+1] - ypt[idx]) / (xpt[idx+1] - xpt[idx]) ;
  } else{
    return 0.0 ;
  }
}

real[] ode_pbtk( real t,      // time
                 real[] y,    // variables
                 real[] theta,
                 real[] x_r,
                 int[] x_i) {

  // parameters
  int N_Cw = x_i[1] ;

  real ku[4] = theta[1:4] ;
  real ke[4] = theta[5:8] ;
  real k1[4] = theta[9:12] ;
  real k2[4] = theta[13:16] ;
  real k3[4] = theta[17:20] ;
  real k4[4] = theta[21:24] ;

  // vector[1+n_met] dydt ;
  real dydt[4] ;

  real tacc = x_r[1] ;
  vector[N_Cw] tp_Cw = to_vector(x_r[2:(N_Cw+1)]) ;
  vector[N_Cw] Cw = to_vector(x_r[(N_Cw+2):(N_Cw+1+N_Cw)]) ;

  k1[1] = 0;
  k2[2] = 0;
  k3[3] = 0;
  k4[4] = 0;

  if(t <= tacc){
    // Accumulation
    dydt[1] = ku[1] * interpolate(t, tp_Cw, Cw) - ke[1] * y[1] +
       k1[1]*y[2] + k2[1]*y[2] + k3[1]*y[3] + k4[1]*y[4] - (k1[1]+k1[2]+k1[3]+k1[4]) * y[1] ;

    dydt[2] = ku[2] * interpolate(t, tp_Cw, Cw) - ke[2] * y[2] +
       k1[2]*y[2] + k2[2]*y[2] + k3[2]*y[3] + k4[2]*y[4] - (k2[1]+k2[2]+k2[3]+k2[4]) * y[2] ;

    dydt[3] = ku[1] * interpolate(t, tp_Cw, Cw) - ke[3] * y[3] +
       k1[3]*y[2] + k2[1]*y[2] + k3[1]*y[3] + k4[1]*y[4] - (k3[1]+k3[2]+k3[3]+k3[4]) * y[3] ;

    dydt[4] = ku[1] * interpolate(t, tp_Cw, Cw) - ke[4] * y[4] +
       k1[4]*y[2] + k2[1]*y[2] + k3[1]*y[3] + k4[1]*y[4] - (k4[1]+k4[2]+k4[3]+k4[4]) * y[4] ;
  } else{
    // Depuration
    dydt[1] = - ke[1] * y[1] +
       k1[1]*y[2] + k2[1]*y[2] + k3[1]*y[3] + k4[1]*y[4] - (k1[1]+k1[2]+k1[3]+k1[4]) * y[1] ;

    dydt[2] = - ke[2] * y[2] +
       k1[2]*y[2] + k2[2]*y[2] + k3[2]*y[3] + k4[2]*y[4] - (k2[1]+k2[2]+k2[3]+k2[4]) * y[2] ;

    dydt[3] = - ke[3] * y[3] +
       k1[3]*y[2] + k2[1]*y[2] + k3[1]*y[3] + k4[1]*y[4] - (k3[1]+k3[2]+k3[3]+k3[4]) * y[3] ;

    dydt[4] = - ke[4] * y[4] +
       k1[4]*y[2] + k2[1]*y[2] + k3[1]*y[3] + k4[1]*y[4] - (k4[1]+k4[2]+k4[3]+k4[4]) * y[4] ;

  }
  return(dydt) ;
}

