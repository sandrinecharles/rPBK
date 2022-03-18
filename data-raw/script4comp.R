# General entries ####
# Empty working directory
rm(list = ls())
# Load required library, if installed
library(tidyverse)
library(crosstalk)
library(deSolve)

# Load parameter values ####
rates <- read.table("data-raw/param4comp.txt", header = TRUE)
# Create vector U, named ku
ku <- rates %>%
  filter(str_detect(rates$param, "ku"))
# Extract elimination rates as a vector
ke <- rates %>%
  filter(str_detect(rates$param, "ke"))
# Extract sigma estimates as a vector
sigma <- rates %>%
  filter(str_detect(rates$param, "sigma"))
##### ATTENTION: the same sigma was used for all compartments #####
# Extract kij rates
kij <- rates %>%
  filter(str_detect(rates$param, "ku", negate = TRUE),
         str_detect(rates$param, "ke", negate = TRUE),
         str_detect(rates$param, "sigma", negate = TRUE))
kij[,c("to","from")] <- str_split_fixed(kij$param, pattern = "", n = 3)[,2:3]
# Extract the nbr of compartments
nc <- nrow(ku)

# Experimental conditions ####
# Exposure concentration
cx <- 11.1
# Accumulation duration
tc <- 7 # in days
# Total experiment duration
tf <- tc + 14

# Simulation from the generic matrix solution ####
# Create matrix E
matrixE <- matrix(NA, ncol = nc, nrow = nc)
for(i in 1:nc){
  for(j in 1:nc){
    if(i==j){
      matrixE[i,j] <- - ke[i, "median"] - sum(kij[kij$from==i, "median"])
    }
    else{
      matrixE[i,j] <- as.numeric(kij %>%
                                   filter(to==i, from==j) %>%
                                   select(median))
    }
  }
}
# Eigen analysis of matrix E
X <- eigen(matrixE)
# Build matrix with eigenvectors as columns
P <- X$vectors
# Create matrix J with eigenvalues on the diagonal
J <- diag(x = X$values)
# Check the Jordan form of matrix E: P %*% J %*% solve(P)

## Accumulation phase ####
# Create a time vector for the accumulation phase
taccu <- seq(0, tc, length.out = 250)
# Simulation of the internal concentration in [0 ; t_acc]
# Initiate the matrix to store simulations
caccu <- matrix(NA, nrow = length(taccu), ncol = nc)
# Initial condition at t = 0
caccu[1,] <- 0
# Loop through the time vector `taccu`
intE_ls = list()
for(i in 2:length(taccu)){
  intJ <- diag(x = (exp(X$values * taccu[i]) - 1) / X$values)
  intE <- P %*% intJ %*% solve(P)
  intE_ls[[i]] = intE
  caccu[i,] <- intE %*% as.vector(ku[,"median"]) * cx / 1000
}

## Depuration phase ####
# Create a time vector for the depuration phase
tdepu <- seq(tc, tf, length.out = 250)
# Simulation of the internal concentration in [t_acc ; t_fin]
# Initiate the matrix to store simulations
cdepu <- matrix(NA, nrow = length(tdepu), ncol = nc)
# Initial condition at t = t_acc
cdepu[1,] <- caccu[nrow(caccu),]
# Loop through the time vector `tdepu`
for(i in 2:length(tdepu)){
  intJ <- diag(x = (1 - exp(- X$values * tc)) *
                 exp(X$values * tdepu[i]) / X$values)
  intE <- P %*% intJ %*% solve(P)
  cdepu[i,] <- intE %*% as.vector(ku[,"median"]) * cx / 1000
}

## Plots of internal concentrations within compartments ####
# Plot window parameters
par(mar = c(4, 5, 0.1, 0.1), mfrow = c(2, 2))
# Generate plots
for(i in 1:nc){
  ymax <- max(caccu[,i], cdepu[,i])
  plot(taccu, caccu[,i], type = "l", las = 1, ylim = c(0, ymax),
       xlim = c(0, tf), las = 1, lwd = 2, xlab = "Time (days)",
       ylab = expression(paste("[Cd] (in ",mu,"g.",g-1," d.w.)")))
  abline(v = tc, lty = 2)
  lines(tdepu, cdepu[,i], xlim = c(tc, tf), lwd = 2)
}

# Simulation from the R-package `deSolve` ####

## Accumulation phase ####
# Write ODE system of equations
maccu <- function(t, y, parms){
  dy1 <- parms[1]*cx-parms[2]*y[1]+
    parms[10]*y[2]+parms[12]*y[3]+parms[14]*y[4]-
    (parms[9]+parms[11]+parms[13])*y[1]
  dy2 <- parms[3]*cx-parms[4]*y[2]+
    parms[9]*y[1]+parms[16]*y[3]+parms[18]*y[4]-
    (parms[10]+parms[15]+parms[17])*y[2]
  dy3 <- parms[5]*cx-parms[6]*y[3]+
    parms[11]*y[1]+parms[15]*y[2]+parms[20]*y[4]-
    (parms[12]+parms[16]+parms[19])*y[3]
  dy4 <- parms[7]*cx-parms[8]*y[4]+
    parms[13]*y[1]+parms[17]*y[2]+parms[19]*y[3]-
    (parms[14]+parms[18]+parms[20])*y[4]
  list(c(dy1, dy2, dy3, dy4))
}
initaccu <- rep(0, nc)
simuaccu <- ode(y = initaccu, times = taccu,
                func = maccu, parms = rates$median, method = "lsoda")

## Depuration phase ####
# Write ODE system of equations
mdepu <- function(t, y, parms){
  dy1 <- -parms[2]*y[1]+parms[10]*y[2]+
    parms[12]*y[3]+parms[14]*y[4]-
    (parms[9]+parms[11]+parms[13])*y[1]
  dy2 <- -parms[4]*y[2]+parms[9]*y[1]+
    parms[16]*y[3]+parms[18]*y[4]-
    (parms[10]+parms[15]+parms[17])*y[2]
  dy3 <- -parms[6]*y[3]+parms[11]*y[1]+
    parms[15]*y[2]+parms[20]*y[4]-
    (parms[12]+parms[16]+parms[19])*y[3]
  dy4 <- -parms[8]*y[4]+parms[13]*y[1]+
    parms[17]*y[2]+parms[19]*y[3]-
    (parms[14]+parms[18]+parms[20])*y[4]
  list(c(dy1, dy2, dy3, dy4))
}
initdepu <- simuaccu[nrow(simuaccu),2:ncol(simuaccu)]
simudepu <- ode(y = initdepu, times = tdepu,
                func = mdepu, parms = rates$median, method = "lsoda")

par(mar = c(4, 5, 0.1, 0.1), mfrow = c(2, 2))
for(i in 2:ncol(simuaccu)){
  ymax <- max(simuaccu[,i], simudepu[,i])/1000
  plot(simuaccu[,1], simuaccu[,i]/1000, type = "l", las = 1, ylim = c(0, ymax),
       xlim = c(0, tf), las = 1, lwd = 2, xlab = "Time (days)",
       ylab = expression(paste("[Cd] (in ",mu,"g.",g-1," d.w.)")))
  abline(v = tc, lty = 2)
  lines(simudepu[,1], simudepu[,i]/1000, lwd = 2)
}

