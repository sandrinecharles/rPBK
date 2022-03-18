# General entries ####
# Empty working directory
rm(list = ls())
# Load required library, if installed
library(tidyverse)
library(crosstalk)
library(deSolve)


# Load parameter values ####
# For Cd
# rates <- read.table("param4comp-gestin-Cd.txt", header = TRUE)
# For Hg
rates <- read.table("data-raw/param4comp-gestin-Hg.txt", header = TRUE)
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
for(i in 2:length(taccu)){
  intJ <- diag(x = (exp(X$values * taccu[i]) - 1) / X$values)
  intE <- P %*% intJ %*% solve(P)
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
       ylab = expression(paste("[Hg] (in ",mu,"g.",g^-1," d.w.)")))
       # ylab = expression(paste("[Cd] (in ",mu,"g.",g^-1," d.w.)")))
  abline(v = tc, lty = 2)
  lines(tdepu, cdepu[,i], xlim = c(tc, tf), lwd = 2)
}
