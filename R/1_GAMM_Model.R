
### Step 1: Fitting the spatiotemporal GAMM model to real abundance data ###


# Setup -------------------------------------------------------------------

library(mgcv)
library(parallel)
library(tweedie)
library(statmod)

source("R/step.R")

# The data ----------------------------------------------------------------

load("data/data.RData")
summary(data)

# See the vignette 'VS.html' for a detailed description of the data.


# GAMM formula ------------------------------------------------------------

# Setting the range parameter for the temporal smooth
Matern <- function(d, rho) exp(-d / rho) * (1 + d / rho)
d <- seq(0, 5, length.out = 100)
plot(d, Matern(d, rho = 1), type = "l", xlab = "Distance [years]", ylab = "Correlation", cex.lab =1.3)
lines(d, Matern(d, rho = 0.5), col = "green")
lines(d, Matern(d, rho = 0.25), col = "red")
legend("topright", legend = paste("rho =", c(1, 0.5, 0.25)), col = c("black", "green", "red"), lwd = 2, bty = "n", cex = 1.1)

# List of predictors:
# dput(names(data))
v <- c("Cropland", "Mosaic_cropland", "Grassland", "Urban", "Forest", "Precip_spring_lag", "Precip_summer_lag", "Precip_winter", "Precip_spring", "Tmax_summer_lag", "Tmin_spring_lag", "Tmin_summer_lag", "Tmin_winter", "Tmin_spring", "Elevation", "Roughness", "Wetness")

f <- formula(paste0("dens ~ ", paste0("s(", v, ", k = k, bs = 'cr')", collapse = " + ")))
f01 <- update(f, ". ~ . + s(plot_id, bs = 're') + s(observer_id, bs = 're') + s(fyear, bs = 're') + s(year, k = 10, bs = 'gp', m = c(3, 0.5)) + s(x, y, k = 30, bs = 'gp')")


# GAMM fitting ------------------------------------------------------------

# nthreads <- detectCores(logical = FALSE) # no. of cores to be used in 'bam'
nthreads <- 8L

# These parameters control for the model smoothness:
k <- 6
gamma <- 3

# Full model
full <- bam(f01, data, family = tw(), discrete = TRUE, nthreads = nthreads, gamma = gamma)
summary(full, re.test = FALSE)

# Backward elimination 
fit <- backward(full, p = 0.05)

# The final model
summary(fit, re.test = FALSE)

# Some diagnostics:
op <- par(mfrow = c(2, 2)); gam.check(fit, rep = 100); par(op)

# Response curves
op <- par(mar = c(5, 3, 1, 1))
plot.gam(fit, scale = 0, scheme = 2, pages = 1, seWithMean = FALSE, residuals = FALSE)
par(op)


# Saving models -----------------------------------------------------------

save(full, file = "data/gamm_full.RData", compress = "xz")
save(fit, file = "data/gamm_selected.RData", compress = "xz")


# What next ---------------------------------------------------------------

# 1. To generate a single instance of a VS, go to the script '2_Virtual_Species.R'.

# 2. To replicate generating VS and VE data from the model, go to the script 'rep_2_Simulation.R' . 

