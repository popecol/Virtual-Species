
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

# List of predictors:
v <- c("Cropland", "Mosaic_cropland", "Grassland", "Urban", "Forest", "Precip_spring_lag", "Precip_summer_lag", "Precip_winter", "Precip_spring", "Tmax_summer_lag", "Tmin_spring_lag", "Tmin_summer_lag", "Tmin_winter", "Tmin_spring", "Elevation", "Roughness", "Wetness")

f <- formula(paste0("dens ~ ", paste0("s(", v, ", k = k, bs = 'cr')", collapse = " + ")))
f01 <- update(f, ". ~ . + s(plot_id, bs = 're') + s(observer_id, bs = 're') + s(fyear, bs = 're') + s(year, k = 10, bs = 'gp') + s(x, y, bs = 'gp', k = 30)")


# GAMM fitting ------------------------------------------------------------

nthreads <- detectCores(logical = FALSE) # no. of cores to be used in 'bam'

# These parameters control for the model smoothness:
k <- 6
gamma <- 3

# Restricting "power" parameter (p) of the Tweedie distribution
tweedie.profile(data$d ~ 1, do.plot = TRUE, verbose = 1, p.vec = seq(1.01, 1.99, len = 6))

# Fitting the full model:
full <- bam(f01, data, family = tw(a = 1.2, b = 1.8), discrete = TRUE, nthreads = nthreads, gamma = gamma)
summary(full, re.test = FALSE)

# Backward elimination 
fit <- backward(full, p = 0.05)

# The final model
summary(fit)

# Some diagnostics:
op <- par(mfrow = c(2, 2)); gam.check(fit, rep = 100); par(op)

# Response curves
op <- par(mar = c(5, 3, 1, 1))
plot.gam(fit, scale = 0, scheme = 2, pages = 1, seWithMean = TRUE, residuals = FALSE)
par(op)


# Saving models -----------------------------------------------------------

save(full, file = "data/gamm_full.RData", compress = "xz")
save(fit, file = "data/gamm_selected.RData", compress = "xz")
