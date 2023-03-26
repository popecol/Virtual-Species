
### Step 4: Fitting calibration functions and mapping frequency distributions ###


# Setup -------------------------------------------------------------------

library(mgcv)
library(data.table)
library(fst)
library(parallel)
library(pbapply)

source("R/setup.R")
source("R/replication.R")

# ncores <- detectCores(logical = FALSE)
ncores <- 4L
cl <- makeCluster(ncores)

clusterEvalQ(cl, {
  library(mgcv)
  library(data.table)
})


# Load data and models -----------------------------------------------------

# The real data
load("data/data.RData")

# The selected model
load("data/gamm_selected.RData")

# The data needed for simulation
load("data/pred_data.RData")

# The simulated data
load("data/vse.RData")
vs_sim <- vse[, 1, ]
ve_sim <- vse[, 2, ]

observed <- data$dens
summary(observed); hist(log1p(observed))


# Fitting Q-Q maps --------------------------------------------------------

# Virtual sampling: VS -> VE
s_hat <- qq_map(vs_sim, ve_sim, m = 1)
summary(s_hat)
xq <- seq(0, max(vs_sim), len = 100)
yq <- predict(s_hat, list(x = xq)); summary(yq)
plot(xq, yq, type = "l", xlab = "VS", ylab = "VE"); grid(); abline(a = 0, b = 1, lty = 2)
# plot(xq[-1], diff(yq), type = "l", xlab = "VS", ylab = "First derivative"); abline(h = 0, lty = 3)

# Mapping simulated VE data using observed data as a target: VE -> VE'
h <- qq_map(ve_sim, observed, m = 1)
summary(h)
xq <- seq(0, max(ve_sim), len = 100)
yq <- predict(h, list(x = xq)); summary(yq)
plot(xq, yq, type = "l", xlab = "VE", ylab = "Observed"); grid(); abline(a = 0, b = 1, lty = 2)
# plot(xq[-1], diff(yq), type = "l", xlab = "VE", ylab = "First derivative"); abline(h = 0, lty = 3)

# Inversed virtual sampling: VE -> VS
inv_s_hat <- qq_map(ve_sim, vs_sim, m = 1)
summary(inv_s_hat)
xq <- seq(0, max(ve_sim), len = 100)
yq <- predict(inv_s_hat, list(x = xq)); summary(yq)
plot(xq, yq, type = "l", xlab = "VE", ylab = "VS"); grid(); abline(a = 0, b = 1, lty = 2)
# plot(xq[-1], diff(yq), type = "l", xlab = "VE", ylab = "First derivative"); abline(h = 0, lty = 3)

# A list defining a composite function t: VS -> VS'
t <- list(
  s_hat = function(x) predict(s_hat, list(x = x)),
  f = function(x) {
    hh <- predict(h, list(x = x))
    q <- quantile(hh, 1 - prevalence(observed))
    threshold(hh, q)
  },
  inv_s_hat = function(x) predict(inv_s_hat, list(x = x))
)


# Simulate from the model -----------------------------------------------
# Second run, to apply calibration functions.

clusterExport(cl, c("pred_data", "observed", "predict.qq_map", "threshold", "prevalence", "s_hat", "h", "inv_s_hat", "t"))

set.seed(1)
vse_prim <- rep_VSE(n = 1000, fit, pred_data, map = t, cl = cl)
# load("data/vse_prim.RData")
dim(vse_prim)

vs_prim <- vse_prim[, 1, ]
ve_prim <- vse_prim[, 2, ]

prevalence(observed)
prevalence(ve_prim)


# Comparing simulated VE and observed data
ma <- max(max(log1p(observed)), max(log1p(ve_prim))) + 0.2
breaks <- seq(0, ma, 1/3)
op <- par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2, 2))
hist_comp(log1p(ve_prim), log1p(observed), breaks = breaks, xname = "Population log-density")
rep_qqplot(observed, ve_prim, cl = cl) # Q-Q plot
rep_hist(observed, ve_prim, mean, breaks = 30, main = "", xlab = "Mean")
rep_hist(observed, ve_prim, sd, breaks = 20, main = "", xlab = "Standard deviation")
# rep_hist(observed, ve_prim, prevalence, breaks = 30, main = "", xlab = "Prevalence")
par(op)


# Saving calibrating functions and simulation results ---------------------

save(vse_prim, file = "data/vse_prim.RData")
save(s_hat, h, inv_s_hat, t, file = "data/qq_fit.RData")

stopCluster(cl)


# What next ---------------------------------------------------------------

# Go to the script 'rep_4_Evaluation.R' to evaluate replicated VE data. 

