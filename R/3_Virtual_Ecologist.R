
### Step 3: Obtaining the virtual ecologist data by sampling from the virtual species ###

# Both original sampling scheme and observational errors are reconstructed. 


# Setup -------------------------------------------------------------------

library(fst)
library(mgcv)
library(data.table)

source("R/pp_check.R")

# Load data and model -----------------------------------------------------

# The data
load("data/data.RData")

# The model
load("data/gamm_selected.RData")

# The VS
vs <- read.fst("data/vs.fst")
vs <- vs[c("id_year", "id", "year", "mu", "mur", "sim_dr")]
vs <- data.table(vs, key = "id_year")


# Replicating the original sampling scheme --------------------------------

ve <- data.table(data, key = "id_year")
ve <- ve[, c("id", "year") := NULL]
ve <- merge(ve, vs, by = "id_year")
ve <- ve[complete.cases(ve), ]

# Observer error
gv <- gam.vcomp(fit)
sd_id_obs <- gv["s(observer_id)", "std.dev"] # estimated SD for observer random intercepts
obs <- sort(unique(ve[["observer_id"]]))
set.seed(123)
obs_r <- data.frame(observer_id = obs, obs_r = rlnorm(length(obs), 0, sd_id_obs))
# summary(obs_r); hist(obs_r$obs_r)
ve <- merge(ve, obs_r, by = "observer_id")
ve[, sim := sim_dr * obs_r]
# summary(ve)


# Comparing observed vs. simulated ----------------------------------------

# Posterior predictive check of log-densities
pp_check(log1p(data$dens), log1p(ve$mu))
pp_check(log1p(data$dens), log1p(ve$mur))
pp_check(log1p(data$dens), log1p(ve$sim_dr))
pp_check(log1p(data$dens), log1p(ve$sim))

# Distributions
breaks <- seq(0, 5, 0.5)
cex.lab <- 1.7
col <- adjustcolor("black", alpha = 0.15)

op <- par(mfrow = c(1, 3), mar = c(5, 5, 1, 1))
hist(log1p(data$dens), breaks = breaks, freq = FALSE, ylim = c(0, 1.1), main = "", xlab = "Observed log-densities", ylab = "Probability density", cex.lab = cex.lab, col = col)
hist(log1p(ve$sim), breaks = breaks, freq = FALSE, ylim = c(0, 1.1), main = "", xlab = "Simulated log-densities", ylab = "Probability density", cex.lab = cex.lab, col = col)
plot(log1p(data$dens), log1p(ve$sim), asp = 1, xlab = "Observed log-densities", ylab = "Simulated log-densities", cex.lab = cex.lab, col = col)
grid(); abline(a = 0, b = 1, lty = 2)
par(op)


# Saving virtual ecologist data -------------------------------------------

save(ve, file = "data/ve.RData")

