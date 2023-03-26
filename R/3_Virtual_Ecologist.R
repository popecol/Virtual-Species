
### Step 3: Obtaining the virtual ecologist data by sampling from the virtual species ###


# Setup -------------------------------------------------------------------

library(fst)
library(mgcv)
library(data.table)

source("R/setup.R")

# Load data and model -----------------------------------------------------

# The data
load("data/data.RData")
observed <- data$dens

# The model
load("data/gamm_selected.RData")

# The VS
vs <- read.fst("data/vs.fst")
vs <- vs[c("id_year", "id", "year", "mu", "sim_dr")]
vs <- data.table(vs, key = "id_year")


# Replicating the original sampling scheme --------------------------------

ve <- data.table(data, key = "id_year")
ve <- ve[, c("id", "year") := NULL]
ve <- merge(ve, vs, by = "id_year")
ve <- ve[complete.cases(ve), ]

# Observer error
obs <- sort(unique(ve[["observer_id"]]))
gv <- quiet(gam.vcomp(fit))
sd_obs  <- gv["s(observer_id)", "std.dev"] # sd for observer identifiers
mu <- location(1, sd_obs)
sigma <- shape(1, sd_obs)

set.seed(1)
obs_r <- data.frame(observer_id = obs, p = rlnorm(length(obs), mu, sigma))
# summary(obs_r); hist(obs_r$p, breaks = 50)
ve <- merge(ve, obs_r, by = "observer_id")
ve[, sim := sim_dr * p]
summary(ve$sim); hist(log1p(ve$sim), breaks = 50)


# Comparing observed vs. simulated ----------------------------------------

# Posterior predictive check of log-densities
pp_check(log1p(observed), log1p(ve$sim))

# Q-Q plot
qqplot(observed, ve$sim, xlab = "Observed", ylab = "Simulated"); grid(); abline(a = 0, b = 1, lty = 2)

# Comparing histograms
breaks <- seq(0, 5, 1/3)
hist_comp(log1p(ve$sim), log1p(observed), breaks = breaks, xname = "Population log-density")
legend(2, 1.5, legend = c("VE", "Observed"), lwd = 2, lty = c(NA, 1), col = col2[2:1], fill = c(col4[2], NA), border = c(col2[2], NA), bty = "n", cex = 1.3, merge = TRUE)


# Saving virtual ecologist data -------------------------------------------
save(ve, file = "data/ve.RData")


# What next ---------------------------------------------------------------

# Go to the script '4_Evaluation.R' for testing the efficiency of reconstructing 
# response curves and population trends. 

