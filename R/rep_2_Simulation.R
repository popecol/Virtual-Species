
### Repeating Step 2 (simulating from the model) and Step 3 (virtual sampling) ###


# Setup -------------------------------------------------------------------

library(mgcv)
library(data.table)
library(fst)
library(parallel)
library(pbapply)

# Custom functions supporting replicated simulations
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
summary(fit, re.test = FALSE)

# Spatiotemporal environmental data
ndata <- read.fst("data/ndata.fst")


# Prepare data for simulations --------------------------------------------

ndat <- data.table(ndata, key = "id_year")
ndat[["fyear"]] <- factor(ndat[["year"]])
ndat[["plot_id"]] <- factor(levels(fit[["model"]][["plot_id"]])[1])
ndat[["observer_id"]] <- factor(levels(fit[["model"]][["observer_id"]])[1])
ndat[["eta"]] <- predict(fit, ndat, type = "link", exclude = c("s(plot_id)", "s(observer_id"))
ndat <- ndat[, c("id_year", "year", "eta")]
dat <- data.table(data, key = "id_year")
dat <- dat[, c("id", "year") := NULL]
pred_data <- merge(dat[, c("id_year", "plot_id", "observer_id")], ndat, by = "id_year")
# load("data/pred_data.RData")

observed <- data$dens
summary(observed); hist(log1p(observed))

clusterExport(cl, c("observed", "predict.qq_map", "threshold"))


# Simulate from the model ------------------------------------------------
# First run, to estimate calibration functions.

set.seed(1)
vse <- rep_VSE(n = 1000, fit, pred_data)
# load("data/vse.RData")
dim(vse)

vs_sim <- vse[, 1, ]
ve_sim <- vse[, 2, ]

prevalence(observed)
prevalence(vs_sim)
prevalence(ve_sim)

# Comparing simulated VE (untransformed) and observed data (takes a while)
op <- par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2, 2))
rep_qqplot(observed, ve_sim, cl = cl) # Q-Q plot
rep_hist(observed, ve_sim, mean, breaks = 30, main = "", xlab = "Mean")
rep_hist(observed, ve_sim, sd, breaks = 30, main = "", xlab = "Standard deviation")
rep_hist(observed, ve_sim, prevalence, breaks = 20, main = "", xlab = "Prevalence")
par(op)

stopCluster(cl)


# Saving simulation results -----------------------------------------------

save(pred_data , file = "data/pred_data.RData")
save(vse, file = "data/vse.RData")


# What next ---------------------------------------------------------------

# Go to the script 'rep_3_Calibration.R' to fit quantile-quantile maps. 

