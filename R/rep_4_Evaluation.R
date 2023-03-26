
### Replicating:
# Step 5: fitting the model to the VE data, and
# Step 6: evaluating the model by comparing its predictions with the VS.


# Setup -------------------------------------------------------------------

library(mgcv)
library(parallel)
library(fst)
library(data.table)
library(tictoc)

source("R/setup.R")
source("R/step.R")
source("R/compare.R")
source("R/replication.R")

# nthreads <- detectCores(logical = FALSE) # no. of cores
nthreads <- 8L


# Load data and models -----------------------------------------------------

# The real data
load("data/data.RData")
observed <- data$dens

# The full model
load("data/gamm_full.RData")

# The selected model
load("data/gamm_selected.RData")

# Q-Q maps
load("data/qq_fit.RData")

# Spatiotemporal environmental data
ndata <- read.fst("data/ndata.fst")

ndata <- data.table(ndata, key = "id_year")
ndata[["fyear"]] <- factor(ndata[["year"]])
ndata[["plot_id"]] <- factor(levels(fit[["model"]][["plot_id"]])[1])
ndata[["observer_id"]] <- factor(levels(fit[["model"]][["observer_id"]])[1])
ndata[["eta"]] <- predict(fit, ndata, type = "link", exclude = c("s(plot_id)", "s(observer_id"))
ndata <- ndata[, c("id_year", "id", "year", "eta")]
id <- sort(unique(ndata[["id"]]))

# Original sampling scheme
ve <- data.table(data, key = "id_year")
ve <- ve[, c("id", "year") := NULL]
obs <- sort(unique(data[["observer_id"]]))
id_yr <- merge(ve, ndata[, "id_year"], by = "id_year")[, 1] # identifiers of id-year combination
save(id_yr, file = "data/id_yr.RData") # for later use

# Distribution properties
trans <- fit[["family"]][["linkinv"]] # inverse link
phi <- summary(fit, re.test = FALSE)$dispersion
xi <- fit[["family"]]$getTheta(TRUE)

gv <- quiet(gam.vcomp(fit))
sd_plot <- gv["s(plot_id)", "std.dev"] # sd for plot identifiers
sd_obs  <- gv["s(observer_id)", "std.dev"] # sd for observer identifiers
mu <- location(1, sd_obs)
sigma <- shape(1, sd_obs)

out <- c("plot_id", "observer_id", "fyear", "x", "y", "year") # variables to exclude
years <- sort(unique(fit[["model"]][["year"]]))
nd <- data.frame(year = years, fyear = factor(years), plot_id = levels(ve$plot_id)[1])

k <- 6 # smoothing parameter
fs <- update(formula(full), ve_prim ~ .)


# Replicating generation of VS --------------------------------------------

# Remark:
# The procedure is computationally demanding: every single iteration takes ~1/2 h.
# However, the simulation results can be downloaded from: https://doi.org/10.5281/zenodo.7732648, 
# (file 'replication.RData').

n <- 2

partial <- vector("list", n)
true_trend <- predicted_trend <- matrix(nrow = length(years), ncol = n)
cover <- numeric(n); cover[] <- NA
lambda <- numeric(n); lambda[] <- NA
lambda_se <- numeric(n); lambda_se[] <- NA
rep_ve <- matrix(NA, nrow = nrow(data), ncol = n)

# tic()

for (i in 1:n)
{
  cat(i, "\n")
  gc()
  
  # Generating VS from the model ----
  vs <- ndata
  id_r <- data.table(id = id, id_r = rnorm(length(id), 0, sd_plot))
  vs <- merge(vs, id_r, by = "id")
  vs[["etar"]] <- vs[["eta"]] + vs[["id_r"]]
  vs[["mur"]] <- trans(vs[["etar"]])
  vs[["vs"]] <- rTweedie(vs[["mur"]], xi, phi)
  vs_vs <- as.numeric(vs[["vs"]])

  # Calibrating VS ----
  for (j in 1:length(t)) vs_vs <- t[[j]](x = vs_vs)
  vs[["vs_prim"]] <- vs_vs
  vs <- vs[, c("id_year", "year", "vs", "vs_prim")]
  
  # Obtaining the VE data by sampling from the VS ----
  vei <- merge(ve, vs, by = "id_year")
  obs_r <- data.frame(observer_id = obs, obs_r = rlnorm(length(obs), mu, sigma))
  vei <- merge(vei, obs_r, by = "observer_id")
  vei[, ve_prim := vs_prim * obs_r]
  rep_ve[, i] <- merge(id_yr, vei, by = "id_year")[["ve_prim"]]
  
  # Step 4: Fitting the spatiotemporal GAMM model to virtual ecologist (VE) data ----
  sim_full <- bam(fs, data = vei, family = tw(), discrete = TRUE, nthreads = nthreads) # full model
  sim_fit <- backward(sim_full, p = 0.05) # backward elimination
  # summary(sim_fit, re.test = FALSE)
  
  # Evaluating the VE model by comparing its predictions with the VS ----
  
  # Response curves:
  pres_full <- part_res(sim_full, out)
  pres_fitted <- part_res(sim_fit, out)
  partial[[i]] <- list(full = pres_full, fitted = pres_fitted)
  
  # Population trend
  ft <- bam(ve_prim ~ year + s(fyear, bs = "re") + s(plot_id, bs = "re"), data = vei, family = tw, discrete = TRUE, nthreads = nthreads)
  pred <- predict(ft, nd, exclude = c("s(plot_id)"), type = "link", se = TRUE)
  pre <- data.frame(year = years, fit = pred$fit, se = pred$se.fit)
  pre <- transform(pre, lci = fit - 1.96 * se, uci = fit + 1.96 * se)
  pre[-1] <- exp(pre[-1])
  pre[-1] <- pre[-1] / pre[-1][1, 1]
  true_d <- tapply(vs$vs_prim, vs$year, mean, na.rm = TRUE)
  true_d <- true_d / true_d[1]
  true_trend[, i] <- true_d
  predicted_trend[, i] <- pre$fit
  su <- summary(ft, re.test = FALSE)
  lambda[i] <- su[["p.coeff"]][["year"]]
  lambda_se[i] <- su[["se"]][["year"]]
  cover[i] <- 100 * mean((true_d < pre$uci) & (true_d > pre$lci))
}

# to <- toc()
# (to$toc - to$tic) / 36 / n # no. of hours per 100 replications

save(partial, true_trend, predicted_trend, cover, lambda, lambda_se, rep_ve, file = "data/replication.RData")


# What next ---------------------------------------------------------------

# The data generated here can be used for evaluating different methodological aspects. 
# Go to the script 'rep_5_Analysis.R' for an exemplary analysis. 

