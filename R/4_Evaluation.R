
### Step 5: Fitting the spatiotemporal GAMM model to virtual ecologist (VE) data ###

### Step 6: Evaluating the VE model by comparing its predictions with the VS data ###


# Setup -------------------------------------------------------------------

library(mgcv)
library(parallel)
library(fst)
library(RColorBrewer)

source("R/step.R")
source("R/compare.R")

# nthreads <- detectCores(logical = FALSE) # no. of cores
nthreads <- 8L


# Load data and model -----------------------------------------------------

# The full model
load("data/gamm_full.RData")

# The selected model
load("data/gamm_selected.RData")

# The VE data
load("data/ve.RData")


# Fitting the original model to the VE data -------------------------------

k <- 6

# Full model
(fs <- update(formula(full), sim ~ .))

sim_full <- bam(fs, data = ve, family = tw(), discrete = TRUE, nthreads = nthreads)
summary(sim_full, re.test = FALSE)

### Backward elimination 
sim_fit <- backward(sim_full, p = 0.05)
summary(sim_fit, re.test = FALSE)

op <- par(mfrow = c(2, 2)); gam.check(sim_fit, rep = 100); par(op)
# summary(sim_fit)

op <- par(mar = c(5, 3, 1, 1))
plot.gam(sim_fit, scale = 0, scheme = 2, pages = 1, seWithMean = FALSE, residuals = FALSE)
par(op)


# Comparing response curves -----------------------------------------------

out <- c("plot_id", "observer_id", "fyear", "x", "y", "year") # variables to exclude from the comparison

# Calculating partial residuals:
reference <- part_res(fit, out) # reference model, i.e. the model used to generate a virtual species
full <- part_res(sim_full, out) # full model fitted to the VE data
fitted <- part_res(sim_fit, out) # final model fitted to the VE data

# Comparing response curves that were used to generate the VS (dashed red lines)
# with fitted curves estimated based on VE data (green).
m <- compare(reference, fitted, scale = 0)
m$coverage # coverage percentages
# m$r # corelation coefficients
# m$RMSE # root mean square error
# m$MAE # mean absolute error

# Correctness of variable selection:
m <- compare(reference, fitted, full)
table(selected = m$vfit, reference = m$vref)


# Comparing population trends ---------------------------------------------

# The 'true' trend:
vs <- read.fst("data/vs.fst")
true_d <- tapply(vs$sim_dr, vs$year, mean, na.rm = TRUE)
true_d <- true_d / true_d[1] # relative abundance
years <- sort(unique(vs$year))
plot(years, true_d, type = "b", xlab = "Year", ylab = "Population abundance", main = "The 'true' trend")

# The trend estimated based on VE data
ft <- bam(sim ~ year + s(fyear, bs = "re") + s(plot_id, bs = "re"), data = ve, family = tw, discrete = TRUE, nthreads = nthreads)
summary(ft, re.test = FALSE)

nd <- data.frame(year = years, fyear = factor(years), plot_id = levels(ve$plot_id)[1])
pred <- predict(ft, nd, exclude = c("s(plot_id)"), type = "link", se = TRUE)
pre <- data.frame(year = years, fit = pred$fit, se = pred$se.fit)
pre <- transform(pre, lci = fit - 1.96 * se, uci = fit + 1.96 * se)
pre[-1] <- exp(pre[-1]) # transforming to the response scale
pre[-1] <- pre[-1] / pre[-1][1, 1] # transforming to relative abundance

op <- par(mar = c(5, 5, 1, 1))
plot(years, true_d, ylim = range(0.35, true_d, pre[c("lci", "uci")]), type = "n", xlab = "Year", ylab = "Relative abundance", cex.lab = 1.3)
abline(h = 1, col = "grey70")
polygon(c(years, rev(years)), c(pre$lci, rev(pre$uci)), border = NA, col = col4[2])
lines(fit ~ years, pre, col = col[2], lty = 1, lwd = 2)
lines(years, true_d, col = col[1], lwd = 2, lty = 2)
legend(2001, 0.5, legend = c("Estimated", "True"), lty = 1:2, col = col[2:1], lwd = 2, bty = "n", cex = 1.1)
par(op)

cor(pre$fit, true_d) # correlation
100 * mean((true_d < pre$uci) & (true_d > pre$lci)) # coverage probability

# Saving models -----------------------------------------------------------

save(sim_full, file = "data/sim_full.RData", compress = "xz")
save(sim_fit, file = "data/sim_selected.RData", compress = "xz")


# What next ---------------------------------------------------------------

# The above procedure can be run many times to obtain distributions of parameters. 
# Go to the script 'rep_2_Simulation.R' to perform this analysis. 

