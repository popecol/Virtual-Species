
### An exemplary analysis of the data from a computational experiment. 
# Step 6: evaluating the model by comparing its predictions with the VS.

# Here, we show how to:
# 1. Compare fitted response curves against the truth.
# 2. Assess the correctness of variable selection.
# 3. Compare fitted population trends against the 'true' trend.
# 4. Test, what is the effect of sampling effort on a trend estimation.


# Setup -------------------------------------------------------------------

library(mgcv)
# library(data.table)

source("R/setup.R")
source("R/compare.R")


# Load data and models -----------------------------------------------------

# The model used as a prototype for generating the VS (a reference)
load("data/gamm_selected.RData")
summary(fit, re.test = FALSE)

# Results of replicated generations of VS from the above model
load("data/replication.RData")

rep_fitted <- lapply(partial, function(x) x[["fitted"]])
rep_full <- lapply(partial, function(x) x[["full"]])
out <- c("plot_id", "observer_id", "fyear", "x", "y", "year") # variables to exclude
reference <- part_res(fit, out) # reference model partial plots (true response curves)
rep_reference <- vector("list", length(rep_fitted))
rep_reference <- lapply(rep_reference, function(x) x <- reference)


# Comparing response curves -----------------------------------------------

# Graphically (Figure 5)
plot_rep(reference, rep_fitted)


# Concordance of reference response curves (those that were used to generate the VS) 
# against fitted ones (those that were obtained from replicated VS generation). 
comp <- mapply(compare, rep_reference, rep_fitted, MoreArgs = list(plot.me = FALSE), SIMPLIFY = FALSE)

rs_coverage <- concordance(comp, "coverage")
rs_r <- concordance(comp, "r")
rs_rmse <- concordance(comp, "RMSE")
rs_mae <- concordance(comp, "MAE")


# Descriptive statistics --------------------------------------------------

# (me <- apply(rs_coverage, 2, median, na.rm = TRUE))
(me <- apply(rs_coverage, 2, mean, na.rm = TRUE))
sort(me)
mean(me)

op <- par(mfrow = c(2, 2))
plot(rs_coverage, xlab = "Coverage [%]")
plot(rs_r, xlab = "Correlation")
plot(rs_rmse, xlab = "RMSE")
plot(rs_mae, xlab = "MAE")
par(op)


# Correctness of variable selection ---------------------------------------

comp_full <- mapply(compare, rep_reference, rep_fitted, rep_full, MoreArgs = list(plot.me = FALSE), SIMPLIFY = FALSE)

vars_selected <- lapply(comp_full, function(x) x$vfit)
vars_reference <- lapply(comp_full, function(x) x$vref)

# No. of variables selected
nvars_selected <- sapply(vars_selected, sum)
range(nvars_selected)
mean(nvars_selected)
quantile(nvars_selected, c(0.025, 0.975))

tab <- mapply(table, vars_selected, vars_reference, SIMPLIFY = FALSE)

# [1, 1] true negatives
# [2, 1] false positives
# [1, 2] false negatives
# [2, 2] true positives

# False positives
FP <- sapply(tab, function(x) 100 * (x[2, 1] / sum(x)))
hist(FP)
range(FP)
mean(FP)
quantile(FP, c(0.025, 0.975))

# False negatives
FN <- sapply(tab, function(x) 100 * (x[1, 2] / sum(x)))
hist(FN)
range(FN)
mean(FN)
quantile(FN, c(0.025, 0.975))

# correct classification rate (%)
ccr <- sapply(tab, function(x) 100 * (x[1, 1] + x[2, 2]) / sum(x))
mean(ccr)
quantile(ccr, c(0.025, 0.975))

# Specificity (true negative rate)
tnr <- sapply(tab, function(x) 100 * (x[1, 1] / (x[1, 1] + x[2, 1])))
mean(tnr)
quantile(tnr, c(0.025, 0.975))

# Sensitivity (true positive rate)
tpr <- sapply(tab, function(x) 100 * (x[2, 2] / (x[2, 2] + x[1, 2])))
mean(tpr)
quantile(tpr, c(0.025, 0.975))


op <- par(mar = c(5, 5, 2, 2))
hist(ccr, xlab = "Correctly selected variables [%]", main = "", cex.lab = 1.3, cex.axis = 1.1, col = col4[2])
par(op)


# Population trends -------------------------------------------------------

year <- sort(unique(fit[["model"]][["year"]]))

# Trend estimates based on replicated VS generation
true_d <- apply(true_trend, 1, mean)

m <- apply(predicted_trend, 1, mean)
ci <- apply(predicted_trend, 1, quantile, probs = c(0.025, 0.975)) # empirical confidence intervals
lci <- ci[1, ]; uci <- ci[2, ]
ylim <- range(predicted_trend)

op <- par(mar = c(5, 5, 1, 1))
plot(year, true_d, xlim = c(2000, 2020), ylim = ylim, type = "n", xlab = "Year", ylab = "Relative abundance", cex.lab = 1.3)
abline(h = 1, col = "grey70")
matlines(year, predicted_trend, type = "l", lty = 1, col = col_shaded[2])
# lines(year, true_d, col = col[2])
lines(year, m, col = col[1], lty = 2, lwd = 2)
lines(year, lci, col = col[2])
lines(year, uci, col = col[2])
# lines(x, ref, col = col[1], lty = 2, lwd = 2)
legend(2000.5, 0.55, legend = c("95% of estimated", "True"), lty = 1:2, col = col[2:1], lwd = c(1, 2), bty = "n", cex = 1.1)
par(op)

# Trend coverage
mean(cover)
quantile(cover, c(0.025, 0.975))

op <- par(mar = c(5, 5, 2, 2))
hist(cover, xlab = "Trend coverage [%]", main = "", cex.lab = 1.3, cex.axis = 1.1, col = col4[2])
par(op)

# The 'true' trend
beta <- 100 * apply(true_trend, 2, function(x) lm(log(x) ~ year)[["coefficients"]][["year"]])
mean(beta)
quantile(beta, c(0.025, 0.975))

# The estimated trend
b <- 100 * apply(predicted_trend, 2, function(x) lm(log(x) ~ year)[["coefficients"]][["year"]])
mean(b)
quantile(b, c(0.025, 0.975))
range(beta)


# Sampling effort ---------------------------------------------------------

# Remark:
# The procedure below is computationally demanding.
# However, the simulation results can be downloaded from: https://doi.org/10.5281/zenodo.7732648, 
# (file 'samp_eff.RData').

nthreads <- 4

# The real data
load("data/data.RData")

# Unique identifiers
load("data/id_yr.RData")

years <- sort(unique(data$year))
plot_id <- levels(data[["plot_id"]])

dat <- merge(id_yr, data, by = "id_year")
f <- ve_prim ~ year + s(fyear, bs = "re") + s(plot_id, bs = "re") # formula

lem <- 10 # no. of different subsample size
# m <- round(seq(10, length(plot_id), length.out = lem))
m <- round(exp(seq(log(10), log(length(plot_id)), length.out = lem)))
m

n1 <- 1
n2 <- 1000

n <- n2 - n1 + 1 # no. of replications
samp_eff_lambda <- samp_eff_lambda_se <- samp_eff_cover <- matrix(NA, lem, n)
samp_eff_predicted_trend <- array(dim = c(length(years), lem, n))

pb <- txtProgressBar(min = 0, max = n, style = 3, char = "+", width = 50)
ii <- 0
for (i in n1:n2)
{
  ii <- ii + 1
  dati <- dat
  dati[["ve_prim"]] <- rep_ve[, i] # substitution with simulated VE data
  true_d <- true_trend[, i] # the true trend (of VS)
  
  for (j in 1:length(m))
  {
    repeat {
      idx <- sample(plot_id, m[j]) # random subsample of m[j] sites
      datj <- droplevels(dati[plot_id %in% idx, ])
      if(all(years %in% unique(datj$year))) break
    }
    ft <- bam(f, datj, family = tw, discrete = TRUE, nthreads = nthreads) # fit the trend model
    id_1 <- factor(levels(datj$plot_id)[1])
    nd <- data.frame(year = years, fyear = factor(years), plot_id = id_1)
    pred <- predict(ft, nd, exclude = c("s(plot_id)"), type = "link", se = TRUE)
    pre <- data.frame(year = years, fit = pred$fit, se = pred$se.fit)
    pre <- transform(pre, lci = fit - 1.96 * se, uci = fit + 1.96 * se)
    pre[-1] <- exp(pre[-1])
    pre[-1] <- pre[-1] / pre[-1][1, 1]
    
    su <- summary(ft, re.test = FALSE)
    samp_eff_predicted_trend[, j, ii] <- pre$fit
    samp_eff_lambda[j, ii] <- su[["p.coeff"]][["year"]]
    samp_eff_lambda_se[j, ii] <- su[["se"]][["year"]]
    samp_eff_cover[j, ii] <- 100 * mean((true_d < pre$uci) & (true_d > pre$lci))
  }
  setTxtProgressBar(pb, ii)
}
close(pb)

save(samp_eff_lambda, samp_eff_lambda_se, samp_eff_cover, samp_eff_predicted_trend, m, file = "data/samp_eff.RData")
