
### Step2: Generating a dynamic virtual species from a spatiotemporal model ###

# Setup -------------------------------------------------------------------

library(fst)
library(mgcv)
library(terra)
library(data.table)
library(RColorBrewer)

source("R/maps.R")

# The fitted model
load("data/gamm_selected.RData")
summary(fit, re.test = FALSE)


# Simulating from the model -----------------------------------------------

# Load spatiotemporal environmental data
# 313368 km^2 * 20 years * 17 environmental variables

vs <- read.fst("data/ndata.fst")
vs <- data.table(vs, key = "id_year")
vs[["fyear"]] <- factor(vs[["year"]])

# Need to provide something for 'predict.bam' with 'discrete=TRUE' to work
vs[["plot_id"]] <- factor(levels(fit[["model"]][["plot_id"]])[1])
vs[["observer_id"]] <- factor(levels(fit[["model"]][["observer_id"]])[1])

# Calculate expected values (on a linear predictor scale)
vs[["eta"]] <- predict(fit, vs, type = "link", exclude = c("s(plot_id)", "s(observer_id"))

# Simulate random intercepts for every square in the study area
sd_plot_id <- gam.vcomp(fit)["s(plot_id)", "std.dev"] # estimated SD for random intercepts
id <- sort(unique(vs[["id"]]))
set.seed(123)
id_r <- data.table(id = id, id_r = rnorm(length(id), 0, sd_plot_id))
# summary(id_r); hist(id_r$id_r)
vs <- merge(vs, id_r, by = "id")

# Add random intercepts to expected values
vs[["etar"]] <- vs[["eta"]] + vs[["id_r"]]

# Transform to the response scale
trans <- fit[["family"]][["linkinv"]]
vs[["mu"]] <- trans(vs[["eta"]])
vs[["mur"]] <- trans(vs[["etar"]])

# Random draws from the process distribution
# The estimated sale parameter (phi):
# phi <- sum((fit$weights * fit$residuals^2)[fit$weights > 0]) / fit$df.residual # Pearson method
phi <- fit$deviance / fit$df.residual # using the deviance method
p <- fit$family$getTheta(TRUE) # the estimated shape parameter
set.seed(123)
vs[["sim_dr"]] <- rTweedie(vs[["mur"]], p, phi)

# Maps --------------------------------------------------------------------

maps(vs, "mu")
maps(vs, "mur")
maps(vs, "sim_dr")
yr <- 2010; maps(vs[vs$year == yr, ], "mu", fact = 1, main = yr)


# Saving virtual species data ---------------------------------------------

vs <- vs[, c("id_year", "id", "year", "fyear", "x", "y", "mu", "mur", "sim_dr")]
write.fst(vs, "data/vs.fst", compress = 100)
