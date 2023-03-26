
### A set of functions supporting replicated simulations of VS and VE data.
# Author: lechu@amu.edu.pl


rep_VSE <- function(n, vs_model, pred_data, map = NULL, cl = NULL, type = "txt")
{
  # This function simulates VS and VE data from the fitted GAMM model.  
  # It is intended to use for obtaining distributional characteristics of simulated data
  # and estimating transformations needed to calibrate these data.
  # 
  # The randomness is introduced at three levels: 
  #   1) while generating random intercepts for sampling plot identifiers, 
  #   2) when drawing a random sample from the Tweedie distribution, 
  #   3) when allocating random sampling error to individual observers. 
  
  # Arguments:
  #         n: number of replications
  #  vs_model: 'gam' or 'bam' model which was used to generate the VS
  # pred_data: data.table with the following variables:
  #               1) id_year: a unique identifier: combination of plot id and year (factor)
  #               2) plot_id: plot identifier (factor)
  #               3) observer_id: observer identifier (factor)
  #               4) year: year (integer)
  #               5) eta: predictions from the 'vs_model' on a link scale with 
  #                       both 'plot_id' and 'observer_id' excluded
  #       map: definition of the composite function to map the simulated VS; it consists
  #            of a list of calibrating functions which are applied sequentially, i.e.
  #            the first function is applied to the simulated VS, the returned result is 
  #            passed to the second function, etc. 
  #        cl: cluster object created by 'makeCluster' (if NULL, no cluster is used)
  #      type: type of the progress bar (use "none" to switch it off)
  # 
  # Returns: an array with simulated VS and VE data
  
  require(pbapply)
  
  op <- options("warn")
  options(warn = -1)
  
  family <- vs_model[["call"]][["family"]]
  if(family[[1]] == "Tweedie") p <- family[["p"]] else
    if(family[[1]] == "tw") p <- vs_model[["family"]]$getTheta(TRUE) else 
      stop("This function works with the Tweedie distribution only.")
  phi <- summary(vs_model, re.test = FALSE)$dispersion
  
  trans <- vs_model[["family"]][["linkinv"]] # inverse link
  id <- sort(unique(pred_data[["plot_id"]]))
  obs <- sort(unique(pred_data[["observer_id"]]))
  
  gv <- quiet(gam.vcomp(vs_model))
  sd_plot <- gv["s(plot_id)", "std.dev"] # sd for plot identifiers
  sd_obs  <- gv["s(observer_id)", "std.dev"] # sd for observer identifiers
  mu <- location(1, sd_obs)
  sigma <- shape(1, sd_obs)

  pbop <- pboptions(type = type)
  w <- pbreplicate(n, {
    id_r <- data.table(plot_id = id, id_r = rnorm(length(id), 0, sd_plot))
    dta <- merge(pred_data, id_r, by = "plot_id")
    dta[["etar"]] <- dta[["eta"]] + dta[["id_r"]]
    dta[["mur"]] <- trans(dta[["etar"]])
    dta[["vs"]] <- rTweedie(dta[["mur"]], p, phi)
    if(!is.null(map)) {
      vs <- dta[["vs"]]
      for (i in 1:length(map)) vs <- map[[i]](vs)
      dta[["vs"]] <- vs
    }
    obs_r <- data.table(observer_id = obs, obs_r = rlnorm(length(obs), mu, sigma))
    dta <- merge(dta, obs_r, by = "observer_id")
    dta[, ve := vs * obs_r]
    cbind(dta[["vs"]], dta[["ve"]])
  }, 
  cl = cl)
  
  pboptions(pbop)
  options(op)
  return(w)
}


qq_map <- function(source, target, ...)
{
  # Fits a quantile map
  
  # Arguments:
  #     source: source variable
  #     target: target variable
  #       ....: arguments passed to mgcv::gam
  
  # Returns: an object of class 'qq_map' that inherits from 'gam'.
  
  convert_2_list <- function(x) {
    # Helper fun: converts a matrix to a list.
    if(is.matrix(x))
      lapply(seq_len(nc), function(i) x[, i])
    else
      rep(list(x), nc)
  }
  
  qq <- function(x, y) list(qqplot(x, y, plot.it = FALSE))
  
  nc <- max(ncol(source), ncol(target), 1)
  source_list <- convert_2_list(source)
  target_list <- convert_2_list(target)
  
  qq_list <- mapply(qq, source_list, target_list)
  qx <- sapply(qq_list, function(i) i$x)
  qy <- sapply(qq_list, function(i) i$y)
  x <- apply(qx, 1, median)
  y <- apply(qy, 1, median)
  obj <- gam(y ~ s(x, ...), subset = x > 0)

  class(obj) <- c("qq_map", "gam")
  return(obj)
}



predict.qq_map <- function(obj, newdata, ...)
{
  # Predict method for the class 'qq_map'.

  mapped <- predict.bam(obj, newdata, ...)
  zeros <- newdata$x == 0
  mapped <- ifelse(zeros, 0, mapped)
  mapped <- threshold(mapped, 0)
  return(as.numeric(mapped))
}



rep_qqplot <- function(x, y, plot.it = TRUE, cl = NULL, type = "txt")
{
  # Q-Q plot for replicated data
  
  # Arguments:
  #       x: observed data
  #       y: replicated data (e.g. returned by the function 'rep_VSE')
  # plot.it: logical. Should the result be plotted?
  #      cl: cluster object created by 'makeCluster'
  #    type: type of the progress bar (use "none" to switch it off)
  # 
  # Invisibly returns a list with observed and simulated quantiles. 
  
  qqp <- function(y) qqplot(x = x, y = y, plot.it = FALSE)[["y"]]
  
  simq <- apply(y, 2, qqp) # simulated quantiles
  obsq <- qqplot(x, y, plot.it = FALSE)[["x"]] # observed quantiles
  
  oppb <- pboptions(type = type)
  q <- pbapply(simq, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE, cl = cl)
  pboptions(oppb)
  
  if(plot.it)
  {
    op <- par(mar = c(4.5, 4.5, 2, 2))
    plot(obsq, q[3, ], type = "n", ylim = range(q), xlab = "Observed quantiles", ylab = "Simulated quantiles", cex.lab = 1.3)
    grid()
    matlines(obsq, simq, col = col_shaded[2])
    lines(obsq, q[1, ], col = col[2], lty = 2)
    lines(obsq, q[3, ], col = col[2], lty = 2)
    # lines(obsq, q[2, ], col = col[2], lwd = 2) # median
    abline(a = 0, b = 1, lty = 2, col = col[1], lwd = 2) # equality line
    par(op)
  }
  invisible(list(x = obsq, y = simq))
}



rep_hist <- function(observed, simulated, FUN = mean, ref = NA, ...)
{
  # Plots a histogram of replicated data and adds a reference.
  
  # Arguments:
  #      observed: observed data
  #     simulated: simulated data
  #           FUN: a function to use as a summary statistic
  #           ref: a reference value (optional): will be used as a reference line instead of FUN(observed)
  #           ...: arguments passed to 'hist'

  if(is.na(ref)) ref <- FUN(observed)
  sim <- apply(simulated, 2, FUN)
  ci <- quantile(sim, c(0.025, 0.975))
  # title <- deparse(substitute(FUN))
  hist(sim, xlim = range(ref, sim), col = col4[2], border = col2[2], cex.lab = 1.3, ...)
  abline(v = ci, lty = 2, col = col[2])
  abline(v = ref, lwd = 2, lty = 2, col = col[1])
}

