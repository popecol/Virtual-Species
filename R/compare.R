
### A set of functions supporting comparison of response curves.
# Author: lechu@amu.edu.pl


part_res <- function(obj, ...) {
  UseMethod("part_res", obj)
}


part_res.gam <- function(obj, out = NULL, ...)
{
  # Calculates components needed to produce partial dependence plots. 
  # Arguments:
  #   obj:  fitted 'gam' or 'bam' object
  #   out:  vector of variables to omit
  #   ...:  additional arguments passed to 'plot.gam'
  # Value: named list
  
  require(mgcv)
  pg <- plot.gam(obj, select = 0, n = 100, ...)
  nam <- names(obj[["var.summary"]])
  id <- !nam %in% out
  nam <- nam[id]
  pg <- pg[id]
  names(pg) <- nam
  return(pg)
}


fill_missing_slots <- function(x, y)
{
  # Returns the list 'x' with slots added from the list 'y' that are missing in 'x'.
  # All supplemented slots contain zeros for 'fit' and 'se' components.
  
  all_names <- union(names(x), names(y))
  absent <- setdiff(all_names, names(x))
  pg_absent <- y[absent]
  zero <- vector("list", length(absent))
  names(zero) <- absent
  zero <- lapply(zero, append, list(fit = rep(0, 100), se = rep(0, 100)))
  pg_absent <- modifyList(pg_absent, zero)
  modifyList(x, pg_absent)
}


ord <- function(x)
{
  # Helper function to set a consistent order of variables
  # Argument: character vector (variable names)
  # Value:    indexes of x ordered according to y (y is defined below)
  
  # The required order:
  y <- c("Cropland", "Mosaic_cropland", "Grassland", "Urban", "Forest", "Precip_spring_lag", "Precip_summer_lag", "Precip_winter", "Precip_spring", "Tmax_summer_lag", "Tmin_spring_lag", "Tmin_summer_lag", "Tmin_winter", "Tmin_spring", "Elevation", "Roughness", "Wetness")
  
  order(match(x, y))
}


compare <- function(reference, fitted, full = NULL, plot.me = TRUE, trans = I, scale = -1)
{
  # Plots reference curves (i.e. those used to generate a virtual species) against fitted ones. 
  # Invisibly returns measures of fit, a list of predictor variables and information needed to calculate the model selection accuracy.
  
  # Arguments:
  #   reference:  list produced by 'part_res' for the virtual species model
  #   fitted:     list produced by 'part_res' for the fitted model
  #   full:       optional list produced by 'part_res' for the full model
  #   plot.me:    plot result?
  #   trans:      function to apply before plotting
  #   scale:      -1 (default): the same y-axis scale for each plot; 0: different y axis for each plot
  
  # Returns a list with components:
  #   variable:   variable name
  #   vref:       logical vector (TRUE if a variable was used to generate a virtual species model)
  #   vfit:       logical vector (TRUE if a variable was selected in the final model)
  #   coverage:   coverage probability: % of the time that the 95% confidence interval of the fitted response curve contains the true value
  #   r:          correlation coeficient between the fitted and the true values of the response curves
  #   RMSE:       root mean squared error (fitted vs. true)
  #   MAE:        mean absolute deviation (fitted vs. true)
  
  coverage <- function(x, y) 100 * mean((y[["fit"]] > x[["lci"]]) & (y[["fit"]] < x[["uci"]]))
  corr <- function(x, y) suppressWarnings(cor(x[["fit"]], y[["fit"]]))
  rmse <- function(x, y) sqrt(mean((y[["fit"]] - x[["fit"]])^2))
  mae <- function(x, y) mean(abs(y[["fit"]] - x[["fit"]]))
  
  all_names <- if(is.null(full)) union(names(reference), names(fitted)) else names(full)
  id <- ord(all_names)
  all_names <- all_names[id]
  
  if(is.null(full)) 
  {
    ref <- fill_missing_slots(reference, fitted)
    fitt <- fill_missing_slots(fitted, reference)
  } else {
    ref <- fill_missing_slots(reference, full)
    fitt <- fill_missing_slots(fitted, full)
  }
  
  fitt <- lapply(fitt, function(x) append(x, list(lci = x$fit - 1.96 * x$se, uci = x$fit + 1.96 * x$se)))
  ref <- ref[ord(names(ref))]
  fitt <- fitt[ord(names(fitt))]
  vref <- all_names %in% names(reference)
  vfit <- all_names %in% names(fitted)
  
  cover <- mapply(coverage, x = fitt, y = ref)
  rr <- mapply(corr, x = fitt, y = ref)
  rms <- mapply(rmse, x = fitt, y = ref)
  ma <- mapply(mae, x = fitt, y = ref)
  
  if(plot.me)
  {
    fitt <- lapply(fitt, function(x) append(x, list(xlab_spaced = gsub("_", " ", x$xlab))))
    ylim <- range(sapply(fitt, function(x) range(x$lci, x$uci)))
    if(scale == -1) 
      fitt <- lapply(fitt, function(x) append(x, list(ylim = ylim))) else
        fitt <- lapply(fitt, function(x) append(x, list(ylim = range(x$lci, x$uci))))
    
    rows <- floor(sqrt(length(all_names)))
    cols <- ceiling(length(all_names) / rows)
    op <- par(mfrow = c(rows, cols), mar = c(5, 3, 1, 1))
    for (i in all_names)
    {
      with(fitt[[i]], plot(x, trans(fit), type = "n", ylim = trans(ylim), xlab = xlab_spaced, ylab = "", cex.lab = 1.3))
      with(fitt[[i]], polygon(c(x, rev(x)), c(trans(lci), trans(rev(uci))), border = NA, col = col4[2]))
      with(fitt[[i]], lines(x, trans(fit), col = col[2], lty = 1, lwd = 2))
      with(ref[[i]], lines(x, trans(fit), col = col[1], lwd = 1, lty = 2))
    }
    par(op)
  }
  return(invisible(list(variable = all_names, vref = vref, vfit = vfit, r = rr, coverage = cover, RMSE = rms, MAE = ma)))
}



concordance <- function(comp, measure = c("coverage", "r", "RMSE", "MAE"))
{
  # Calculates concordance measures between the fitted and the true values of the response curves. 
  # This function is intended to use when generating a multiple instances of a virtual species. 
  
  # Arguments:
  #   comp:     list produced by 'mapply' when applying a 'compare' function to reference and fiited lists.
  #   measure:  type of concordance measure (coverage probability, correlation coefficient, mean square root error or mean absolute error)

  # Returns: an object of class 'concordance' (which inherits from a 'data frame'). 
  # Every column in this frame refers to a predictor variable and rows represent replications of a VS.

  all_names <- unique(unlist(sapply(comp, function(x) x[["variable"]])))
  nvars <- length(all_names)
  
  d <- lapply(comp, function(x) x[[measure]])
  
  w <- vector("list", nvars)
  names(w) <- all_names
  for (i in 1:nvars)
  {
    v <- all_names[i]
    idx <- sapply(d, function(x) is.numeric(try(x[[v]], silent = TRUE)))
    w[[i]] <- sapply(d[idx], function(x) x[[v]])
  }
  idx <- names(w) %in% names(reference)
  w <- w[idx]
  w <- do.call(cbind, w)
  if(measure == "coverage") w <- ifelse(w > 0, w, NA)
  w <- data.frame(w)
  class(w) <- c("concordance", "data.frame")
  return(w)
}


plot.concordance <- function(obj, FUN = median, ...)
{
  # A 'plot' method for the class 'concordance'.  
  # Calculates quantiles (95% and 50%) for statistics measuring similarity between the
  # fitted and the reference response curves. Then, plots the result. 
  # Invisibly returns a data frame containing respective quantiles for every predictor variable.
  
  # Arguments:
  #   obj:   object of the class 'concordance'. 
  #   FUN:   type of an average (e.g. median or mean) passed to 'apply'. 
  #   ...:   additional arguments passed to 'plot' function.

  q <- data.frame(t(apply(obj, 2, quantile, probs = c(0.025, 0.25, 0.75, 0.975), na.rm = TRUE)))
  m <- apply(obj, 2, FUN, na.rm = TRUE)
  
  xlim <- range(q, na.rm = TRUE)
  nam <- names(obj)
  nam <- nam[ord(nam)]
  nv <- length(nam)
  nam_spaced = gsub("_", " ", nam)
  
  op <- par(mar = c(5, 8, 2, 2))
  plot(m, 1:nv, type = "n", xlim = xlim, ylim = c(0.5, 10.5), ylab = "", cex.lab = 1.3, yaxt = "n", ...)
  axis(2, at = 1:nv, labels = rev(nam_spaced), las = 1, cex.axis = 1.1)
  grid(nx = NULL, ny = NA)
  abline(h = 1:nv, col = "lightgray", lty = "dotted")
  segments(q$X2.5., 1:nv, q$X97.5., col = col2[2], lwd = 2)
  segments(q$X25., 1:nv, q$X75., col = col2[2], lwd = 4)
  points(m, 1:nv, col = "white", cex = 1.7, pch = 16)
  points(m, 1:nv, col = col[2], cex = 1.3, pch = 16)
  par(op)
  
  return(invisible(q))
}


plot_rep <- function(reference, rep_fitted)
{
  # Plots response curves for all predictor variables from a replicated instances of a VS (with the true overlaid).
  
  # Arguments:
  #   reference:   object produced by the 'part_res' function for the reference model.
  #   rep_fitted:  as above, but in a form of a list containing results of multiple calls to 'part_res' for fitted models. 
  # Returns: nothing, just a plot. 
  
  all_names <- names(reference)
  id <- ord(all_names)
  all_names <- all_names[id]
  
  # Setting the layout
  rows <- floor(sqrt(length(all_names)))
  cols <- ceiling(length(all_names) / rows)
  op <- par(mfrow = c(rows, cols), mar = c(5, 3, 1, 1))
  
  for (i in all_names)
  {
    is_i <- sapply(rep_fitted, function(x) grep(i, names(x)))
    idx <- which(sapply(is_i, length) > 0)[1]
    x <- rep_fitted[[idx]][[i]][["x"]]
    y <- sapply(rep_fitted, function(x) x[[i]][["fit"]])
    if(is.list(y))
    {
      idx <- !sapply(y, is.null)
      y <- y[idx]
      y <- do.call(cbind, y)
    }
    
    ref <- reference[[i]][["fit"]]
    # m <- apply(y, 1, mean)
    ci <- data.frame(t(apply(y, 1, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))))
    ylim <- range(y)
    xlab <- gsub("_", " ", i)
    
    plot(x, ref, type = "n", xlab = xlab, ylab = "", cex.lab = 1.3, ylim = ylim)

    yp <- y[, sample(1:dim(y)[2], 250)]
    
    matlines(x, yp, type = "l", lty = 1, col = col_shaded[2])
    
    lines(x, ci$X2.5., col = col[2])
    lines(x, ci$X97.5., col = col[2])
    # lines(x, ci$X25., col = col[2])
    # lines(x, ci$X75., col = col[2])
    # lines(x, ci$X50., col = col[2])
    # lines(x, m, col = col[2])
    lines(x, ref, col = col[1], lty = 2, lwd = 2)
  }
  par(op)
}


quiet <- function(x) {
  # Helper fun: returns function output without printed messages
  # Author: Hadley Wickham
  
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}
