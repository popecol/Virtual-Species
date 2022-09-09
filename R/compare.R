
# A set of functions supporting comparison of response curves.
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
  y <- c("Cropland", "Mosaic cropland", "Grassland", "Urban", "Forest", "Precip spring lag", "Precip summer lag", "Precip winter", "Precip spring", "Tmax summer lag", "Tmin spring lag", "Tmin summer lag", "Tmin winter", "Tmin spring", "Elevation", "Roughness", "Wetness")
  
  order(match(x, y))
}


compare <- function(reference, fitted, full = NULL, plot.me = TRUE, trans = I, scale = -1)
{
  # Plots reference curves (i.e. those used to generate a virtual species) against fitted ones. 
  # Invisibly returns measures of fit and model selection accuracy.
  
  # Arguments:
  #   reference:  list produced by 'part_res' for the virtual species model
  #   fitted:     list produced by 'part_res' for the fitted model
  #   full:       optional list produced by 'part_res' (full model)
  #   plot.me:    plot result?
  #   trans:      function to apply before plotting
  #   scale:      -1 (default): the same y-axis scale for each plot; 0: different y axis for each plot
  
  # Returns a list with components:
  #   variable:   variable name
  #   vref:       logical vector (TRUE if a variable was used to generate a virtual species model)
  #   vfit:       logical vector (TRUE if a variable was selected in the final model)
  #   r:          correlation coeficient between the fitted and the true values of the response curves
  #   coverage:   coverage probability: % of the time that the 95% confidence interval of the fitted response curve contains the true value

  corr <- function(x, y) suppressWarnings(cor(x[["fit"]], y[["fit"]])) # correlation coefficient
  coverage <- function(x, y) 100 * mean((y[["fit"]] > x[["lci"]]) & (y[["fit"]] < x[["uci"]])) # coverage probability
  
  all_names <- if(is.null(full)) union(names(reference), names(fitted)) else names(full)
  id <- ord(all_names)
  all_names <- all_names[id]
  
  if(is.null(full)) {
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
  
  if(plot.me) {
    require(RColorBrewer)
    col <- brewer.pal(5, "RdYlGn")[c(1, 5)]
    col4 <- adjustcolor(col, alpha = 1/4)
    
    ylim <- range(sapply(fitt, function(x) range(x$lci, x$uci)))
    if(scale == -1) 
      fitt <- lapply(fitt, function(x) append(x, list(ylim = ylim))) else
        fitt <- lapply(fitt, function(x) append(x, list(ylim = range(x$lci, x$uci))))
    
    rows <- floor(sqrt(length(all_names)))
    cols <- ceiling(length(all_names) / rows)
    op <- par(mfrow = c(rows, cols), mar = c(5, 3, 1, 1))
    for (i in all_names)
    {
      with(fitt[[i]], plot(x, trans(fit), type = "n", ylim = trans(ylim), xlab = xlab, ylab = "", cex.lab = 1.3))
      with(fitt[[i]], polygon(c(x, rev(x)), c(trans(lci), trans(rev(uci))), border = NA, col = col4[2]))
      with(fitt[[i]], lines(x, trans(fit), col = col[2], lty = 1, lwd = 2))
      with(ref[[i]], lines(x, trans(fit), col = col[1], lwd = 1, lty = 2))
    }
    par(op)
  }
  
  return(invisible(list(variable = all_names, vref = vref, vfit = vfit, r = rr, coverage = cover)))
}
