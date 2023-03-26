
### Setup and a set of helper functions.
# Author: lechu@amu.edu.pl

require(RColorBrewer)
col <- brewer.pal(5, "RdYlGn")[c(1, 5)]
col2 <- adjustcolor(col, alpha = 1/2)
col4 <- adjustcolor(col, alpha = 1/4)
col10 <- adjustcolor(col, alpha = 1/10)
col_shaded <- adjustcolor(col, alpha = 1/35)


prevalence <- function(x) sum(x > 0) / length(x)


threshold <- function(x, y)
{
  # Assignes zeros to all values in x that are not greater than y.
  x[x <= y] <- 0
  x
}


location <- function(m, s) {
  # Calculates the mean of the log-normal distribution.
  # Arguments:
  #   m: arithmetic mean
  #   s: arithmetic standard deviation
  
  log(m^2 / sqrt(s^2 + m^2))
} 


shape <- function(m, s) {
  # Calculates the standard deviation of the log-normal distribution.
  # Arguments:
  #   m: arithmetic mean
  #   s: arithmetic standard deviation
  
  sqrt(log(1 + (s^2 / m^2)))
}


hist_comp <- function(x, y, xname = NA, cex.lab = 1.3, ...)
{
  # Plots a histogram of x and adds an outline of y as a reference.
  # Can be used to visually compare histograms.
  
  # Arguments:
  #     x: a vector of values for which the histogram is desired
  #     y: a vector of values for the reference histogram
  # xname: variable name (optional)
  #   ...: arguments passed to 'hist'
  
  hy <- hist(y, plot = FALSE, ...)
  mid <- hy$mids
  mid <- mid - (mid[2] - mid[1]) / 2
  ydens <- hy$density
  mid <- c(0, mid); ydens <- c(0, ydens)
  
  hx <- hist(x, plot = FALSE, ...)
  # hx <- hist(x, plot = FALSE)
  xdens <- hx$density
  hx$breaks = hy$breaks
  if(!is.na(xname)) hx$xname <- xname
  
  plot(hx, freq = FALSE, ylim = range(xdens, ydens), col = col4[2], border = col2[2], main = "", ylab = "Probability density", cex.lab = cex.lab)
  lines(mid, ydens, type = "s", lwd = 2, col = col2[1])
}


quiet <- function(x) {
  # Helper fun: returns function output without printed messages
  # Author: Hadley Wickham
  
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}


maps <- function(newdata, vname, range = NULL, bias = 2, probs = 0.9999, fact = 3, nmaps = 4, ...)
{
  # Makes a series of maps showing spatio-temporal data

  # Arguments:
  # newdata   data.frame or data.table in a long format, containing square id-s (id), 
  #           coordinates in EPSG:2180 projection (x, y), time (year),
  #           and a variable "vname"
  # vname     variable to be shown on maps
  # bias      "bias" parameter passed to the function RColorBrewer::colorRampPalette
  # probs     probability for the function "quantile" to set the upper range for the legend
  # fact      "fact" parameter passed to the function terra::aggregate
  # nmaps     no. of maps to be drawn
  # ...       additional arguments passed to terra::plot
  
  require(RColorBrewer)
  require(data.table)
  require(terra)
  
  pal <- rev(brewer.pal(11, "RdYlGn"))
  colors <- colorRampPalette(pal, bias = bias)(256)
  
  if(is.data.table(newdata)) DT <- newdata[, c("id", "x", "y", "year", ..vname)] else
    DT <- data.table(newdata[c("id", "x", "y", "year", vname)])
  years <- sort(unique(DT[["year"]]))
  
  idxy <- DT[year == years[1], ]
  idxy <- idxy[, c("year", vname) := NULL]
  idxy[, `:=`(x = x * 1e3, y = y * 1e3)]
  
  dat <- DT[, c("x", "y") := NULL]
  w <- dcast(dat, id ~ year, value.var = vname)
  w <- merge(idxy, w, by = "id")
  w <- data.frame(w)
  r <- rast(w[-1], type = "xyz", crs = "epsg:2180")
  names(r) <- years
  if(fact == 1) r1 <- r else
    r1 <- aggregate(r, fact = fact, fun = "mean", na.rm = TRUE)
  rv <- as.vector(r1)
  
  if (is.null(range))
  {
    if (probs < 1) zakres <- c(min(rv, na.rm = TRUE), quantile(rv, na.rm = TRUE, probs)) else
      zakres <- range(rv, na.rm = TRUE)
  } else
    zakres <- range
  
  y <- round(seq(1, nlyr(r1), len = min(nmaps, nlyr(r1))))
  # plot(r1, y, col = colors, pax = list(labels = FALSE), range = zakres, fun = function() lines(pl, col = "grey40"), ...)
  plot(r1, y, col = colors, pax = list(labels = FALSE), range = zakres, ...)
  return(invisible(r1))
}



pp_check <- function(observed, simulated)
  # Posterior predictive check
  
{
  bpval <- mean(observed > simulated)
  ra <- range(observed, simulated)
  op <- par(mar = c(4.5, 4.5, 5.5, 2))
  plot(observed, simulated, asp = 1, xlim = ra, ylim = ra, xlab = "Actual Dataset", ylab = "Simulated Dataset", main = paste("Posterior Predictive Check", "\n", "Bayesian p-value =", round(bpval, 2)), col = adjustcolor("black", alpha = 1/4))
  abline(0, 1)
  grid()
  par(op)
  return(invisible(bpval))
}
