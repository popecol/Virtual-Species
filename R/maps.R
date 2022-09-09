maps <- function(newdata, vname, range = NULL, bias = 2, probs = 0.9999, fact = 3, nmaps = 4, ...)
{
  # Makes a series of maps showing spatio-temporal data
  # Author: lechu@amu.edu.pl
  
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
