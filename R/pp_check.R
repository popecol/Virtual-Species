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
