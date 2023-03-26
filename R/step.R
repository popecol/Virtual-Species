
# Backward model selection for 'gam' objects.
# Author: lechu@amu.edu.pl


backward <- function(obj, p = 0.05)
{
  su <- summary(obj, re.test = FALSE)
  A <- to.remove(su)
  while (any(A$p.value > p))
  {
    idx <- which.max(A$p.value)
    pval <- prettyNum(A$p.value)[idx]
    dv <- rownames(A)[idx]
    cat("Removing", dv, "     p-value =", pval, "\n")
    dterm <- paste0(". ~ . - ", paste0(substr(dv, 1, nchar(dv) - 1), ", k = k, bs = 'cr')"))
    newf <- update(formula(obj), dterm)
    gc(verbose = FALSE)
    obj <- update(obj, newf)
    su <- summary(obj, re.test = FALSE)
    A <- to.remove(su)
  }
  cat("All p-values below", p, "\n")
  return(obj)
}


to.remove <- function(su) {
  # Helper function extracting p-values from a 'summary.gam' object, 
  # excluding random factors, time and space.
  
  A <- data.frame(su$s.table)
  A <- A[!grepl("(plot_id|observer_id|year|fyear|(x,y))", rownames(A)), ]
  A <- transform(A, edf = round(edf, 4), p.value = round(p.value, 4))
}

