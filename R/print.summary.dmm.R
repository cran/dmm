print.summary.dmm <-
function(x, ...)
# print.summary.dmm() - format a summary.dmm object for printing
{
  cat("Call:\n")
  print(x$call)

  cat("\nCoefficients fitted by OLS for fixed effects:\n\n")
  for(i in 1:length(x$btables)) {
    print(x$btables[[i]],digits=x$digits)
    cat("\n")
  }
  cat("\nComponents partitioned by DME from residual var/covariance after OLS-fixed-effects fit:\n\n")
  for(i in 1:length(x$ctables)) {
    print(x$ctables[[i]],digits=x$digits)
    cat("\n")
  }

  if(x$fixedgls) {
    cat("\nCoefficients fitted by GLS for fixed effects:\n\n")
    for(i in 1: length(x$gbtables)){
      print(x$gbtables[[i]],digits=x$digits)
      cat("\n")
    }
    cat("\nComponents partitioned by DME from residual var/covariance after GLS-fixed-effects fit:\n\n")
    for(i in 1:length(x$gctables)) {
      print(x$gctables[[i]],digits=x$digits)
      cat("\n")
    }
  }
}
