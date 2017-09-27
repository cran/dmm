print.csummary.dmm <-
function(x, ...)
# print.csummary.dmm() - format a csummary.dmm object for printing
{
  cat("Call:\n")
  print(x$call)

  cat("\nComponents partitioned by DME from residual var/covariance after OLS-b fit:\n\n")
  for(i in 1:length(x$ctables)) {
    print(x$ctables[[i]],digits=x$digits)
    cat("\n")
  }

  if(x$gls) {
    cat("\nComponents partitioned by DME from residual var/covariance after GLS-b fit:\n\n")
    for(i in 1:length(x$gctables)) {
      print(x$gctables[[i]],digits=x$digits)
      cat("\n")
    }
  }
}
