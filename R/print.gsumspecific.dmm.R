print.gsumspecific.dmm <-
function(x, ...)
# print.gsumspecific.dmm() - format a gsumspecific.dmm object for printing
{
  cat("Call:\n")
  print(x$call)

    cat("\nComponents partitioned by DME from residual var/covariance after OLS-b fit:\n\n")
  for(ic in 1:length(x$csumlist)) {
    cat("\nSpecific class: ",names(x$csumlist[ic]),"\n")
    print(x$csumlist[[ic]],digits=x$digits)
    cat("\n")
  }

  if(x$gls) {
    cat("\nComponents partitioned by DME from residual var/covariance after GLS-b fit:\n\n")
    for(ic in 1:length(x$gcsumlist)) {
      cat("\nSpecific class: ",names(x$gcsumlist[ic]),"\n")
      print(x$gcsumlist[[ic]],digits=x$digits)
      cat("\n")
    }
  }
}
