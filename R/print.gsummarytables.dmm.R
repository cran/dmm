print.gsummarytables.dmm <-
function(x,...)
# print.gsummarytables.dmm()  -  format a gsummarytables.dmm object for printing
{
  cat("Proportion of phenotypic var/covariance to each component:\n")
  for(i in 1: length(x$ftables)) {
    print(x$ftables[[i]],digits=x$digits)
    cat("\n")
  }
  cat("Correlation corresponding to each var/covariance component:\n")
  for(i in 1: length(x$rtables)) {
    print(x$rtables[[i]],digits=x$digits)
    cat("\n")
  }
  cat("Phenotypic var/covariance from summing components:\n")
  for(i in 1: length(x$ptables)) {
    print(x$ptables[[i]],digits=x$digits)
    cat("\n")
  }
}
