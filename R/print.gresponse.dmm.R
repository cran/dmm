print.gresponse.dmm <-
function(x, ...)
# print.gresponse.dmm() - format a gresponse.dmm object for printing
{
  cat("Call:\n")
  print(x$call)
  cat("\nPredicted response to selection using component(s) ",x$response,"\n")

  cat("\nGenetic selection differentials achieved by given psd:\n\n")
  print(x$gsd,digits=x$digits)
  cat("\n")
 
  cat("\nGenetic selection differentials achieved by a unit psd on each trait:\n\n")
  print(x$ugsd,digits=x$digits)
  cat("\n")

  cat("Directional selection gradient:\n\n")
  print(x$dsg,digits=x$digits)
  cat("\n")
 
  cat("Genetic var/covariance matrix:\n\n")
  print(x$gcov,digits=x$digits)
  cat("\n")
 
  cat("Phenotypic var/covariance matrix:\n\n")
  print(x$pcov,digits=x$digits)
  cat("\n")
  
  cat("Phenotypic selection differentials:\n\n")
  print(x$psd,digits=x$digits)
  cat("\n")
}
