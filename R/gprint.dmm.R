gprint.dmm <-
function(x, traitset="all",fixedgls=F, ...)
#  gprint.dmm() - short description of genetic parameters from a dmm fitted model object
{
  cat("Call:\n")
  print(x$call)
  cat("Fixed formula:\n")
  print(x$fixform)
  cat("Cohort formula:\n")
  print(x$cohortform)
  cat("Traits:\n")
  if(traitset[1] == "all"){
    traits <- dimnames(x$b)[[2]][1:ncol(x$b)]
  }
  else {
    traits <- traitset
  }
  print(traits)
  traitpairs <- permpaste(traits)

  if(!is.null(x$specific)) {  # specific
    cat("Specific classes:\n")
    print(names(x$specific))
    for(kc in names(x$specific)) {
      cat("\nSpecific class:",kc,"\n")
      cat("Var/covariance components as a proportion of OLS-fixed-effects phenotypic var/covariance:\n")
      print(x$specific[[kc]]$fraction[ ,traits])
      cat("Correlations corresponding to OLS-fixed-effects var/covariance components:\n")
      print(x$specific[[kc]]$correlation[ ,traitpairs])
      cat("Phenotypic var/covariance OLS-fixed-effects:\n")
      print(x$specific[[kc]]$phenotypic.variance[traits,traits])
    }
 
    if(fixedgls) {
      cat("Var/covariance components as a proportion of GLS-fixed-effects phenotypic var/covariance:\n")
      print(x$gls$specific[[kc]]$fraction[ ,traits])
      cat("Correlations corresponding to GLS-fixed-effects var/covariance components:\n")
      print(x$gls$specific[[kc]]$correlation[ ,traitpairs])
      cat("Phenotypic var/covariance GLS-fixed-effects:\n")
      print(x$gls$specific[[kc]]$phenotypic.variance[traits,traits])
    }
  }

  else {  # nonspecific
    cat("Var/covariance components as a proportion of OLS-fixed-effects phenotypic var/covariance:\n")
    print(x$fraction[ ,traits])
    cat("Correlations corresponding to OLS-fixed-effects var/covariance components:\n")
    print(x$correlation[ ,traitpairs])
    cat("Phenotypic var/covariance OLS-fixed-effects:\n")
    print(x$phenotypic.variance[traits,traits])
 
    if(fixedgls) {
      cat("Var/covariance components as a proportion of GLS-fixed-effects phenotypic var/covariance:\n")
      print(x$gls$fraction[ ,traits])
      cat("Correlations corresponding to GLS-fixed-effects var/covariance components:\n")
      print(x$gls$correlation[ ,traitpairs])
      cat("Phenotypic var/covariance GLS-fixed-effects:\n")
      print(x$gls$phenotypic.variance[traits,traits])
    }
  }
}
