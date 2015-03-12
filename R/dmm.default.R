dmm.default <-
function(mdf,fixform = Ymat ~ 1,components=c("VarE(I)","VarG(Ia)"),cohortform=NULL,posdef=T,gls=F,glsopt=list(maxiter=200,bdamp=0.8,stoptol=0.01), dmeopt="qr",ncomp.pcr="rank",relmat="inline",dmekeep=F,dmekeepfit=F){
# dmm.default()
  cat("Dyadic mixed model fit for datafile:",substitute(mdf)," \n")
  outlist <- dmesolve(mdf,fixform,components,cohortform,posdef,gls,glsopt,dmeopt,ncomp.pcr,relmat,dmekeep,dmekeepfit)
  outlist$call <- match.call()
  class(outlist) <- "dmm"
  return(outlist)
}
