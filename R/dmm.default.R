dmm.default <-
function(mdf,fixform = Ymat ~ 1,components=c("VarE(I)","VarG(Ia)"),specific.components=NULL,cohortform=NULL,posdef=T,gls=F,glsopt=list(maxiter=200,bdamp=0.8,stoptol=0.01), dmeopt="qr",ncomp.pcr="rank",relmat="inline",dmekeep=F,dmekeepfit=F,traitspairwise=F,traitsblockwise=F,...)
{
# dmm.default()
  if(traitspairwise && traitsblockwise) {
    stop("dmm: not a valid option combination\n")
  }
  if(!traitspairwise && !traitsblockwise) {
    cat("Dyadic mixed model fit for datafile:",substitute(mdf)," \n")
    outlist <- dmesolve(mdf,fixform,components,specific.components,cohortform,posdef,gls,glsopt,dmeopt,ncomp.pcr,relmat,dmekeep,dmekeepfit)
#   outlist$call <- match.call()
    attributes(outlist) <- c(attributes(outlist),list(call = match.call()))
    class(outlist) <- "dmm"
  }
  else if(traitspairwise) {
    cat("Traitspairwise mixed model fit for datafile:",substitute(mdf)," \n")
    outlist <- dmm_array(mdf,fixform,components,specific.components,cohortform,posdef,gls,glsopt,dmeopt,ncomp.pcr,relmat,dmekeep,dmekeepfit)
    attributes(outlist) <- c(attributes(outlist),list(call = match.call()))
#   class(outlist) <- "dmmarray"
  }
  else if(traitsblockwise) {
    cat("Traitsblockwise mixed model fit for datafile:",substitute(mdf)," \n")
    outlist <- dmm_blockarray(mdf,fixform,components,specific.components,cohortform,posdef,gls,glsopt,dmeopt,ncomp.pcr,relmat,dmekeep,dmekeepfit,...)
    attributes(outlist) <- c(attributes(outlist),list(call = match.call()))
#   class(outlist) <- "dmmblockarray"
  }
# outlist$call <- match.call()
  return(outlist)
}
