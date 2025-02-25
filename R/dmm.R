dmm <-
function(mdf,fixform = Ymat ~ 1,components=c("VarE(I)","VarG(Ia)"),specific.components=NULL,cohortform=NULL,posdef=T,fixedgls=F,fixedglsopt=list(maxiter=200,bdamp=0.8,stoptol=0.01),dmefglsopt=list(maxiter=100,bdamp=0.8,stoptol=0.001), dmeopt="qr",ncomp.pcr="rank",relmat="inline",dmekeep=F,dmekeepfit=F,traitspairwise=F,traitsblockwise=F,...) 
{
# dmm()
   UseMethod("dmm")
}
