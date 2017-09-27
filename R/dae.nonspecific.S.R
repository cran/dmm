dae.nonspecific.S <-
function(zpre1,zpost1,zpre2,zpost2,zop,mmat,componentname,cnames,cnamesie,emat,vmat,icol,iecol,gls)
# dae.nonspecific.S() - expectations for a nonspecific component
#             zpre,zposr are zi or zm or zc
#             special case with two expectations anded together
#             no rel matrix for either expectation
{

  nszpre1 <- zpre1[["NS"]]
  nszpost1 <- zpost1[["NS"]]
  nszpre2 <- zpre2[["NS"]]
  nszpost2 <- zpost2[["NS"]]
#   zaz <- matrix(0,nrows(emat),nrows(emat)
    zaz <- eval(parse(text=paste("(nszpre1 %*% t(nszpost1))", zop, " (nszpre2 %*% t(nszpost2))"," + 0",sep="")))
    emat[,icol] <- as.vector(mmat %*% zaz %*% mmat) # one col of W matrix
    if(gls) {
      vmat[,icol] <- as.vector(zaz) # one col of V matrix
    }
    cnames[icol] <- componentname # name for this col
    icol <- icol + 1

   daelist <- list(cnames=cnames,cnamesie=cnamesie,emat=emat,vmat=vmat,icol=icol,iecol=iecol)
   return(daelist)
}
