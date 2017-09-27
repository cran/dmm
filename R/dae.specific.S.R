dae.specific.S <-
function(zpre1,zpost1,zpre2,zpost2,zop,mmat,factorno,componentname,effnames,effcodes,effnandc,comcodes,varcodes,cnames,cnamesie,emat,vmat,icol,iecol,gls,ctable)
# dae.specific.S() - expectations for a specific component
#               zpre, zpost are zi or zm or zc
#               returns incremented icol which is next vacant col of W
#               special case with two expectations ANDed together
#               no rel matrix for either expectation
{

#   zaz <- matrix(0,nrows(emat),nrows(emat)

  nfactorcodes <- length(effcodes[[factorno]])
  k <- 0  # counts estimable components
  kie <- 0  # counts inestimable components - pos in comcodes[[fno]][]
  if(any(!is.na(match(componentname,ctable$allindenv)))) {  # cross-class-cov not estimable
    for ( i in 1:nfactorcodes) {
      for( j in 1:nfactorcodes) {
        if ( i == j) {   # within class estimable
          k <- k + 1
          kie <- kie + 1
          specificcomponentname <- paste(varcodes[[factorno]][k],componentname,sep=":")
          szpre1 <- zpre1[[effnandc[[factorno]][i]]]
          szpost1 <- zpost1[[effnandc[[factorno]][j]]]
          szpre2 <- zpre2[[effnandc[[factorno]][i]]]
          szpost2 <- zpost2[[effnandc[[factorno]][j]]]
          zaz <- eval(parse(text=paste("(szpre1 %*% t(szpost1))", zop, "(szpre2 %*% t(szpost2))"," + 0",sep="")))
          emat[,icol] <- as.vector(mmat %*% zaz %*% mmat) # one col of W matrix
          if(gls) {
            vmat[,icol] <- as.vector(zaz) # one col of V matrix
          }
          cnames[icol] <- specificcomponentname # name for this col
          icol <- icol + 1
       }
       else { # cross class inestimable - just setup names
          kie <- kie + 1
          specificcomponentname <- paste(comcodes[[factorno]][kie],componentname,sep=":")
          cnamesie[iecol] <- specificcomponentname
          iecol <- iecol + 1
        }  # end cross class inestimable
      }
    }
  } # end cross-class cov not estimable

  else {   # cross-class-cov estimable
    for ( i in 1:nfactorcodes) {
      for( j in 1:nfactorcodes) {
        k <- k + 1
        specificcomponentname <- paste(comcodes[[factorno]][k],componentname,sep=":") 
        szpre1 <- zpre1[[effnandc[[factorno]][i]]]
        szpost1 <- zpost1[[effnandc[[factorno]][j]]]
        szpre2 <- zpre2[[effnandc[[factorno]][i]]]
        szpost2 <- zpost2[[effnandc[[factorno]][j]]]
        zaz <- eval(parse(text=paste("(szpre1 %*% t(szpost1))", zop, "(szpre2 %*% t(szpost2))"," + 0",sep="")))
        emat[,icol] <- as.vector(mmat %*% zaz %*% mmat) # one col of W matrix
        if(gls) {
          vmat[,icol] <- as.vector(zaz) # one col of V matrix
        }
        cnames[icol] <- specificcomponentname # name for this col
        icol <- icol + 1
      }
    }
  }  # # end cross-class cov estimable

   daelist <- list(cnames=cnames,cnamesie=cnamesie,emat=emat,vmat=vmat,icol=icol,iecol=iecol)
   return(daelist)
}
