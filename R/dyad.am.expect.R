dyad.am.expect <-
function(am, fixedgls,dmeopt,mmat){
# dyad.am.expect()
# evaluate parts of dyadic model equation  including emat and emat.qr
# am is an ante-model object
# am$zi, am$zm,am$zc are a lists of Z matrices
# am$rel is a list of  relationship matrices
# returns a list object containing emat, emat.qr, and cnames(col names for emat)
#
# Note: This code sets the order of components in siga[,]
#
# setup m matrix
# mmat <- diag(am$n) - am$x %*% ginv(am$x)
# cat("mmat:\n")
# print(mmat)
#
  dae <- list(cnames="VarE(I)",cnamesie=vector(mode="character",length=0),emat=matrix(0,am$n * am$n, am$v),vmat=matrix(0,am$n * am$n, am$v),icol=1,iecol=1)
# initial values
# cnames <- "VarE(I)"  # character vector
# emat <- matrix(0,am$n * am$n, am$v)
# vmat <- matrix(0,am$n * am$n, am$v)
# icol <- 1
  nsf <- length(am$specific.components)

  ctable <- make.ctable()

# nonspecific
  if(length(am$components) > 0 ) {
    for (kc in 1: length(am$components)) {
      index <- match(am$components[kc],ctable$all)
      pre <- ctable$allzpre[index]  
      post <- ctable$allzpost[index]  
      prel <- ctable$allrel[index]
      if((pre != "S") & (post != "S") & (prel != "I") ) { # normal cases
        zpre <- eval(parse(text=paste("am$",pre,sep="")))
        zpost <- eval(parse(text=paste("am$",post,sep="")))
        rel <- eval(parse(text=paste("am$rel$",prel,sep="")))
        dae <- dae.nonspecific(zpre,rel,zpost,mmat,am$components[kc],dae$cnames,dae$cnamesie,dae$emat,dae$vmat,dae$icol,dae$iecol,fixedgls,dmeopt)
      }
      else if((pre != "S") & (post != "S") & (prel == "I") ){ # cases with no rel matrix
        zpre <- eval(parse(text=paste("am$",pre,sep="")))
        zpost <- eval(parse(text=paste("am$",post,sep="")))
        dae <- dae.nonspecific.I(zpre,zpost,mmat,am$components[kc],dae$cnames,dae$cnamesie,dae$emat,dae$vmat,dae$icol,dae$iecol,fixedgls,dmeopt)
      }
      else if((pre == "S") & (post == "S") & (prel == "S") ) { # cases with 2 parts AND'ed
        indexs <- match(am$components[kc],ctable$cohort)
        pre1 <- ctable$cohortzpre1[indexs]
        pre2 <- ctable$cohortzpre2[indexs]
        post1 <- ctable$cohortzpost1[indexs]
        post2 <- ctable$cohortzpost2[indexs]
        op <- ctable$cohortop[indexs]
        zpre1 <- eval(parse(text=paste("am$",pre1,sep="")))
        zpost1 <- eval(parse(text=paste("am$",post1,sep="")))
        zpre2 <- eval(parse(text=paste("am$",pre2,sep="")))
        zpost2 <- eval(parse(text=paste("am$",post2,sep="")))
        zop <- paste(" ",op," ",sep="")
        dae <- dae.nonspecific.S(zpre1,zpost1,zpre2,zpost2,zop,mmat,am$components[kc],dae$cnames,dae$cnamesie,dae$emat,dae$vmat,dae$icol,dae$iecol,fixedgls,dmeopt)
      }
      else {   
        stop("Expectation for component not recognised:\n")
      }
    }
  }


  
# specific

  if(nsf > 0) {
    for(kf in 1:length(am$effnames)) {

      for(lc in 1:length(am$specific.components[[kf]])){
        index <- match(am$specific.components[[kf]][lc],ctable$all)
        pre <-  ctable$allzpre[index]
        post <- ctable$allzpost[index]
        prel <- ctable$allrel[index]

        if((pre != "S") & (post != "S") & (prel != "I") ) {  # normal cases
          zpre <- eval(parse(text=paste("am$",pre,sep="")))
          zpost <- eval(parse(text=paste("am$",post,sep="")))
          rel <- eval(parse(text=paste("am$rel$",prel,sep="")))
          dae <- dae.specific(zpre, rel, zpost, mmat, kf ,am$specific.components[[kf]][lc], am$effnames,am$effcodes,am$effnandc,am$comcodes,am$varcodes,dae$cnames,dae$cnamesie,dae$emat,dae$vmat,dae$icol,dae$iecol,fixedgls,dmeopt,ctable)
        }
        else if((pre != "S") & (post != "S") & (prel == "I") ) {  # cases with no rel matrix
          zpre <- eval(parse(text=paste("am$",pre,sep="")))
          zpost <- eval(parse(text=paste("am$",post,sep="")))
          dae <- dae.specific.I(zpre, zpost, mmat, kf ,am$specific.components[[kf]][lc], am$effnames,am$effcodes,am$effnandc,am$comcodes,am$varcodes,dae$cnames,dae$cnamesie,dae$emat,dae$vmat,dae$icol,dae$iecol,fixedgls,dmeopt,ctable)
        }
        else if((pre == "S") & (post == "S") & (prel == "S") ) {  # cases with 2 parts ANDed
      indexs <- match(am$specific.components[[kf]][lc],ctable$cohort)
      pre1 <- ctable$cohortzpre1[indexs]
      pre2 <- ctable$cohortzpre2[indexs]
      post1 <- ctable$cohortzpost1[indexs]
      post2 <- ctable$cohortzpost2[indexs]
      op <- ctable$cohortop[indexs]
      zpre1 <- eval(parse(text=paste("am$",pre1,sep="")))
      zpost1 <- eval(parse(text=paste("am$",post1,sep="")))
      zpre2 <- eval(parse(text=paste("am$",pre2,sep="")))
      zpost2 <- eval(parse(text=paste("am$",post2,sep="")))
      zop <- paste(" ",op," ",sep="")
      dae <- dae.specific.S(zpre1,zpost1,zpre2,zpost2,zop,mmat,kf,am$specific.components[[kf]][lc],am$effnames,am$effcodes,am$effnandc,am$comcodes,am$varcodes,dae$cnames,dae$cnamesie,dae$emat,dae$vmat,dae$icol,dae$iecol,fixedgls,dmeopt,ctable)
        }
        else {
          stop("Expectation for component not recognised:\n")
        }
      }
    }
  }
#
# resize emat cols <- dae$icol which is less than am$v if VarE(I) specific
# reset am$v
  cat("No of components defined = ",am$v,"\n")
  am$v <- dae$icol - 1
  cat("No of components estimable = ",am$v,"\n")
  emat <- matrix(dae$emat[,1:am$v],am$n * am$n, am$v)
  dimnames(emat) <- list(NULL,dae$cnames)

  vmat <- NULL
  if(fixedgls | dmeopt == "fgls") {
    vmat <- matrix(dae$vmat[,1:am$v],am$n * am$n, am$v)
    dimnames(vmat) <- list(NULL,dae$cnames)
  }

# summarize emat
    cat("Checking dyadic model equations:\n")
#   cat("emat:\n")
#   print(dae$emat)
#   emat.sum <- apply(emat,2,sum)
#   cat("column.sums:\n")
#   print(emat.sum)
    emat.mean <- apply(emat,2,mean)
#   cat("column.means:\n")
#   print(emat.mean)
    emat.var <- apply(emat,2,var)
#   cat("column.variances:\n")
#   print(emat.var)
    emat.cor <- cor(emat,dae$emat)
    colnames(emat.cor) <- rownames(emat.cor)
#   cat("column.correlations:\n")
#   print(emat.cor)

#  do QR transform on emat
  emat.qr <- qr(emat)

# check emat for E'E positive definite - ie emat.qr$rank should be am$v
  fullrank <- T
  if(emat.qr$rank != am$v) {
    cat(" Rank of DME matrix = ",emat.qr$rank," no of components = ",am$v,"\n")
    if(dmeopt != "pcr") {
      cat("Dyadic model equations not of full rank: either omit some components or try dmeopt='pcr' \n")
      fullrank <- F
      cat("Check outputobject$dme.correl to see which components are confounded:\n")
    }
  }

# do qr on vmat  - only need if fixedgls=T
  if(fixedgls | dmeopt == "fgls") {
    vmat.qr <- qr(dae$vmat)
  }
  else {
    vmat.qr <- NULL
  }

  explist <- list(emat=emat, emat.qr=emat.qr,
       vmat=vmat, vmat.qr=vmat.qr,emat.mean=emat.mean,
       emat.var=emat.var,emat.cor=emat.cor, newv=am$v,
       cnames=dae$cnames,cnamesie=dae$cnamesie,
       fullrank=fullrank) 
  return(explist)
}
