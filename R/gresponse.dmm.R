gresponse.dmm <-
function(dmmobj,traitset="all",gls=F,psd=rep(1,length(traits)),response="Ia",digits=3, ...)
# gresponse.dmm() - predict genetic change given phenotypic selection differential
{
     if(traitset[1] == "all"){
    traits <- dimnames(dmmobj$b)[[2]][1:ncol(dmmobj$b)]
  }
  else {
    traits <- traitset
  }
  traitpairs <- permpaste(traits)
  l <- length(traits)
  alltraitpairs <- permpaste(dimnames(dmmobj$b)[[2]])

# individual additive responses
  if(response == "Ia") {
    if(!gls) {
      p <- dmmobj$phenotypic.variance[traits,traits]
      g <- matrix(dmmobj$variance.components["VarG(Ia)",traitpairs],l,l)
    }
    else if(gls) {
      p <- dmmobj$gls$phenotypic.variance[traits,traits]
      g <- matrix(dmmobj$gls$variance.components["VarG(Ia)",traitpairs],l,l)
    }
    
    dimnames(p) <- dimnames(dmmobj$phenotypic.variance[traits,traits])
    dimnames(g) <- dimnames(p)
    dimnames(g)[[1]] <- fixpaste(dimnames(p)[[1]],"Ia")
    dimnames(g)[[2]] <- dimnames(g)[[1]]
    psd <- as.matrix(psd,l,1)

    dimnames(psd)[[1]] <- dimnames(p)[[1]]
    gsd <- g %*% solve(p) %*% psd  # genetic selection differential
    ugsd <- g %*% solve(p)  # gsd for unit psd on each trait
    dsg <- solve(p) %*% psd  # directional selection gradient
    retobj <- list(gsd=gsd,ugsd=ugsd,dsg=dsg,gcov=g,pcov=p,psd=psd,digits=digits,
                    response=response)
  }


# individual additive and maternal additive responses
  if(response == "Ia+Ma"){
    if(!gls) {
      p <- dmmobj$phenotypic.variance[traits,traits]
      gii <- matrix(dmmobj$variance.components["VarG(Ia)",traitpairs],l,l)
      gim <- matrix(dmmobj$variance.components["CovG(Ia,Ma)",traitpairs],l,l)
      gmi <- matrix(dmmobj$variance.components["CovG(Ma,Ia)",traitpairs],l,l)
      gmm <- matrix(dmmobj$variance.components["VarG(Ma)",traitpairs],l,l)
    }
    else if(gls) {
      p <- dmmobj$gls$phenotypic.variance[traits,traits]
      gii <- matrix(dmmobj$gls$variance.components["VarG(Ia)",traitpairs],l,l)
      gim <- matrix(dmmobj$gls$variance.components["CovG(Ia,Ma)",traitpairs],l,l)
      gmi <- matrix(dmmobj$gls$variance.components["CovG(Ma,Ia)",traitpairs],l,l)
      gmm <- matrix(dmmobj$gls$variance.components["VarG(Ma)",traitpairs],l,l)
    }

    dimnames(p) <- dimnames(dmmobj$phenotypic.variance[traits,traits])
    dimnames(gii) <- dimnames(p)
    dimnames(gim) <- dimnames(p)
    dimnames(gmi) <- dimnames(p)
    dimnames(gmm) <- dimnames(p)
    dimnames(gii)[[1]] <- fixpaste(dimnames(p)[[1]],"Ia")
    dimnames(gii)[[2]] <- dimnames(gii)[[1]]
    dimnames(gmm)[[1]] <- fixpaste(dimnames(p)[[1]],"Ma")
    dimnames(gmm)[[2]] <- dimnames(gmm)[[1]]
    dimnames(gim)[[1]] <- dimnames(gii)[[1]]
    dimnames(gim)[[2]] <- dimnames(gmm)[[2]]
    dimnames(gmi)[[1]] <- dimnames(gmm)[[1]]
    dimnames(gmi)[[2]] <- dimnames(gii)[[2]]
    psd <- as.matrix(psd,l,1)
    dimnames(psd)[[1]] <- dimnames(p)[[1]]

    gpart <- cbind(rbind(gii,gmi),rbind(gim,gmm))
    r <- rbind(diag(l), 0.5 * diag(l))
    gsd <- gpart %*% r %*% solve(p) %*% psd  # genetic selection differential
    ugsd <- gpart %*% r %*% solve(p)  # gsd for unit psd on each trait
    dsg <- solve(p) %*% psd  # directional selection gradient
    retobj <- list(gsd=gsd,ugsd=ugsd,dsg=dsg,gcov=gpart,pcov=p,psd=psd,digits=digits,
                   response=response)
  }

  retobj$call <- match.call()
  class(retobj) <- "gresponse.dmm"
  return(retobj)
}
