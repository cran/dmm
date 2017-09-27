comtopar.specific <-
function(v, l, thisvc,vara, vthisvc, sethisvc, ctable,ic,phencovclass,longrownames,siga){
# comtopar.specific() - special comtopar for when some components are specific
#            - calls covtopar for each animal var/cov component
#            - to get hsq (fraction) and correlation for component
#	     - if component is a covariance, correlation uses appropriate variance
#              components in denominator - these variance components must be fitted
#            - sigt is total animal var/cov - ie phenotypic var/cov matrix - l x l
#            - thisvc is all v components var/cov's - v x l^2
#            - also computes Rp - phenotypic correlation matrix
#            - also computes SE's
#            - v = no of components, l = no of traits
#            - in outlist - labelc is vector of component names - v+1
#                         - correa is component correlations - v+1 x l^2
#                         - fracta is component variance fractions - v+1 x l

# adjust vthisvc to ensure posdef
# vthisvc <- nearPD(vthisvc,ensureSymmetry=T)$mat

# calculate sigt - sum of components
  sigt <- matrix(0,l,l, dimnames=dimnames(vara))
  sigt <- matrix(apply(thisvc,2,sum),l,l,dimnames=dimnames(vara))
# sampling variances and SE's for sigt
  vsigt <- matrix(0,l,l,dimnames=dimnames(vara))
  for (il in 1:l) {
    ib <- (il-1)*l+il
    for(jl in 1:l) {
      jb <- (jl-1)*l+jl
#     vsigt[ib,jb] <- vt(v,l,vthisvc,il,jl)
      vsigt[il,jl] <- vt(v,l,vthisvc,il,jl)
    }
  }
  sesigt <- sqrt(vsigt)

# components
  correa <- matrix(0,v+1,l^2, dimnames=list(c(rownames(thisvc),"VarP(I)"), colnames(thisvc)))
  fracta <- matrix(0,v+1,l,dimnames=list(c(rownames(thisvc),"VarP(I)"), colnames(sigt)))
  varcomp <- matrix(0,v+1,l^2,dimnames=list(c(rownames(thisvc),"VarP(I)"), colnames(thisvc)))
  sevarcomp <- matrix(0,v+1,l^2,dimnames=list(c(rownames(thisvc),"VarP(I)"), colnames(thisvc)))
  labelc <- rep(NULL,v+1)

  for(i in 1:v) {
    covi <- matrix(thisvc[i, ], l, l)
    varcomp[i, ] <- thisvc[i, ]
    sevarcomp[i, ] <- sethisvc[i, ]

# case 1
     if(!is.cecov(rownames(thisvc)[i],ctable$allcov)) {  # not a crosseffect cov
      # noncecov component - use covtopar() if w/n class specific
      #                    - use crossclasscovtopar() if cross class specific
      if(!is.specific(longrownames[i]) | is.withinclass(longrownames[i])) {  # within class or nonspecific
        pari <- covtopar(covi, sigt)

      }  # end w/n class var

# case 2
      else {  # cross class covariance
        # need to find corresponding variances
        varlist <- match.vars(longrownames[i])
        var1 <- siga[varlist$var1name,]
        var2 <- siga[varlist$var2name,]
#       var1 <- vc[[varlist$var1class]][varlist$genericname,]
#       var2 <- vc[[varlist$var2class]][varlist$genericname,]
        pari <- crossclasscovtopar(covi,sigt,var1,var2)

      }  # end cross class var
    } # end var component

# case 3
    else {
      # covariance component - use crosseffectcovtopar()
      if(!is.specific(longrownames[i]) | is.withinclass(longrownames[i])) { # within-class variance or nonspecific

        varlist <- match.crosseffect.vars(longrownames[i])
        var1 <- siga[varlist$var1name, ]
        var2 <- siga[varlist$var2name, ]
#       var1 <- vc[[varlist$var1class]][varlist$genericname1, ]
#       var2 <- vc[[varlist$var2class]][varlist$genericname2, ]
        pari <- crossclasscovtopar(covi,sigt,var1,var2)

      } # end within class 

# case 4
      else { # cross class cross effect cov
        # note if a cross effect cov is cross class, both the correaponding variances must be cross class for the same specific effect
        varlist <- match.crosseffect.vars(longrownames[i])
        var1 <- siga[varlist$var1name, ]
        var2 <- siga[varlist$var2name, ]
#       var1 <- vc[[varlist$var1class]][varlist$genericname1, ]
#       var2 <- vc[[varlist$var2class]][varlist$genericname2, ]
        pari <- crossclasscovtopar(covi,sigt,var1,var2)
      }  # end cross class cross effect cov

    }  # end else cov component

#   save stuff from pari for ith component
    correa [i, ] <- pari$corre
    fracta[i, ] <- pari$fract
    labelc[i] <- dimnames(thisvc)[[1]][i]
  }  # end for i
  
# phenotypic variance component
  pari <- covtopar(sigt, sigt)
  correa[v+1, ] <- pari$corre
  fracta[v+1, ] <- pari$fract
  varcomp[v+1, ] <- matrix(sigt,1,l^2)
  sevarcomp[v+1, ] <- as.vector(sesigt)
  labelc[v+1] <- "VarP(I)"

# SE's
  sep.list <- separ(varcomp, vthisvc, v,l, fracta, correa)

  outlist <- list(component=labelc,
          phencovclass=phencovclass, component.longnames=longrownames,
          correlation=correa, correlation.variance=sep.list$vcorre,
          correlation.se=sqrt(sep.list$vcorre),
          fraction=fracta, fraction.variance=sep.list$vfract,
          fraction.se=sqrt(sep.list$vfract),
          variance.components=varcomp,variance.components.se=sevarcomp,
          phenotypic.variance=sigt, phenotypic.variance.se=sesigt,
          observed.variance=vara)
  class(outlist) <- "dmm"
  return(outlist)
}
