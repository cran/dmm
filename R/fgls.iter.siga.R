fgls.iter.siga <-
function(am, start.siga,mmat,dyad.explist,evec,dmefglsopt,dmeopt,ctable){
# fgls.iter.siga() - iterate and compute new siga from DME by fgls
# Multivariate version
    stopcrit <- 10
#   stoptol <- glsopt$stoptol *0.1
    stoptol <- dmefglsopt$stoptol 
    maxiter <- dmefglsopt$maxiter
    maxiter <- 100
    bdamp <- dmefglsopt$bdamp 
    newsiga <- start.siga
    count <- 0
    degf <- am$n - am$k
    bias <- am$n/degf
    ipkminv <- ginvipk(am$n,am$n)
    ipkminvblock <- kronecker(diag(am$l * am$l),ipkminv)  # blocked for multiv
#   ipkminvblock <- kronecker(diag(am$l),ipkminv)  # blocked for multiv
    xblock <- kronecker(diag(am$l), am$x) # xblock is lxl blocks each nxk
    hblock <- xblock %*% ginv(t(xblock) %*% xblock) %*% t(xblock)
    mmatblock <- diag(am$n * am$l ) - hblock # mmat has to be n*l because of hmat and V

#   hmat <- am$x %*% ginv(t(am$x) %*% am$x) %*% t(am$x)
#   hmat <- kronecker (diag(am$l), hmat)
#   mmatblock <- diag(am$n * am$l ) - kronecker (diag(am$l), hmat)

#   cat("mmatblock:\n")
#   print(mmatblock)
    evecblock <- matrix(evec,am$n * am$n *am$l * am$l,1)  # evecblock is vector n*n*l*l

    ematblock <- kronecker(diag(am$l * am$l),dyad.explist$emat)
#cat("evec:\n")
#print(evec)
#cat("evecblock:\n")
#print(evecblock)

    while (stopcrit > stoptol && count < maxiter) {
        cat("Iteration(fgls) round: ", count, "\n")
        oldsiga <- newsiga
#    compute V matrix (nl x nl)
        v <- expect.v(am,oldsiga,dyad.explist) # v is lxl blocks each nxn
#       v <- matrix(nearPD(v,ensureSymmetry=T)$mat,am$n * am$l,am$n * am$l)
#       cat("V matrix:\n")
#       print(v)
#       vinv <- ginv(v)
#    solve for siga
        newsiga.list <- fgls.update(am, v, mmatblock, ematblock, evecblock, ipkminvblock)
        newsiga <- newsiga.list$siga
#  cat("newsiga:\n")
#  print(newsiga)
#       dimnames(newsiga) <- dimnames(oldsiga)
#       newsiga <- siga.posdef(newsiga, am, ctable) # test posdef pre update
#    damp the siga update
        newsiga <- oldsiga + (newsiga - oldsiga) * bdamp
#       cat("Updated siga:\n")
#       print(newsiga)

#  check updated siga posdef
#       newsiga <- siga.posdef(newsiga, am, ctable)
#       cat("Updated siga made positive definite:\n")
#       print(newsiga)

#    look at stopcrit
        sumdev <- 0
        for (ll in 1:(am$l*am$l)) {
            for (i in 1:am$v) {
                sumdev <- sumdev + abs(newsiga[i, ll] - oldsiga[i, ll])
                # (newsiga - oldsiga) is always (new-old)*bdamp - ie part of the difference
            }
        }
        stopcrit <- sumdev/(am$l * am$l * am$v)
        count <- count + 1
        cat("Round = ",count," Stopcrit = ",stopcrit,"\n")
    }
#     end of iteration
      cat("Iteration(fgls) completed - count = ",count,"\n")
#   convergence check
      if(count == maxiter){
        cat("Failed to converge (fgls)\n")
        return()
      }
      cat("Convergence achieved (fgls)\n")

#   SE of siga 
    degfd <- am$n * am$n - am$v
#     vsiga <- kronecker(vard, solve(crossprod(qr.R(newemat.qr))), make.dimnames=T)
      vsiga <- newsiga.list$vsiga  # may be univariate ? OK now
      diagvsiga <- matrix(diag(vsiga),am$v,am$l * am$l)
      dimnames(diagvsiga) <- dimnames(start.siga)
cat("diagvsiga:\n")
print(diagvsiga)
      diagvsiga <- siga.posdef(diagvsiga, am, ctable)
cat("diagvsiga posdef:\n")
print(diagvsiga)
#cat("vsiga:\n")
#print(vsiga)
#     sesiga <- matrix(sqrt(diag(vsiga)), am$v, am$l * am$l, dimnames=dimnames(start.siga))
      sesiga <- sqrt(diagvsiga)
cat("sesiga:\n")
print(sesiga)
#    dimnames
     dimnames(newsiga) <- dimnames(start.siga)
     dimnames(sesiga) <- dimnames(start.siga)

#  determine "fit" for fgls
     fitted.values <- dyad.explist$emat %*% newsiga
     resid.values <- evec - fitted.values
     dme.fit.fgls <- list(y=evec, coefficients=newsiga, residuals=resid.values , fitted.values=fitted.values )

    parlist <- c(list(siga=newsiga,sesiga=sesiga, vsiga=vsiga), dme.fit.fgls) #  vard.ols=vard, msr.ols=msr, msrdf=degf, vara.ols=msa)
    return(parlist)
}
