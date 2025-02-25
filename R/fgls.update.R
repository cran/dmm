fgls.update <-
function(am,v,mmat,emat,evec,ipkminv)
# fgls.update()
# fgls estimate of siga given oldsiga to estimate V
# multivariate version
{
#cat("mmat:\n")
#print(mmat)
#cat("V:\n")
#print(v)
  bmvm <- mmat %*% v %*% mmat
#cat("bmvm:\n")
#print(bmvm)
#   bkb <- kronecker(blockvec(bmvm,am$n),t(blockvec(bmvm,am$n))) # vec(B-1) %x%t(vec(B-1)) # version 3
#   bkbinv <- ginv(bkb)
#   cat("bkb:\n")
#   print(bkb)
#   cat("bkb inverted:\n")
#   print(bkbinv)
  bmvminv <- ginv(bmvm)
# bkinv <- kronecker(bmvminv,bmvminv)  # version 1
# bkinv <- kronecker(matrix(t(bmvminv),am$n * am$n * am$l * am$l,1),t(matrix(t(bmvminv),am$n * am$n * am$l * am$l,1))) # vec(B-1) %x% t(vec(B-1))  # version 2
  bkinv <- kronecker(blockvec(bmvminv,am$n),t(blockvec(bmvminv,am$n))) # vec(B-1) %x% t(vec(B-1)) # version 3
#    bkinv <- bkbinv  # kludge to use bkb in kronecker then invert
#cat("bmvminv:\n")
#print(bmvminv)
#cat("blockvec bmvminv:\n")
#print(blockvec(bmvminv,am$n))
#stop("stop:\n")
#cat("bkinv:\n")
#print(bkinv)
#cat("ipkminv:\n")
#print(nrow(ipkminv))
#print(ncol(ipkminv))
#print(ipkminv)
  omegainv <- ipkminv %*% bkinv
#cat("omegainv:\n")
#print(diag(omegainv))
# cat("emat:\n")
# print(emat)
  xtom <- t(emat) %*% omegainv
# cat("xtom:\n")
# print(xtom)
  xtomx <- xtom %*% emat
# cat("xtomx:\n")
# print(xtomx)
# xtomx <- matrix(nearPD(xtomx,ensureSymmetry=T)$mat,am$v * am$l * am$l, am$v * am$l * am$l)
  xtomxi <- ginv(xtomx)
  dimnames(xtomxi) <- dimnames(xtomx)
#cat("xtomxi:\n")
#print(xtomxi)
# cat("evec:\n")
# print(evec)
  xtomy <- xtom %*% evec
#cat("xtomy:\n")
#print(xtomy)
  blocksiga <- xtomxi %*% xtomy
# cat("blocksiga in fgls.update:\n")
# print(blocksiga)
  
  siga <- matrix(blocksiga,am$v,am$l*am$l)
 cat("siga in fgls.update:\n")
   print(siga)
# vsiga <- solve(t(emat) %*% emat) %*% siga
# vsiga <- diag(am$l*am$l) %x% xtomxi
  vsiga <- xtomxi
#cat("vsiga in fgls.update:\n")
#print(vsiga)
  
# outlist <- list(siga=siga,vsiga=xtomxi)
  outlist <- list(siga=siga,vsiga=vsiga)
  return(outlist)
}
