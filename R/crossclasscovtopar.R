crossclasscovtopar <-
function(covx, covt, var1, var2){
# crossclasscovtopar() - covs to parameters for cross-class cases
#                       - same trait and cross trait cases
#                       - var1 and var2 are 'values' not names
#             - need to copy var1 and var2 out of siga before call this routine
  n <- ncol(covx)  # n traits
  if(ncol(covt) != n) {
    stop("Covx and Covt must be same size in crosseffectcovtopar():\n")
  }
  fract <- rep(0,n)
  corre <- matrix(0,n,n)
  for(j in 1:n) {
    jj <- (j-1)*n + j
    for(i in 1:n) {
      ii <- (i-1)*n + i
      if(i == j) {
       fract[i] <- covx[i,i]/covt[i,i]
      }
      able <- var1[jj] * var2[ii]
      if(!is.na(able)) {
        if(able > 0) {
          corre[i,j] <- covx[i,j]/sqrt(able)
        }
        else {
          corre[i,j] <- 0
        }
      }
      else {
        corre[i,j] <- NA
      }
    }
  }
  return(list(corre=corre, fract=fract))
}
