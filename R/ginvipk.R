ginvipk <-
function(m,n)
# inverse of [I(n*n) + commutation matrix or order m x n]
{
 K <- matrix(0,m*n,m*n)
 for (i in 1 : m){
  for (j in 1 : n){
    K[i + m*(j - 1), j + n*(i - 1)] <- 0.25
  }
 }
 for (i in 1 : (m*n)){
    K[i,i] <- K[i,i] + 0.25
 }
return(K)
}
