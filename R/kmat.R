kmat <-
function(m,n)
# commutation matrix or order m x n
{
 K <- matrix(0,m*n,m*n)
 for (i in 1 : m){
  for (j in 1 : n){
    K[i + m*(j - 1), j + n*(i - 1)] <- 1
  }
 }
return(K)
}
