colstack <-
function(A,b)
# block stack b copies of matrix A columnwise
{
  B <- A[,rep(seq(ncol(A)),b)]
  return(as.matrix(B))
}
