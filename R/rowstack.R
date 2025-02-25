rowstack <-
function(A,b)
# block stack b copies of matrix A rowwise
{
  B <- A[rep(seq(nrow(A)),b),]
  return(as.matrix(B))
}
