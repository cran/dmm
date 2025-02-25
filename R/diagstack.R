diagstack <-
function(A,b)
# stack b copies of matrix A diagonally
{
  B <- kronecker(diag(b),A,make.dimnames=T)
  return(B)
}
