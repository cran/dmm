blockvec <-
function(B,blocksize)
# unblock a matrix B in order 1,1  1,2  2,1  2,2 ..
# and put blocks into a 'column vector of blocks' in above order
# an make result a matrix
# n must be an exact divisor of nrow (b) and ncol(b)
# B must be square
{
  nblocks <- (nrow(B)/blocksize) * (ncol(B)/blocksize)
  blockspercol <- ncol(B)/blocksize
  colsout <- ncol(B)/sqrt(nblocks)
  rowsout <- colsout * nblocks
  A <- matrix(0,rowsout,colsout)
  for (i in 1: blockspercol) {
    for (j in 1:blockspercol){
      for(k in 1:blocksize) {
        for(l in 1:blocksize){
#         ka <-  k + (i-1) * blocksize + (j-1) * blocksize
          ka <-  k + (i-1) * blocksize * blockspercol + (j-1) * blocksize
          kb <-  (i-1) * blocksize + k
          la <-  l
          lb <-  (j-1) * blocksize + l
      A[ka,la] <- B[kb,lb]
        }
      }
    }
  }
  return(A)
}
