selfpaste <-
function(x)
# selfpaste() - paste elements of vector x with themselves
{
  n <- 0
  selfx <- rep(0,length(x)^2)
  for(i in 1: length(x)){
    n <- n + 1
    selfx[n] <- paste(x[i],":",x[i],sep="", collapse=NULL)
  }
  return(selfx)
}
