duplicated_first <-
function(x){
# duplicated_first() - sets T for the first of each set of duplicate entries
  y <- rep(F,length(x))
  y[match(x[duplicated(x)],x)] <- T
  return(y)
}
