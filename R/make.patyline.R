make.patyline <-
function(df)
# make.patyline  -  make a paternal Y-chromosome line code
{
  yrel <- am.pyrel(df)
  m <- nrow(yrel)
  lineno <- rep(0,m)
  nextlineno <- 0
  for(i in 1:m) {
    nextlineno <- nextlineno + 1
    if(lineno[i] == 0) {
      lineno[i] <- nextlineno
    }

    for(j in i:m) {
      if(yrel[i,j] == 1) {
        if(lineno[j] == 0){
          lineno[j] <- nextlineno
        }
      }
    }
  }
  return(lineno)
}
