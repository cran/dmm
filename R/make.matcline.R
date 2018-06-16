make.matcline <-
function(df)
# make.matcline  -  make a cytoplasmic maternal line code
{
  crel <- am.mcrel(df)
  m <- nrow(crel)
  lineno <- rep(0,m)
  nextlineno <- 0
  for(i in 1:m) {
    nextlineno <- nextlineno + 1
    if(lineno[i] == 0) {
      lineno[i] <- nextlineno
    }

    for(j in i:m) {
      if(crel[i,j] == 1) {
        if(lineno[j] == 0){
          lineno[j] <- nextlineno
        }
      }
    }
  }
  return(lineno)
}
