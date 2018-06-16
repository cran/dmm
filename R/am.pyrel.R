am.pyrel <-
function(df){
  # am.pyrel()
  #compute paternal Y-chromosome relationship matrix from a data frame
  m <- length(df$Id)   # no of animals in data frame
  yrel <- diag(m)
  for(i in 1:m) {
    pat <- df$SId[i]
    # off diagonal elements
    if (i > 1) {
      if(!is.na(pat)){
        yrel[i,pat] <- 1
        yrel[pat,i] <- 1
        for(j in 1:(i-1)) {
          yrel[j,i] <- yrel[j,pat]
          yrel[i,j] <- yrel[j,i]
        }
      }
    }
  }
  return(yrel)
}
