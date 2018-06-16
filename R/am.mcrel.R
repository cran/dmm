am.mcrel <-
function(df){
  # am.mcrel()
  #compute maternal cytoplasmic relationship matrix from a data frame
  m <- length(df$Id)   # no of animals in data frame
  arel <- diag(m)
  for(i in 1:m) {
    mat <- df$DId[i]
    # off diagonal elements
    if (i > 1) {
      if(!is.na(mat)){
        arel[i,mat] <- 1
        arel[mat,i] <- 1
        for(j in 1:(i-1)) {
          arel[j,i] <- arel[j,mat]
          arel[i,j] <- arel[j,i]
        }
      }
    }
  }
  return(arel)
}
