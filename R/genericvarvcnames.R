genericvarvcnames <-
function(x)
#  genericvarvcnames()  - find the variance component part of a specific vc name
#                x is a vector of specific varvc names
#        need to take before 2nd : and after last :
{
  gvarvcname <- rep(0,length(x))
  for( i in 1:length(x)) {
    colons <- which(strsplit(x[i],"")[[1]] == ":")
    if(length(colons) == 0) {
      gvcname[i] <- x[i]
    }
    else {
      lastcolon <- colons[length(colons)]
      bracket <- which(strsplit(x[i],"")[[1]] == "(")
      if ( lastcolon > bracket) {
        lastcolon <- colons[length(colons)-1]
      }
      secondcolon <- colons[2]
      gvarvcname[i] <- paste(substr(x[i],1,secondcolon),substr(x[i],lastcolon+1,nchar(x[i])),sep="")
    }
  }
  return(gvarvcname)
}
