pedcheck <-
function(df){
# pedcheck()  -  check Id, SId, DId in dataframe df are valid for dmm
# Id
  d <- diff(as.numeric(df$Id)) 
  if(length(unique(d)) != 1) {
    stop("Id's must be unique:\n")
  }
  if(as.numeric(df$Id[1]) != 1) {
    stop("Id's must start at 1:\n")
  }
  if(d[1] != 1) {
    stop("Id's must be an arithmetic sequence:\n")
  }
# SId
  if(any(is.na(match(df$SId[!is.na(df$SId)],df$Id)))){
    stop("SId's must occur as an Id in the dataframe:\n")
  }
# DId
  if(any(is.na(match(df$DId[!is.na(df$DId)],df$Id)))){
    stop("DId's must occur as an Id in the dataframe:\n")
  }
  return(0)
}
