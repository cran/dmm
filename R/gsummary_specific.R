gsummary_specific <-
function(dmmobj,traitset,componentset,bytrait,gls,digits,...)
# gsummary_specific() - summarize component estimates for specific case
{
    if(traitset[1] == "all"){
      traitpairs <- dimnames(dmmobj$siga)[[2]]
      traits <- traitpairstotraits(traitpairs)
    }
    else {
      traits <- traitset
      traitpairs <- permpaste(traits)
    }
    l <- length(traits)
    alltraitpairs <- dimnames(dmmobj$siga)[[2]]
 
    if(componentset[1] == "all") {
      oneclass <- names(dmmobj$specific)[1]
      components <- dimnames(dmmobj$specific[[oneclass]]$variance.components)[[1]]
    }
    else {
      components <- componentset
    }
  nc <- length(names(dmmobj$specific))
  csumlist <- vector("list",nc)
  ic <- 1
  for(kc in names(dmmobj$specific)) {
    csumlist[[ic]] <- make.gsummarytables(dmmobj$specific[[kc]],traitset,componentset,bytrait,gls,digits, ...)
    ic <- ic + 1
  }
  names(csumlist) <- names(dmmobj$specific)
  retobj <- list(csumlist=csumlist,traits=traits,components=components,bytrait=bytrait,gls=gls,digits=digits)

  if(gls) {
    nc <- length(names(dmmobj$gls$specific))
    gcsumlist <- vector("list",nc)
    ic <- 1
    for(kc in names(dmmobj$gls$specific)) {
      gcsumlist[[ic]] <- make.gsummarytables(dmmobj$gls$specific[[kc]],traitset,componentset,bytrait,gls,digits, ...)
      ic <- ic + 1
    }
    names(gcsumlist) <- names(dmmobj$gls$specific)
    retobj <- list(csumlist=csumlist,gcsumlist=gcsumlist,traits=traits,components=components,bytrait=bytrait,gls=gls,digits=digits)
  }

  retobj$call <- match.call()
  class(retobj) <- "gsumspecific.dmm"
  return(retobj)
}
