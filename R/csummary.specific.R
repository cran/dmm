csummary.specific <-
function(object,traitset,componentset,bytrait,gls,digits,...)
# csummary.specific() - summarize component estimates for specific case
{
    if(traitset[1] == "all"){
      traitpairs <- dimnames(object$siga)[[2]]
      traits <- traitpairstotraits(traitpairs)
    }
    else {
      traits <- traitset
      traitpairs <- permpaste(traits)
    }
    l <- length(traits)
    alltraitpairs <- dimnames(object$siga)[[2]]
 
    if(componentset[1] == "all") {
      oneclass <- names(object$specific)[1]
      components <- dimnames(object$specific[[oneclass]]$variance.components)[[1]]
    }
    else {
      components <- componentset
    }
  nc <- length(names(object$specific))
  csumlist <- vector("list",nc)
  ic <- 1
  for(kc in names(object$specific)) {
    csumlist[[ic]] <- make.csummarytables(object$specific[[kc]],traitset,componentset,bytrait,gls,digits, ...)
    ic <- ic + 1
  }
  names(csumlist) <- names(object$specific)
  retobj <- list(csumlist=csumlist,traits=traits,components=components,bytrait=bytrait,gls=gls,digits=digits)

  if(gls) {
    nc <- length(names(object$gls$specific))
    gcsumlist <- vector("list",nc)
    ic <- 1
    for(kc in names(object$gls$specific)) {
      gcsumlist[[ic]] <- make.csummarytables(object$gls$specific[[kc]],traitset,componentset,bytrait,gls,digits, ...)
      ic <- ic + 1
    }
    names(gcsumlist) <- names(object$gls$specific)
    retobj <- list(csumlist=csumlist,gcsumlist=gcsumlist,traits=traits,components=components,bytrait=bytrait,gls=gls,digits=digits)
  }

  retobj$call <- match.call()
  class(retobj) <- "csumspecific.dmm"
  return(retobj)
}
