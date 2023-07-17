csummary.dmm <-
function(object,traitset="all",componentset="all",bytrait=T,gls=F,digits=3, ...)
# csummary.dmm()  - make component summary tables for a dmm object
{

  if(!is.null(object$specific)) {  # specific
    retobj <- csummary_specific(object,traitset,componentset,bytrait,gls,digits,...)
    class(retobj) <- "csumspecific.dmm"
    return(retobj)
  } #  end if specific

  else {  # nonspecific
    if(traitset[1] == "all"){
#     traits <- dimnames(object$b)[[2]][1:ncol(object$b)]
      traitpairs <- dimnames(object$variance.components)[[2]]
      traits <- traitpairstotraits(traitpairs)
    }
    else {
      traits <- traitset
      traitpairs <- permpaste(traits)
    }
    l <- length(traits)
    alltraitpairs <- dimnames(object$variance.components)[[2]]
 
    if(componentset[1] == "all") {
      components <- dimnames(object$variance.components)[[1]]
    }
    else {
      components <- componentset
    }

  if(bytrait) {

    ctables <- vector("list",l*l)  # one table per traitpair
    count <- 0
    for(i in traits) {
      for(j in traits) {
        traitpair <- paste(i,":",j,sep="",collapse=NULL)
        ij <- match(traitpair,alltraitpairs)
        count <- count + 1
        ci95lo <- object$variance.components[components,ij] - 1.96 * object$variance.components.se[components,ij]
        ci95hi <- object$variance.components[components,ij] + 1.96 * object$variance.components.se[components,ij]
        ctable <- data.frame(Traitpair=alltraitpairs[ij],
                     Estimate=object$variance.components[components,ij],
                     StdErr=object$variance.components.se[components,ij],
                     CI95lo=ci95lo,CI95hi=ci95hi,row.names=components)
        ctables[[count]] <- ctable
      }
    }
  }
  else {  # not bytrait

    ctables <- vector("list",length(components))  # one table per component
    count <- 0
    for(i in components) {
       count <- count + 1
       ci95lo <- object$variance.components[i,traitpairs] - 1.96 * object$variance.components.se[i,traitpairs]
       ci95hi <- object$variance.components[i,traitpairs] + 1.96 * object$variance.components.se[i,traitpairs]
       ctable <- data.frame(Component=i,
                    Estimate=object$variance.components[i,traitpairs],
                    StdErr=object$variance.components.se[i,traitpairs],
                    CI95lo=ci95lo,CI95hi=ci95hi,row.names=traitpairs)
       ctables[[count]] <- ctable
    }
  }
  retobj <- list(ctables=ctables,traits=traits, components=components, bytrait=bytrait, gls=gls, digits=digits)

  if(gls) {
  if(bytrait) {

    gctables <- vector("list",l*l)  # one table per traitpair
    count <- 0
    for(i in traits) {
      for(j in traits) {
        traitpair <- paste(i,":",j,sep="",collapse=NULL)
        ij <- match(traitpair,alltraitpairs)
        count <- count + 1
        ci95lo <- object$gls$variance.components[components,ij] - 1.96 * object$gls$variance.components.se[components,ij]
        ci95hi <- object$gls$variance.components[components,ij] + 1.96 * object$gls$variance.components.se[components,ij]
        ctable <- data.frame(Traitpair=alltraitpairs[ij],
                     Estimate=object$gls$variance.components[components,ij],
                     StdErr=object$gls$variance.components.se[components,ij],
                     CI95lo=ci95lo,CI95hi=ci95hi,row.names=components)
        gctables[[count]] <- ctable
      }
    }
  }
  else {  # not bytrait

    gctables <- vector("list",length(components))  # one table per component
    count <- 0
    for(i in components) {
       count <- count + 1
       ci95lo <- object$gls$variance.components[i,traitpairs] - 1.96 * object$gls$variance.components.se[i,traitpairs]
       ci95hi <- object$gls$variance.components[i,traitpairs] + 1.96 * object$gls$variance.components.se[i,traitpairs]
       ctable <- data.frame(Component=i,
                    Estimate=object$gls$variance.components[i,traitpairs],
                    StdErr=object$gls$variance.components.se[i,traitpairs],
                    CI95lo=ci95lo,CI95hi=ci95hi,row.names=traitpairs)
       gctables[[count]] <- ctable
    }
  }
  retobj <- list(ctables=ctables,gctables=gctables,traits=traits, components=components, bytrait=bytrait, gls=gls, digits=digits)
  }

  retobj$call <- match.call()
  class(retobj) <- "csummary.dmm"
  return(retobj)
  }  #end else nonspecific
}
