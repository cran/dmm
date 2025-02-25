summary.dmm <-
function(object,traitset="all",componentset="all",bytrait=T,fixedgls=F,digits=3, ...)
# summary.dmm() - make summary tables for a dmm object
{
  # check object is a dmm object fully populated
  if(is.null(object$b)) {
    stop("summary.dmm: dmm object does not contain item b:\n:")
  }
  if(is.null(object$siga)) {
    stop("summary.dmm: dmm object does not contain item siga:\n")
  }

   if(traitset[1] == "all"){
    traits <- dimnames(object$b)[[2]][1:ncol(object$b)]
  }
  else {
    traits <- traitset
  }
  traitpairs <- permpaste(traits)
  l <- length(traits)
  alltraitpairs <- permpaste(dimnames(object$b)[[2]])
 
  if(componentset[1] == "all") {
    components <- dimnames(object$siga)[[1]]
  }
  else {
    components <- componentset
  }
  coefs <- dimnames(object$b)[[1]]

  if(bytrait) {
    btables <- vector("list",l)  # one table per trait
    count <- 0
    for(j in traits) {
      count <- count + 1
      ci95lo <- object$b[ ,j] - 1.96 * object$seb[ ,j]
      ci95hi <- object$b[ ,j] + 1.96 * object$seb[ ,j]
      btable <- data.frame(Trait=j, Estimate=object$b[ ,j],
                      StdErr=object$seb[ ,j],CI95lo=ci95lo,CI95hi=ci95hi)
      # save btables as one element of a list of tables
      btables[[count]] <- btable
    }

    ctables <- vector("list",l*l)  # one table per traitpair
    count <- 0
    for(i in traits) {
      for(j in traits) {
        traitpair <- paste(i,":",j,sep="",collapse=NULL)
        ij <- match(traitpair,alltraitpairs)
        count <- count + 1
        ci95lo <- object$siga[components,ij] - 1.96 * object$sesiga[components,ij]
        ci95hi <- object$siga[components,ij] + 1.96 * object$sesiga[components,ij]
        ctable <- data.frame(Traitpair=alltraitpairs[ij],
                   Estimate=object$siga[components,ij],
                   StdErr=object$sesiga[components,ij],
                   CI95lo=ci95lo,CI95hi=ci95hi,row.names=components)
        ctables[[count]] <- ctable
      }
    }
  }  # end if bytrait

  else {  # not bytrait
    btables <- vector("list",length(coefs))  # one table per coefficient
    count <- 0
    for(j in coefs) {
      count <- count + 1
      ci95lo <- object$b[j,traits] - 1.96 * object$seb[j,traits]
      ci95hi <- object$b[j,traits] + 1.96 * object$seb[j,traits]
      btable <- data.frame(Coefficient=j, Estimate=object$b[j,traits],
                      StdErr=object$seb[j,traits],CI95lo=ci95lo,CI95hi=ci95hi,
                      row.names=traits)
      # save btables as one element of a list of tables
      btables[[count]] <- btable
    }

    ctables <- vector("list",length(components))  # one table per component
    count <- 0
    for(i in components) {
       count <- count + 1
       ci95lo <- object$siga[i,traitpairs] - 1.96 * object$sesiga[i,traitpairs]
       ci95hi <- object$siga[i,traitpairs] + 1.96 * object$sesiga[i,traitpairs]
       ctable <- data.frame(Component=i,
                    Estimate=object$siga[i,traitpairs],
                    StdErr=object$sesiga[i,traitpairs],
                    CI95lo=ci95lo,CI95hi=ci95hi,row.names=traitpairs)
       ctables[[count]] <- ctable
    }
  }
  retobj <- list(btables=btables, ctables=ctables,traits=traits, components=components, bytrait=bytrait, fixedgls=fixedgls, digits=digits)

  if(fixedgls) {
  if(bytrait) {
    gbtables <- vector("list",l)  # one table per trait
    count <- 0
    for(j in traits) {
      count <- count + 1
      ci95lo <- object$gls$b[ ,j] - 1.96 * object$gls$seb[ ,j]
      ci95hi <- object$gls$b[ ,j] + 1.96 * object$gls$seb[ ,j]
      btable <- data.frame(Trait=j, Estimate=object$gls$b[ ,j],
                      StdErr=object$gls$seb[ ,j],CI95lo=ci95lo,CI95hi=ci95hi)
      # save btables as one element of a list of tables
      gbtables[[count]] <- btable
    }

    gctables <- vector("list",l*l)  # one table per traitpair
    count <- 0
    for(i in traits) {
      for(j in traits) {
        traitpair <- paste(i,":",j,sep="",collapse=NULL)
        ij <- match(traitpair,alltraitpairs)
        count <- count + 1
        ci95lo <- object$gls$siga[components,ij] - 1.96 * object$gls$sesiga[components,ij]
        ci95hi <- object$gls$siga[components,ij] + 1.96 * object$gls$sesiga[components,ij]
        ctable <- data.frame(Traitpair=alltraitpairs[ij],
                     Estimate=object$gls$siga[components,ij],
                     StdErr=object$gls$sesiga[components,ij],
                     CI95lo=ci95lo,CI95hi=ci95hi,row.names=components)
        gctables[[count]] <- ctable
      }
    }
  }
  else {  # not bytrait
    gbtables <- vector("list",length(coefs))  # one table per coefficient
    count <- 0
    for(j in coefs) {
      count <- count + 1
      ci95lo <- object$gls$b[j,traits] - 1.96 * object$gls$seb[j,traits]
      ci95hi <- object$gls$b[j,traits] + 1.96 * object$gls$seb[j,traits]
      btable <- data.frame(Coefficient=j, Estimate=object$gls$b[j,traits],
                 StdErr=object$gls$seb[j,traits],CI95lo=ci95lo,CI95hi=ci95hi,
                 row.names=traits)
      # save btables as one element of a list of tables
      gbtables[[count]] <- btable
    }

    gctables <- vector("list",length(components))  # one table per component
    count <- 0
    for(i in components) {
       count <- count + 1
       ci95lo <- object$gls$siga[i,traitpairs] - 1.96 * object$gls$sesiga[i,traitpairs]
       ci95hi <- object$gls$siga[i,traitpairs] + 1.96 * object$gls$sesiga[i,traitpairs]
       ctable <- data.frame(Component=i,
                    Estimate=object$gls$siga[i,traitpairs],
                    StdErr=object$gls$sesiga[i,traitpairs],
                    CI95lo=ci95lo,CI95hi=ci95hi,row.names=traitpairs)
       gctables[[count]] <- ctable
    }
  }
  retobj <- list(btables=btables,ctables=ctables,gbtables=gbtables,gctables=gctables,traits=traits, components=components, bytrait=bytrait, fixedgls=fixedgls, digits=digits)
  }

  retobj$call <- match.call()
  class(retobj) <- "summary.dmm"
  return(retobj)
}
