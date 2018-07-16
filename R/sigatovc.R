sigatovc <-
function(siga,vsiga,sesiga,am,nsf)
# sigatovc()  -  map siga[,] components into a list vc 
#                with one element per phencovclass
{

    # make a vector of phencovclasses ( classes for phenotypic covariance)
    phencovclasses <- phenclasses(am$comcodes)

    # setup list of sets of components- one element for each phencovclass
    #                                 - each element a class version of siga
    vc <- vector("list",length=length(phencovclasses) )
    var.vc <- vector("list",length=length(phencovclasses) )
    se.vc <- vector("list",length=length(phencovclasses) )

    for (ic in 1 : length(phencovclasses) ) {
      icno <- 0  # subscript in  extracted subset
      invc <- rep(F,nrow(siga))
      navc <- rep(F,nrow(siga))
      colonp <- which(strsplit(phencovclasses[ic],"")[[1]] == ":")
      for(vcno in 1:nrow(siga)) {  # look at each estimated component in turn
        # does this component belong to current phencovclass
        colons <- which(strsplit(rownames(siga)[vcno],"")[[1]] == ":")
        if(length(colons) == 0 ) {  # ie a nonspecific component
          icno <- icno + 1
          invc[vcno] <- T
          naspec <- rep(F,nsf)
          if(nsf > 0) {
            for(kf in 1:nsf) {
              class1 <- substr(phencovclasses[ic],colonp[(kf-1)*3+1]+1,colonp[(kf-1)*3+2]-1)
              if(kf < nsf) {
               end2 <- colonp[(kf-1)*3+3]-1
              }
              else if ( kf == nsf) {
               end2 <- nchar(phencovclasses[ic])
              }
              class2 <- substr(phencovclasses[ic],colonp[(kf-1)*3+2]+1,end2)
  
              if(class1 == class2) {
                naspec[kf] <- T
              }
              else {
                naspec[kf] <- F
              }
            }
          }
          navc[vcno] <- all(naspec)
        }

        else {   # a specific component
          # check if rownames(siga)[ic] chars 1 to 3rd : -1 , matches phencopvclasses[ic] chars 1 to 3rd : -1, or 3rd : + 1 to end ( 2 factor case)
          compclass <- substr(rownames(siga)[vcno],1,colons[3]-1)
          phenclassset <- rep(" ",nsf)
          if (nsf > 1) {
            end <- colonp[3] - 1
          }
          else if(nsf == 1) {
            end <- nchar(phencovclasses[1])
          }
          phenclassset[1] <- substr(phencovclasses[ic],1,end)
          if(nsf > 1) {
            for(kf in 2:nsf) {
              if(kf < nsf) {
                end <- colons[kf*3]-1
              }
              else if(kf == nsf) {
                end <- nchar(phencovclasses[ic])
              }
              phenclassset[kf] <- substr(phencovclasses[ic],colonp[(kf-1)*3]+1,end)
            }
          }
          if(any(!is.na(match(phenclassset,compclass)))) {
            icno <- icno + 1
            invc[vcno] <- T
            navc[vcno] <- T
          }
        }
      }
      vc[[ic]] <- matrix(siga[invc,],icno,ncol(siga))
      var.vc[[ic]] <- matrix(vsiga[rep(invc,ncol(siga)),rep(invc,ncol(siga))],icno*ncol(siga),icno*ncol(siga))
      se.vc[[ic]] <- matrix(sesiga[invc, ],icno,ncol(sesiga))
      for ( exvcno in 1:icno) {
        navcindex <- which(invc)[exvcno]
        if (!navc[navcindex]) {  # need a mapping of icno to vcno
        vc[[ic]][exvcno,] <- rep(NA,ncol(siga))
        varindex <- rep(0,am$l * am$l)
        i <- 0
        for ( tpno in 1 : (am$l * am$l)) {
          i <- i + 1
          varindex[i] <- (tpno - 1) * icno + exvcno
        }
        var.vc[[ic]][varindex,varindex] <- NA
        se.vc[[ic]][exvcno,] <- rep(NA,ncol(sesiga))
        }
      }
      rownames(vc[[ic]]) <- rownames(siga)[invc]
      colnames(vc[[ic]]) <- colnames(siga)
      rownames(se.vc[[ic]]) <- rownames(vc[[ic]])
      colnames(se.vc[[ic]]) <- colnames(vc[[ic]])
      rownames(var.vc[[ic]]) <- rownames(vsiga[rep(invc,ncol(siga)),rep(invc,ncol(siga))])
      colnames(var.vc[[ic]]) <- rownames(var.vc[[ic]])
    }  # end loop over ic
    names(vc) <- phencovclasses
    names(se.vc) <- names(vc)
    names(var.vc) <- names(vc)

#   save long and short rownames
    rownames.vc.long <- vector("list",length=length(phencovclasses))
    rownames.vc.short <- vector("list",length=length(phencovclasses))
    rownames.se.vc.long <- vector("list",length=length(phencovclasses))
    rownames.se.vc.short <- vector("list",length=length(phencovclasses))
    rownames.var.vc.long <- vector("list",length=length(phencovclasses))
    rownames.var.vc.short <- vector("list",length=length(phencovclasses))
    for ( ic in 1:length(phencovclasses)) {
      rownames.vc.long[[ic]] <- rownames(vc[[ic]])
      rownames.vc.short[[ic]] <- genericvcnames(rownames(vc[[ic]]))
      rownames.se.vc.long[[ic]] <- rownames(se.vc[[ic]])
      rownames.se.vc.short[[ic]] <- genericvcnames(rownames(se.vc[[ic]]))
      rownames.var.vc.long[[ic]] <- rownames(var.vc[[ic]])
      rownames.var.vc.short[[ic]] <- genericvarvcnames(rownames(var.vc[[ic]]))
    }
    names(rownames.vc.long) <- phencovclasses
    names(rownames.vc.short) <- phencovclasses
    names(rownames.se.vc.long) <- phencovclasses
    names(rownames.se.vc.short) <- phencovclasses
    names(rownames.var.vc.long) <- phencovclasses
    names(rownames.var.vc.short) <- phencovclasses
#
#   set rownames to short version
    for ( ic in 1:length(phencovclasses)) {
      rownames(vc[[ic]])  <- rownames.vc.short[[ic]]
      rownames(se.vc[[ic]]) <- rownames.se.vc.short[[ic]]
      rownames(var.vc[[ic]]) <- rownames.var.vc.short[[ic]]
      colnames(var.vc[[ic]]) <- rownames.var.vc.short[[ic]]
    }
# 

  outlist <- list(vc=vc,se.vc=se.vc,var.vc=var.vc,phencovclasses=phencovclasses,rownames.vc.short=rownames.vc.short,rownames.vc.long=rownames.vc.long,rownames.se.vc.short=rownames.se.vc.short,rownames.se.vc.long=rownames.se.vc.long,rownames.var.vc.short=rownames.var.vc.short,rownams.var.vc.long=rownames.var.vc.long)
  return(outlist)
}
