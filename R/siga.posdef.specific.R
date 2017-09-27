siga.posdef.specific <-
function(siga, am, ctable,varopt="both",covopt="bound"){
#  siga.posdef.specific() - make each matrix which is a row of siga positive definite (if a var either nonspecific or w/n class)
#                         - make rows which are a cross-class var ( ie a cov) make each 2x2 cov matrix pd
#                         - for rows which are cross-effect and w/n class cov's make each 2x2 cov matrix pd so correlation in bounds
#                         - for rows which are cross-effect and cross-class make each 2x2 matrix pd
    smalleig <- 1e-09

#  vars - nonspecific or within class
    for (i in 1:nrow(siga)) {
     if(!is.cecov(rownames(siga)[i],ctable$allcov)) {  # not a crosseffect cov
       # check  to see if is within class or cross class component
        if(!is.specific(rownames(siga)[i]) | is.withinclass(rownames(siga)[i])){  # nonspecific or within-class 
          # variance - check if pd
          if(any(is.na(siga[i,]))) {
            thisvcna <- T
          }
          else {
            thisvcna <- F
          }
          if(!thisvcna) {
            if(varopt == "eigen" || varopt == "both") {
              eigvc <- eigen(matrix(siga[i, ], am$l, am$l),symmetric=TRUE)
              neg <- FALSE
              for (j in 1:am$l) {
                if (eigvc$values[j] <= 0) {
                    eigvc$values[j] <- smalleig
                    neg <- TRUE
                }
              }
              if (neg) {
                siga[i, ] <- matrix(eigvc$vectors %*% diag(eigvc$values, 
                  nrow = length(eigvc$values)) %*% solve(eigvc$vectors), 
                  1, am$l * am$l)
              }
            }
            if(varopt == "nearPD" || varopt == "both"){
             pdmat <- nearPD(matrix(siga[i, ], am$l, am$l),ensureSymmetry=T)
            #overwrite siga
             siga[i, ] <- matrix(pdmat$mat, 1, am$l * am$l)
            }
          } # end NA check on siga
        }  # end if within or nonspecific
      }  # end if var
    } # end for

# vars  specific cross-class
    for  (i in 1:nrow(siga)) {
      if(!is.cecov(rownames(siga)[i],ctable$allcov)) {  # not a crosseffect cov
       # check to see if is within class or cross class component
        if(is.crossclass(rownames(siga)[i])) { # cross class variance

          # need to find corresponding variances
          varlist <- match.vars(rownames(siga)[i])
          var1 <- siga[varlist$var1name,]
          var2 <- siga[varlist$var2name,]
#         var1 <- vc[[varlist$var1class]][varlist$genericname,]
#         var2 <- vc[[varlist$var2class]][varlist$genericname,]

         for( j in 1 : am$l) {
           jj <- (j-1)*am$l + j
           for(k in 1 : am$l) {
             kk <- (k-1)*am$l + k
             jk <- (j-1) * am$l + k
             if(covopt == "nearPD") {
               covmat <- matrix(c(var1[jj],siga[i,jk],siga[i,jk],var2[kk]),2,2)
               covmat <- nearPD(covmat,keepDiag=T)$mat
               # overwrite siga[i,jk],siga[c1,jj],siga[c2,kk]
               siga[i,jk] <- covmat[1,2]
             }
             else if(covopt == "bound"){
               clim <- sqrt(var1[jj]*var2[kk])
               if(siga[i,jk] < 0) {
                 siga[i,jk] <- max(siga[i,jk],-clim)
               }
               else if(siga[i,jk] > 0) {
                 siga[i,jk] <- min(siga[i,jk],clim)
               }
             }
           }
         }
        }  # end if cross-class-cov
      }  # end if is a var
    }  # end for

#  cross-effect cov - nonspecific or within class
    for (i in 1:nrow(siga)) {
      if(is.cecov(rownames(siga)[i],ctable$allcov)) {  # is a crosseffect cov
        # cross-effect covariance - keep correlation in bounds
        if(!is.specific(rownames(siga)[i]) | is.withinclass(rownames(siga)[i])){  # nonspecific or within  class

          varlist <- match.crosseffect.vars(rownames(siga)[i])
          var1 <- siga[varlist$var1name,]
          var2 <- siga[varlist$var2name,]
#         var1 <- vc[[varlist$var1class]][varlist$genericname1,]
#         var2 <- vc[[varlist$var2class]][varlist$genericname2,]
          c1 <- varlist$var1name
          c2 <- varlist$var2name

         for( j in 1 : am$l) {
           jj <- (j-1)*am$l + j
           for(k in 1 : am$l) {
             kk <- (k-1)*am$l + k
             jk <- (j-1) * am$l + k
             if(covopt == "nearPD") {
               covmat <- matrix(c(siga[c1,jj],siga[i,jk],siga[i,jk],siga[c2,kk]),2,2)
               covmat <- nearPD(covmat,keepDiag=T)$mat
               # overwrite siga[i,jk],siga[c1,jj],siga[c2,kk]
               siga[i,jk] <- covmat[1,2]
#              siga[c1,jj] <- covmat[1,1]
#              siga[c2,kk] <- covmat[2,2]
             }
             else if(covopt == "bound"){
               clim <- sqrt(siga[c1,jj]*siga[c2,kk])
               if(siga[i,jk] < 0) {
                 siga[i,jk] <- max(siga[i,jk],-clim)
               }
               else if(siga[i,jk] > 0) {
                 siga[i,jk] <- min(siga[i,jk],clim)
               }
             }
           }
         }
        } # end w/n class cov or nonspecific
      }  # end if is it a cov
    }  # end for

#  cross-class cross-effect cov
    for (i in 1:nrow(siga)) {
      if(is.cecov(rownames(siga)[i],ctable$allcov)) {  # is a crosseffect cov
        if(is.crossclass(rownames(siga)[i])) {  # is cross-class and a cov

          varlist <- match.crosseffect.vars(rownames(siga)[i])
          var1 <- siga[varlist$var1name,]
          var2 <- siga[varlist$var2name,]
#         var1 <- vc[[varlist$var1class]][varlist$genericname1,]
#         var2 <- vc[[varlist$var2class]][varlist$genericname2,]

         for( j in 1 : am$l) {
           jj <- (j-1)*am$l + j
           for(k in 1 : am$l) {
             kk <- (k-1)*am$l + k
             jk <- (j-1) * am$l + k
             if(covopt == "nearPD") {
               covmat <- matrix(c(var1[jj],siga[i,jk],siga[i,jk],var2[kk]),2,2)
               covmat <- nearPD(covmat,keepDiag=T)$mat
               # overwrite siga[i,jk],siga[c1,jj],siga[c2,kk]
               siga[i,jk] <- covmat[1,2]
             }
             else if(covopt == "bound"){
               clim <- sqrt(var1[jj]*var2[kk])
               if(siga[i,jk] < 0) {
                 siga[i,jk] <- max(siga[i,jk],-clim)
               }
               else if(siga[i,jk] > 0) {
                 siga[i,jk] <- min(siga[i,jk],clim)
               }
             }
           }
         }
        } # end if cross class cross effect cov
      } # end if cross-effect cov
    } # end for i

    return(siga)
}
