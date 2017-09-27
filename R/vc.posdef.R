vc.posdef <-
function(vc, ic, rownames.vc.long, siga, am, ctable,varopt="both",covopt="bound"){
#  vc.posdef() - make each matrix which is a row of vc positive definite (if a var)
    smalleig <- 1e-09
    thisvc <- vc[[ic]]

#  vars - nonspecific or within class
    for (i in 1:nrow(thisvc)) {
      if(any(!is.na(match(rownames(thisvc)[i],ctable$allvar)))){  # is it a var?
       # check names(thisvc) to see if is within class or cross class component
       # only within class components are variances
        if(is.var(rownames.vc.long[[ic]][i]) ) {  #  within-class variance or nonspecific
          # variance - check if pd
          if(any(is.na(thisvc[i,]))) {
            thisvcna <- T
          }
          else {
            thisvcna <- F
          }
          if(!thisvcna) {
            if(varopt == "eigen" || varopt == "both") {
              eigvc <- eigen(matrix(thisvc[i, ], am$l, am$l),symmetric=TRUE)
              neg <- FALSE
              for (j in 1:am$l) {
                if (eigvc$values[j] <= 0) {
                    eigvc$values[j] <- smalleig
                    neg <- TRUE
                }
              }
              if (neg) {
                thisvc[i, ] <- matrix(eigvc$vectors %*% diag(eigvc$values, 
                  nrow = length(eigvc$values)) %*% solve(eigvc$vectors), 
                  1, am$l * am$l)
              }
            }
            if(varopt == "nearPD" || varopt == "both"){
             pdmat <- nearPD(matrix(thisvc[i, ], am$l, am$l),ensureSymmetry=T)
            #overwrite thisvc
             thisvc[i, ] <- matrix(pdmat$mat, 1, am$l * am$l)
            }
          } # end NA check on thisvc
        }  # end if within or nonspecific
      }  # end if var
    } # end for
    vc[[ic]] <- thisvc

# vars  specific cross-class
    for  (i in 1:nrow(thisvc)) {
      if(any(!is.na(match(rownames(thisvc)[i],ctable$allvar)))){  # is it a var?
       # check names(thisvc) to see if is within class or cross class component
       # only within class components are variances
        if(!is.var(rownames.vc.long[[ic]][i]) ) {  #  cross-class variance

          # need to find corresponding variances
          varlist <- match.vars(rownames.vc.long[[ic]][i])
#         var1 <- siga[varlist$var1name,]
#         var2 <- siga[varlist$var2name,]
          var1 <- vc[[varlist$var1class]][varlist$genericname,]
          var2 <- vc[[varlist$var2class]][varlist$genericname,]

         for( j in 1 : am$l) {
           jj <- (j-1)*am$l + j
           for(k in 1 : am$l) {
             kk <- (k-1)*am$l + k
             jk <- (j-1) * am$l + k
             if(covopt == "nearPD") {
               covmat <- matrix(c(var1[jj],thisvc[i,jk],thisvc[i,jk],var2[kk]),2,2)
               covmat <- nearPD(covmat,keepDiag=T)$mat
               # overwrite thisvc[i,jk],thisvc[c1,jj],thisvc[c2,kk]
               thisvc[i,jk] <- covmat[1,2]
             }
             else if(covopt == "bound"){
               clim <- sqrt(var1[jj]*var2[kk])
               if(thisvc[i,jk] < 0) {
                 thisvc[i,jk] <- max(thisvc[i,jk],-clim)
               }
               else if(thisvc[i,jk] > 0) {
                 thisvc[i,jk] <- min(thisvc[i,jk],clim)
               }
             }
           }
         }
        }  # end if cross-class-cov
      }  # end if is a var
    }  # end for
    vc[[ic]] <- thisvc

#  cross-effect cov - nonspecific or within class
    for (i in 1:nrow(thisvc)) {
      if(!any(!is.na(match(rownames(thisvc)[i],ctable$allvar)))){  # is it a cov?
        # cross-effect covariance - keep correlation in bounds
        if(is.var(rownames.vc.long[[ic]][i])) { # w/n class or nonspecific
          if(rownames(thisvc)[i] == "CovG(Ia,Ma)") {
            # work out the c1x & c2y nos for Ia and Ma - ie for VarG(Ia) & VarG(Ma)
            c1 <- match("VarG(Ia)", rownames(thisvc))
            c2 <- match("VarG(Ma)", rownames(thisvc))
          }
          if(rownames(thisvc)[i] == "CovG(Ma,Ia)") {
            c1 <- match("VarG(Ma)", rownames(thisvc))
            c2 <- match("VarG(Ia)", rownames(thisvc))
          }
          if(rownames(thisvc)[i] == "CovG(Id,Md)") {
            c1 <- match("VarG(Id)", rownames(thisvc))
            c2 <- match("VarG(Md)", rownames(thisvc))
          }
          if(rownames(thisvc)[i] == "CovG(Md,Id)") {
            c1 <- match("VarG(Md)", rownames(thisvc))
            c2 <- match("VarG(Id)", rownames(thisvc))
          }
          if(rownames(thisvc)[i] == "CovG(Ia:a,Ma:a)") {
            c1 <- match("VarG(Ia:a)", rownames(thisvc))
            c2 <- match("VarG(Ma:a)", rownames(thisvc))
          }
          if(rownames(thisvc)[i] == "CovG(Ma:a,Ia:a)") {
            c1 <- match("VarG(Ma:a)", rownames(thisvc))
            c2 <- match("VarG(Ia:a)", rownames(thisvc))
          }
          if(rownames(thisvc)[i] == "CovG(Ia:d,Ma:d)") {
            c1 <- match("VarG(Ia:d)", rownames(thisvc))
            c2 <- match("VarG(Ma:d)", rownames(thisvc))
          }
          if(rownames(thisvc)[i] == "CovG(Ma:d,Ia:d)") {
            c1 <- match("VarG(Ma:d)", rownames(thisvc))
            c2 <- match("VarG(Ia:d)", rownames(thisvc))
          }
          if(rownames(thisvc)[i] == "CovG(Id:d,Md:d)") {
            c1 <- match("VarG(Id:d)", rownames(thisvc))
            c2 <- match("VarG(Md:d)", rownames(thisvc))
          }
          if(rownames(thisvc)[i] == "CovG(Md:d,Id:d)") {
            c1 <- match("VarG(Md:d)", rownames(thisvc))
            c2 <- match("VarG(Id:d)", rownames(thisvc))
          }
          if(rownames(thisvc)[i] == "CovGs(Ia,Ma)") {
            c1 <- match("VarGs(Ia)", rownames(thisvc))
            c2 <- match("VarGs(Ma)", rownames(thisvc))
          }
          if(rownames(thisvc)[i] == "CovGs(Ma,Ia)") {
            c1 <- match("VarGs(Ma)", rownames(thisvc))
            c2 <- match("VarGs(Ia)", rownames(thisvc))
          }
          if(rownames(thisvc)[i] == "CovE(I,M)") {
            c1 <- match("VarE(I)", rownames(thisvc))
            c2 <- match("VarE(M)", rownames(thisvc))
          }
          if(rownames(thisvc)[i] == "CovE(M,I)") {
            c1 <- match("VarE(M)", rownames(thisvc))
            c2 <- match("VarE(I)", rownames(thisvc))
          }
          if(rownames(thisvc)[i] == "CovE(I,M&!C)") {
            c1 <- match("VarE(I)", rownames(thisvc))
            c2 <- match("VarE(M&!C)", rownames(thisvc))
          }
          if(rownames(thisvc)[i] == "CovE(M&!C,I)") {
            c1 <- match("VarE(M&!C)", rownames(thisvc))
            c2 <- match("VarE(I)", rownames(thisvc))
          } 
          if(rownames(thisvc)[i] == "CovE(I,M&C)") {
            c1 <- match("VarE(I)", rownames(thisvc))
            c2 <- match("VarE(M&C)", rownames(thisvc))
          }
          if(rownames(thisvc)[i] == "CovE(M&C,I)") {
            c1 <- match("VarE(M&C)", rownames(thisvc))
            c2 <- match("VarE(I)", rownames(thisvc))
          } 
         for( j in 1 : am$l) {
           jj <- (j-1)*am$l + j
           for(k in 1 : am$l) {
             kk <- (k-1)*am$l + k
             jk <- (j-1) * am$l + k
             if(covopt == "nearPD") {
               covmat <- matrix(c(thisvc[c1,jj],thisvc[i,jk],thisvc[i,jk],thisvc[c2,kk]),2,2)
               covmat <- nearPD(covmat,keepDiag=T)$mat
               # overwrite thisvc[i,jk],thisvc[c1,jj],thisvc[c2,kk]
               thisvc[i,jk] <- covmat[1,2]
#              thisvc[c1,jj] <- covmat[1,1]
#              thisvc[c2,kk] <- covmat[2,2]
             }
             else if(covopt == "bound"){
               clim <- sqrt(thisvc[c1,jj]*thisvc[c2,kk])
               if(thisvc[i,jk] < 0) {
                 thisvc[i,jk] <- max(thisvc[i,jk],-clim)
               }
               else if(thisvc[i,jk] > 0) {
                 thisvc[i,jk] <- min(thisvc[i,jk],clim)
               }
             }
           }
         }
        } # end w/n class cov or nonspecific
      }  # end if is it a cov
    }  # end for
    vc[[ic]] <- thisvc

#  cross-class cross-effect cov
    for (i in 1:nrow(thisvc)) {
      if(!any(!is.na(match(rownames(thisvc)[i],ctable$allvar)))){ # is it a cov?
        if(!is.var(rownames.vc.long[[ic]][i]) ) {  #  cross-class cov 

          varlist <- match.crosseffect.vars(rownames.vc.long[[ic]][i])
#         var1 <- siga[varlist$var1name,]
#         var2 <- siga[varlist$var2name,]
          var1 <- vc[[varlist$var1class]][varlist$genericname1,]
          var2 <- vc[[varlist$var2class]][varlist$genericname2,]

         for( j in 1 : am$l) {
           jj <- (j-1)*am$l + j
           for(k in 1 : am$l) {
             kk <- (k-1)*am$l + k
             jk <- (j-1) * am$l + k
             if(covopt == "nearPD") {
               covmat <- matrix(c(var1[jj],thisvc[i,jk],thisvc[i,jk],var2[kk]),2,2)
               covmat <- nearPD(covmat,keepDiag=T)$mat
               # overwrite thisvc[i,jk],thisvc[c1,jj],thisvc[c2,kk]
               thisvc[i,jk] <- covmat[1,2]
             }
             else if(covopt == "bound"){
               clim <- sqrt(var1[jj]*var2[kk])
               if(thisvc[i,jk] < 0) {
                 thisvc[i,jk] <- max(thisvc[i,jk],-clim)
               }
               else if(thisvc[i,jk] > 0) {
                 thisvc[i,jk] <- min(thisvc[i,jk],clim)
               }
             }
           }
         }
        } # end if cross class cross effect cov
      } # end if cross-effect cov
    } # end for i

    vc[[ic]] <- thisvc
    return(vc)
}
