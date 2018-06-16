am.zandrel <-
function (mdf,df, k, l, x, y, cohortparts, components,specific.components,
          relmat,ctable) {
  # am.zandrel()
  # construct the animal model matrices z[r],rel[r] from a data frame (df)
  #  or construct z[] but extract rel[] from a list (mdf)
  # setup am object as a list
  #
  # must specify k = no fixed effects, l = no traits
  #  cohortparts = vector of df ccol nos for parts of cohort defining common env
  #  components = vector listing non-specific animal components to be fitted
  #  specific.components =  components specific to one or more factors
  #  relmat = flags whether inline or withdf
  #
  # determines m = no individuals in df
  #            n = no of individuals with data
  # return all in list object called am

  m <- length(df$Id)  # no of individuals in pedigree
  n <- nrow(y)        # no of individuals with y data and X codes

#  check pedigree in dataframe valid
  if(pedcheck(df) > 0){
   stop("dmm(): pedigree not valid:\n")
  }

  if(is.null(cohortparts)) {
    cohortcode <- NULL
    ncohortcodes <- 1
  }
  else {
    cohortcode <- make.cohort(df,cohortparts)
    cohortcodes <- unique(cohortcode)  # list of codes
    cohortcodes <- cohortcodes[!is.na(cohortcodes)]  # NA's removed
    ncohortcodes <- length(cohortcodes)
#    ncohortcodes <- length(unique(cohortcode)) counts NA !!
    cat("ncohortcodes = ",ncohortcodes,"\n")
  }
# setup R matrices - for e, a, m, ...  effects
  rel <- list(e=NULL,a=NULL,aa=NULL,d=NULL,ad=NULL,dd=NULL,s=NULL)

# setup logical miss[] - flags individuals in pedigree but no Y data
  miss <- match(dimnames(df)[[1]], dimnames(y)[[1]])
#

# non-specific and specific mappings into pedigree
# All specific factors must be in df, but dont have to be in the fixed model
  nnsc <- length(components) # no of non-specific components
  nsf <- length(specific.components) # no of factors with specific components
  cat("No of factors with specific components:",nsf,"\n")
  effnandc <- vector("list",length=nsf) # effect names and codes pasted
  effnames <- names(specific.components) # ie "Sex", "Age", ..
  effcodes <- vector("list",length=nsf)  # effect codes
  comcodes <- vector("list",length=nsf)  # var/cov component codes
  varcodes <- vector("list",length=nsf)  # var component codes
  nsc <- 0   # no of specific components across all factors
  nsvc <- 0  # no of specific var/cov components across all factors
  if(nsf > 0) {
    for (i in 1:nsf) {
      effcodes[[i]] <- levels(df[,effnames[i]])
      comcodes[[i]] <- combpaste(effnames[i],permpaste(effcodes[[i]]))
      varcodes[[i]] <- combpaste(effnames[i],selfpaste(effcodes[[i]]))
      effnandc[[i]] <- combpaste(effnames[i],effcodes[[i]])
      nsc <- nsc + length(effcodes[[i]])
      nsvc <- nsvc + length(comcodes[[i]]) * length(specific.components[[i]])  # includes cross-factor covs
    }
  names(effcodes) <- effnames
  names(comcodes) <- effnames
  names(effnandc) <- effnames
# cat("effcodes:\n")
# print(effcodes)
# cat("comcodes:\n")
# print(comcodes)
# cat("effnandc:\n")
# print(effnandc)
  }
  nc <- nnsc + nsc  # total no of components fitted
  cat("No of non-specific components partitioned:",nnsc,"\n")
  cat("No of factors with specific components:",nsf,"\n")
  cat("No of specific variance components partitioned (per component):",nsc,"\n")
  cat("No of specific variance and covariance components partitioned (per component):",nsvc,"\n")

# Individual effects
#  Z matrices (1 for all I's, 1 for each level of each specific factor)
# Always make individual Z matrices
  nzi <- 1 + nsc
  zi <- vector("list",length = nzi)
  zinames <- c("NS",unlist(effnandc))
  names(zi) <- zinames # name of "NS" is blank

# construct each element ( matrix) of zi and add it to list zi

# zi$NS - incidence of all individuals measured in pedigree - (n x m)
#     - always do z$NS  ( NS = non-specific)
  zi$NS <- matrix(0,n,m)
  inrow <- 0.0  # inrow counts individuals with no missing data(or X codes)
  for (i in 1:m) {   # i indexes individuals by row in df
    if (!is.na(miss[i])) {
      inrow <- inrow + 1   # inrow indexes individuals by row in x,y,z
      for (j in 1:m) {
        zi$NS[inrow,j] <- 0.0
      }
      zi$NS[inrow,i] <- 1.0
    }
  }
#
# setup zi$... for sex specific and fixed effect specific cases
# only do this if is at least one specific factor
  if(nsf > 0) {
    count <- 0
    for(kf in 1:length(effnames)) {
      if(any(!is.na(match(specific.components[[kf]],ctable$ind)))) {
        count <- count + 1
      }
    }
    if(count > 0 ) {
      kl <- 1  # kl is subscript for zinames 
      for (kf in 1:length(effnames)) {   # kf is factor name
        for (lc in 1:length(effcodes[[kf]])) {   # lc is codes for factor kf
          kl <- kl + 1
#         zi[[zinames[kl]]] <- matrix(0,n,m)
          zi[[kl]] <- matrix(0,n,m)
          inrow <- 0
          for(i in 1:m) {
            if (!is.na(miss[i])) {
              inrow <- inrow + 1
              if(df[i,effnames[kf]] == effcodes[[kf]][lc]) {
#               zi[[zinames[kl]]][inrow,i] <- 1
                zi[[kl]][inrow,i] <- 1
              }
            }
          }
        }
      }
    }
  }



# Maternal effects
# DId must be in the df
# zm - incidence of dams of measured individuals in pedigree - ( n x m)
  nzm <- 1 + nsc
  zm <- vector("list",length = nzm)
  names(zm) <- zinames
# Only make zm$NS if a non-specific maternal component present
 if(any(!is.na(match(components,ctable$mat)))){  # non-specific maternal component check
  inrow <- 0.0
  zm$NS <- matrix(0,n,m)
# miss <- match(dimnames(df)[[1]], dimnames(y)[[1]])
  for(i in 1:m) {
    if(!is.na(miss[i])) {
      inrow <- inrow + 1   # inrow indexes individuals by row in x,y,z
      for (j in 1:m) {
        zm$NS[inrow,j] <- 0.0
      }
      if(!is.na(df$DId[i])){
        zm$NS[inrow,df$DId[i]] <- 1.0
      }
    }
  }
 }

# Setup zm$... for sex-specific and fixed effect specific cases
# Only do this if there is a specific maternal component present
  if( nsf > 0) {
    count <- 0
    for(kf in 1:length(effnames)){
      if(any(!is.na(match(specific.components[[kf]],ctable$mat)))){
        count <- count + 1
      }
    }
    if(count > 0) {
      kl <- 1  # kl is subscript for zinames
      for (kf in 1:length(effnames)) {   # kf is factor name
        for (lc in 1:length(effcodes[[kf]])) {   # lc is codes for factor kf
          kl <- kl + 1
          zm[[kl]] <- matrix(0,n,m)
          inrow <- 0
          for(i in 1:m) {
            if (!is.na(miss[i])) {
              inrow <- inrow + 1
              if(df[i,effnames[kf]] == effcodes[[kf]][lc]) { # if Ind has effcode
                if(!is.na(df$DId[i])){ # if dam is present
                  zm[[kl]][inrow,df$DId[i]] <- 1  # set Ind's row and dam's col 1
                }
              }
            }
          }
        }
      }
    }
  }


# Common cohort env effect ( same cohort)
# All cohort factors must be in df and in cohort formula
# zc - incidence of measured individuals in cohorts (n x ncohorts)
  nzc <- 1 + nsc
  zc <- vector("list",length = nzc)
  names(zc) <- zinames
# Only make zc$NS if a non-specific cohort component present
  if(any(!is.na(match(components,ctable$cohort)))){  # specific cohort component check
    cohortcode <- as.factor(cohortcode)
#  aov() method - discarded due to problems with nrow(zc)
#   ce.aov <- aov(df$Ymat ~ -1 + cohortcode,x=T)
#   zc <- ce.aov$x
# match method - similar to zi and zm
    inrow <- 0
    zc$NS <- matrix(0,n,ncohortcodes)
#   miss <- match(dimnames(cohortcode)[[1]], dimnames(y)[[1]])
  # miss gives mth cohortcode[] matches nth y[] value
    for(i in 1:m){
      if(!is.na(miss[i])) {
        inrow <- inrow + 1 # inrow indexes individuals by row in x,y,z
        for(j in 1:ncohortcodes) {
          zc$NS[inrow,j] <- 0
        }
#       if(!is.na(cohortcode[i])){
          whichcohort <- match(cohortcode[i], levels(cohortcode))
          zc$NS[inrow,whichcohort] <- 1.0
#       }
      }
    }

    if(nrow(zi$NS) != nrow(zc$NS)) {
      stop("these must be equal:\n")
    }
  }

# Setup zc$... for sex-specific and fixed effect specific cases
# Only do this if there is a specific cohort component present
  if(nsf > 0) {
    count <- 0
    for(kf in 1:length(effnames)){
      if(any(!is.na(match(specific.components[[kf]],ctable$cohort)))){
        count <- count + 1
      }
    }
    if(count > 0) {
      kl <- 1  # kl is subscript for zinames
      for (kf in 1:length(effnames)) {   # kf is factor name
        for (lc in 1:length(effcodes[[kf]])) {   # lc is codes for factor kf
          kl <- kl + 1
          zc[[kl]] <- matrix(0,n,ncohortcodes)
          inrow <- 0
          for(i in 1:m) {
            if (!is.na(miss[i])) {
              inrow <- inrow + 1
              if(df[i,effnames[kf]] == effcodes[[kf]][lc]) { # if Ind has effcode
                if(!is.na(cohortcode[i])){ # if cohortcode is present
                  whichcohort <- match(cohortcode[i],levels(cohortcode))
                  zc[[kl]][inrow,whichcohort] <- 1  # set Ind's row and cohort's col 1
                }
              }
            }
          }
        }
      }
    }
  }

# Maternal cytoplasmic lines
# zy - incidence of measured individuals in maternal cytoplasmic lines (n x nmatclines)
  nzy <- 1 + nsc
  zy <- vector("list",length = nzy)
  names(zy) <- zinames
# Only make zy$NS if a non-specific cytoplasmic component present
  if(any(!is.na(match(components,ctable$matcg)))){  # specific matc component check
    # check if matccodes are on df as a column called "MLine"
    if(!is.null(df$MLine)) {
      matccode <- df$MLine
    }
    else{
      matccode <- make.matcline(df)
    }
    matccodes <- unique(matccode)
    nmatccodes <- length(matccodes)
    cat("no of maternal line codes = ",nmatccodes,"\n")
# match method - similar to zi and zm
    inrow <- 0
    zy$NS <- matrix(0,n,nmatccodes)
#   miss <- match(dimnames(matccode)[[1]], dimnames(y)[[1]])
  # miss gives mth matccode[] matches nth y[] value
    for(i in 1:m){
      if(!is.na(miss[i])) {
        inrow <- inrow + 1 # inrow indexes individuals by row in x,y,z
        for(j in 1:nmatccodes) {
          zy$NS[inrow,j] <- 0
        }
#       if(!is.na(matccode[i])){
          whichmatcline <- match(matccode[i], matccodes)
          zy$NS[inrow,whichmatcline] <- 1.0
#       }
      }
    }

    if(nrow(zi$NS) != nrow(zy$NS)) {
      stop("these must be equal:\n")
    }
  }

# Setup zy$... for sex-specific and fixed effect specific cases
# Only do this if there is a specific cytoplasmic component present
  if(nsf > 0) {
    count <- 0
    for(kf in 1:length(effnames)){
      if(any(!is.na(match(specific.components[[kf]],ctable$matcg)))){
        count <- count + 1
      }
    }
    if(count > 0) {
      #check if matccode is in df as column called "MLine"
      if(!is.null(df$MLine)) {
        matccode <- df$MLine
      }
      else {
        matccode <- make.matcline(df)
      }
      matccodes <- unique(matccode)
      nmatccodes <- length(matccodes)
      cat("no of maternal line codes = ",nmatccodes,"\n")

      kl <- 1  # kl is subscript for zinames
      for (kf in 1:length(effnames)) {   # kf is factor name
        for (lc in 1:length(effcodes[[kf]])) {   # lc is codes for factor kf
          kl <- kl + 1
          zy[[kl]] <- matrix(0,n,nmatccodes)
          inrow <- 0
          for(i in 1:m) {
            if (!is.na(miss[i])) {
              inrow <- inrow + 1
              if(df[i,effnames[kf]] == effcodes[[kf]][lc]) { # if Ind has effcode
                if(!is.na(matccode[i])){ # if matccode is present
                  whichmatcline <- match(matccode[i],matccodes)
                  zy[[kl]][inrow,whichmatcline] <- 1  # set Ind's row and matcline's col 1
                }
              }
            }
          }
        }
      }
    }
  }


# Paternal Y-chromosome  lines
# zp - incidence of measured individuals in paternal Y-chromosome lines (n x npatylines)
  nzp <- 1 + nsc
  zp <- vector("list",length = nzp)
  names(zp) <- zinames
# Only make zp$NS if a non-specific paternal-Y component present
  if(any(!is.na(match(components,ctable$patyg)))){  # specific paty component check
    # check if patycodes are on df as a column called "PLine"
    if(!is.null(df$PLine)) {
      patycode <- df$PLine
    }
    else{
      patycode <- make.patyline(df)
    }
    patycodes <- unique(patycode)
    npatycodes <- length(patycodes)
    cat("no of paternal line codes = ",npatycodes,"\n")
# match method - similar to zi and zm
    inrow <- 0
    zp$NS <- matrix(0,n,npatycodes)
#   miss <- match(dimnames(patycode)[[1]], dimnames(y)[[1]])
  # miss gives mth patycode[] matches nth y[] value
    for(i in 1:m){
      if(!is.na(miss[i])) {
        inrow <- inrow + 1 # inrow indexes individuals by row in x,y,z
        for(j in 1:npatycodes) {
          zp$NS[inrow,j] <- 0
        }
#       if(!is.na(patycode[i])){
          whichpatyline <- match(patycode[i], patycodes)
          zp$NS[inrow,whichpatyline] <- 1.0
#       }
      }
    }

    if(nrow(zi$NS) != nrow(zp$NS)) {
      stop("these must be equal:\n")
    }
  }

# Setup zp$... for sex-specific and fixed effect specific cases
# Only do this if there is a specific paternal Y-chromosome component present
  if(nsf > 0) {
    count <- 0
    for(kf in 1:length(effnames)){
      if(any(!is.na(match(specific.components[[kf]],ctable$patyg)))){
        count <- count + 1
      }
    }
    if(count > 0) {
      #check if patycode is in df as column called "PLine"
      if(!is.null(df$PLine)) {
        patycode <- df$PLine
      }
      else {
        patycode <- make.patyline(df)
      }
      patycodes <- unique(patycode)
      npatycodes <- length(patycodes)
      cat("no of maternal line codes = ",npatycodes,"\n")

      kl <- 1  # kl is subscript for zinames
      for (kf in 1:length(effnames)) {   # kf is factor name
        for (lc in 1:length(effcodes[[kf]])) {   # lc is codes for factor kf
          kl <- kl + 1
          zp[[kl]] <- matrix(0,n,npatycodes)
          inrow <- 0
          for(i in 1:m) {
            if (!is.na(miss[i])) {
              inrow <- inrow + 1
              if(df[i,effnames[kf]] == effcodes[[kf]][lc]) { # if Ind has effcode
                if(!is.na(patycode[i])){ # if patycode is present
                  whichpatyline <- match(patycode[i],patycodes)
                  zp[[kl]][inrow,whichpatyline] <- 1  # set Ind's row and patyline's col 1
                }
              }
            }
          }
        }
      }
    }
  }



#
# Setup relationship matrices
#
# Additive genetic effects - need rel$a matrix
  count <- 0
  if(nsf > 0) {
    for(kf in 1:length(effnames)) {
      if(any(!is.na(match(specific.components[[kf]],ctable$addg)))){
        count <- count + 1
      }
    }
  }
  if (count > 0 || any(!is.na(match(components,ctable$addg)))) {
    # rel$a - additive relationship matrix
   if(relmat == "inline") {
    rel$a <- am.arel(df)
   }
   else if(relmat == "withdf"){
    rel$a <- as.matrix(mdf$rel$a)
   }
  }

  
# Residual or E effect - use I for residual rel matrix
# rel$e - residual relationship matrix - usually I
# Always make rel$e
  if(relmat == "inline"){
    rel$e <- diag(m)
  }
  else if(relmat == "withdf") {
    rel$e <- as.matrix(mdf$rel$e)
  }

# Non additive genetic effects
  # dominance
  count <- 0
  if(nsf > 0) {
    for(kf in 1:length(effnames)) {
      if(any(!is.na(match(specific.components[[kf]],ctable$domg)))){
        count <- count + 1
      }
    }
  }
  if (count > 0 || any(!is.na(match(components,ctable$domg)))) {
    rel$d <- as.matrix(mdf$rel$d)
  }
  # epistasis-additive
  count <- 0
  if(nsf > 0) {
    for(kf in 1:length(effnames)) {
      if(any(!is.na(match(specific.components[[kf]],ctable$epiaddg)))){
        count <- count + 1
      }
    }
  }
  if (count > 0 || any(!is.na(match(components,ctable$epiaddg)))) {
    rel$aa <- as.matrix(mdf$rel$aa)
  }
  # epistasis-dominance
  count <- 0
  if(nsf > 0) {
    for(kf in 1:length(effnames)) {
      if(any(!is.na(match(specific.components[[kf]],ctable$epidomg)))){
        count <- count + 1
      }
    }
  }
  if (count > 0 || any(!is.na(match(components,ctable$epidomg)))) {
    rel$ad <- as.matrix(mdf$rel$ad)
    rel$dd <- as.matrix(mdf$rel$dd)
  }
  # sex-linked
  count <- 0
  if(nsf > 0) {
    for(kf in 1:length(effnames)) {
      if(any(!is.na(match(specific.components[[kf]],ctable$sexlinaddg)))){
        count <- count + 1
      }
    }
  }
  if (count > 0 || any(!is.na(match(components,ctable$sexlinaddg)))) {
     rel$s <- as.matrix(mdf$rel$s)
  }


# construct the am list object
  am <- list(m=m,n=n,k=k,l=l,v=nnsc + nsvc,x=x,y=y,zi=zi,zm=zm,zc=zc,zy=zy,zp=zp,rel=rel,components=components,specific.components=specific.components,zinames=zinames,effnames=effnames,effcodes=effcodes,effnandc=effnandc,comcodes=comcodes,varcodes=varcodes,nnsc=nnsc,nsc=nsc,nc=nc,nsvc=nsvc)
 
  return(am)
}
