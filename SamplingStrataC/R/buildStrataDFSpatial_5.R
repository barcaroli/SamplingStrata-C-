# ----------------------------------------------------
# Function to produce the "strata" dataframe
# starting from the available sampling frame
# using De Gruijter et al formula 2 (2019)
# It works only with optimizeStrataSpatial
# that is equal to optimizeStrata2 but calls this 
# new function
# Author: Giulio Barcaroli
# Date: April 2019
# ----------------------------------------------------
buildStrataDFSpatial <- function(dataset,
                                  fitting=c(1),
                                  range=c(0),
                                  kappa=3,
                                  progress=FALSE,
                                  verbose=FALSE) {
#---------------------------------------------
  # covar <- ifelse(is.nan(covar),0,covar)
  # covar <- as.matrix(covar,nrow=length(pred),ncol=length(pred))
  # pred <- as.data.frame(pred)
  # covar <- as.data.frame(covar)
  # dist <- sqrt((outer(dataset$LON,dataset$LON,"-"))^2+(outer(dataset$LAT,dataset$LAT,"-"))^2)
  #---------------------------------------------
  colnames(dataset) <- toupper(colnames(dataset))
  # if (is.factor(dataset$DOMAINVALUE)) levels(dataset$DOMAINVALUE) <- levels(droplevels(dataset$DOMAINVALUE))
  nvarX <- length(grep("X", names(dataset)))
  nvarY <- length(grep("Y", names(dataset)))
  if (verbose == TRUE) {
    cat("\nComputations are being done on population data\n")
  }
  
    #---------------------------------------------------------     
    numdom <- length(unique(dataset$DOMAINVALUE))
    stratatot <- NULL
    # create progress bar
    if (progress == TRUE) pb <- txtProgressBar(min = 0, max = numdom, style = 3)
    # begin domains cycle
    dataset$DOMAINVALUE <- as.numeric(dataset$DOMAINVALUE)
    for (d in unique(dataset$DOMAINVALUE)) {
      if (progress == TRUE) Sys.sleep(0.1)
      # update progress bar
      if (progress == TRUE) setTxtProgressBar(pb, d)
      dom <- d
		  domain <- dataset[dataset$DOMAINVALUE == dom, ]
      listX <- NULL
      namesX <- NULL
      for (i in 1:nvarX) {
        name <- paste("X", i, sep = "")
        namesX <- cbind(namesX, name)
        if (i < nvarX) 
            listX <- paste(listX, "domain$X", i, ",", sep = "") else listX <- paste(listX, "domain$X", i, sep = "")
      }
      listM <- NULL
      listS <- NULL
      for (i in 1:nvarY) {
          listM <- paste(listM, "M", i, ",", sep = "")
          listS <- paste(listS, "S", i, ",", sep = "")
      }
      stmt <- paste("domain$STRATO <- as.factor(paste(",listX,",sep='*'))", sep = "")
      eval(parse(text = stmt))
      for (i in 1:nvarY) {
        # stmt <- paste("var <- dataset$VAR",i,sep="")
        # eval(parse(text = stmt))
        STRATO <- NULL
        Y <- NULL
        stmt <- paste("Y <- domain$Y", i, "[!is.na(domain$Y", 
            i, ")]", sep = "")
        eval(parse(text = stmt))
        stmt <- paste("STRATO <- domain$STRATO[!is.na(domain$Y", 
            i, ")]", sep = "")
        eval(parse(text = stmt))
        STRATO <- factor(STRATO)
        # Computation of M and S  --------------------------
        dataset$STRATO <- STRATO
        stmt <- paste("M",i," <- NULL",sep="")
        eval(parse(text=stmt))
        stmt <- paste("S",i," <- NULL",sep="")
        eval(parse(text=stmt))
        stmt <- paste("Y <- domain$Y", i, "[!is.na(domain$Y", 
                      i, ")]", sep = "")
        eval(parse(text = stmt))
        stmt <- paste("M", i, " <- tapply(Y,STRATO,mean) ", sep = "")
        eval(parse(text = stmt))
        # stmt <- paste("zz <- outer(dataset$Y",i,",dataset$Y",i,",'-')^2",sep="")
        # eval(parse(text = stmt))
        l.split <- split(dataset, dataset$STRATO, drop = TRUE)
        sd <- sapply(l.split, function(df) stdev_Rcpp(df,i,fitting[i],range[i],kappa))
        # sd <- as.numeric(mclapply(l.split, function(df) stdev_Rcpp(df,i,fitting,range,kappa)))
        stmt <- paste("S", i, " <- sd ", sep = "")
        eval(parse(text = stmt))
        # for (j in (1:length(levels(STRATO)))) {
        #   strat <- levels(STRATO)[j]
        #   stmt <- paste("zz <- outer(dataset$Y",i,",dataset$Y",i,",'-')^2",sep="")
        #   eval(parse(text = stmt))
        #   sd <- stdev(zz,dist,var,strat,dataset,fitting,range,kappa)
        #   stmt <- paste("S",i,"[",j,"] <- sd",sep="")
        #   eval(parse(text=stmt))
        # }
        # ------------------------------------------------------------------------
        stmt <- paste("stratirid <- unlist(attr(M", i, ",'dimnames'))", 
            sep = "")
        eval(parse(text = stmt))
        strati <- data.frame(X1 = levels(domain$STRATO), stringsAsFactors = TRUE)
        stmt <- paste("m <- data.frame(cbind(X1=stratirid,X2=M", 
            i, "), stringsAsFactors = TRUE)", sep = "")
        eval(parse(text = stmt))
        m <- merge(strati, m, by = c("X1"), all = TRUE)
        m$X2 <- as.character(m$X2)
        m$X2 <- as.numeric(m$X2)
        m$X2 <- ifelse(is.na(m$X2), 0, m$X2)
        stmt <- paste("M", i, " <- m$X2", sep = "")
        eval(parse(text = stmt))
        stmt <- paste("s <- data.frame(cbind(X1=stratirid,X2=S", 
            i, "), stringsAsFactors = TRUE)", sep = "")
        eval(parse(text = stmt))
        s <- merge(strati, s, by = c("X1"), all = TRUE)
        s$X2 <- as.character(s$X2)
        s$X2 <- as.numeric(s$X2)
        s$X2 <- ifelse(is.na(s$X2), 0, s$X2)
        stmt <- paste("S", i, " <- s$X2", sep = "")
        eval(parse(text = stmt))
      }
        N <- as.numeric(table(domain$STRATO))
        STRATO <- domain$STRATO
        COST <- rep(1, length(levels(domain$STRATO)))
        CENS <- rep(0, length(levels(domain$STRATO)))
        DOM1 <- rep(as.character(dom), length(levels(domain$STRATO)))
        stmt <- paste("strata <- as.data.frame(cbind(STRATO=levels(STRATO),N,", 
            listM, listS, "COST,CENS,DOM1), stringsAsFactors = TRUE)")
        eval(parse(text = stmt))
        for (i in 1:nvarX) {
            stmt <- paste("strata$X", i, " <- rep(0, length(levels(domain$STRATO)))", 
                sep = "")
            eval(parse(text = stmt))
        }
        strata$STRATO <- as.character(strata$STRATO)
        for (i in 1:nrow(strata)) {
            strata[i, c(namesX)] <- unlist(strsplit(strata$STRATO[i], 
                "\\*"))
        }
        stratatot <- rbind(stratatot, strata)
    }  # end domain cycle
    if (progress == TRUE) close(pb)
    colnames(stratatot) <- toupper(colnames(stratatot))
    stratatot$DOM1 <- as.factor(stratatot$DOM1)
    # write.table(stratatot, "strata.txt", quote = FALSE, sep = "\t",
    #             dec = ".", row.names = FALSE)
    # stratatot <- read.delim("strata.txt")
    # unlink("strata.txt")
    options("scipen"=100)
    indx <- sapply(stratatot, is.factor)
    stratatot[indx] <- lapply(stratatot[indx], function(x) as.numeric(as.character(x)))
    for (j in (1:nvarX)) {
      stmt <- paste("stratatot$X",j," <- as.numeric(stratatot$X",j,")",sep="")
      eval(parse(text=stmt))
    }
      for (j in (1:nrow(stratatot))) {
        stmt <- paste("stratatot$M",i,"[j] <- ifelse(stratatot$M",i,"[j] == 0,0.000000000000001,stratatot$M",i,"[j])",sep="")
        eval(parse(text=stmt))
      }
    # }
    # if (writeFiles == TRUE )
    # write.table(stratatot, "strata.txt", quote = FALSE, sep = "\t", 
    #     dec = ".", row.names = FALSE)
    # stratatot <- read.delim("strata.txt")
    if (verbose == TRUE) {
      cat("\nNumber of strata: ",nrow(stratatot))
      cat("\n... of which with only one unit: ",sum(stratatot$N==1))
    }
    return(stratatot)
}
