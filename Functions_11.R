########## log marginal likelihood:
log.marginal.likelihood <- function(s2, n, b="default", hypothesis, nsim=1e5) {
 
  J <- length(n)
  if (all(b=="default")) {bF <- 2/n; bB <- (1+1/J)/n} else {bB <- bF <- b}
  nubullet <- sum(bB*n)-J
  tau2bullet <- sum(bB*(n-1)*s2)/nubullet  
  complement <- F
  nearorder <- F

  if (hypothesis=="null") {hypothesis <- paste(as.character(1:length(s2)), collapse="=")}
  if (hypothesis=="order") {hypothesis <- paste(as.character(1:length(s2)), collapse="<")}
  if (hypothesis=="unconstrained") {hypothesis <- paste(as.character(1:length(s2)), collapse=",")}
  if (hypothesis=="near-order") {hypothesis <- paste(as.character(1:length(s2)), collapse="<"); nearorder <- T}
  if (grepl("near", hypothesis)) {hypothesis <- gsub("[near ()]", "", hypothesis); nearorder <- T}      
  if (hypothesis=="complement") {hypothesis <- paste(as.character(1:length(s2)), collapse="<"); complement <- T}
  if (grepl("not", hypothesis)) {
    complement <- T
    hypothesis <- gsub("[not ()]", "", hypothesis)
    if(grepl("r", hypothesis)) hypothesis <- unlist(strsplit(hypothesis, split="r"))
  }
  
  #if (any(duplicated(as.numeric(unlist(strsplit(gsub("[ =<,()]", "", hypothesis), split=""))))) && !complement) {stop("Each variance may only appear once in each hypothesis.")}
                  
  equalities <- gsub("[ ()]", "", hypothesis)      
  equalities <- unlist(strsplit(equalities, split="[,<]"))
  equalities <- lapply(as.list(unique(equalities)), function(x) as.numeric(unlist(strsplit(x, split="="))))  
   
  s2list <- lapply(equalities, function(x) s2[x])
  nlist <- lapply(equalities, function(x) n[x])
  bFlist <- lapply(equalities, function(x) bF[x])
  m <- sapply(equalities, length)
  M <- length(m)                                                            
                                                                                  
  dfF <- dfFb <- SSF <- SSFb <- rep(NA, M)                                         
  for (i in 1:M) {                                                                                  
    dfF[i] <- sum(nlist[[i]])-m[i]
    dfFb[i] <- sum(bFlist[[i]]*nlist[[i]])-m[i]
    SSF[i] <- sum((nlist[[i]]-1)*s2list[[i]])
    SSFb[i] <- sum(bFlist[[i]]*(nlist[[i]]-1)*s2list[[i]])
  }
  dfB <- dfF+nubullet
  SSB <- SSF+nubullet*tau2bullet

  logmbx <- c(0, 0, 0)
  names(logmbx) <- c("FBF", "BBF", "aFBF")

  logmbx[c(1, 3)] <- 1/2*log(prod(unlist(bFlist)))-1/2*sum((1-unlist(bFlist))*unlist(nlist))*log(pi)+sum(lgamma(dfF/2)-lgamma(dfFb/2)-1/2*dfF*log(SSF)+1/2*dfFb*log(SSFb))
  logmbx[2] <- -1/2*log(prod(unlist(nlist)))-1/2*(sum(unlist(nlist))-J)*log(pi)+1/2*M*nubullet*log(nubullet*tau2bullet)-M*lgamma(nubullet/2)+sum(lgamma(dfB/2)-1/2*dfB*log(SSB))

  if (complement) {
  
    inequalities <- lapply(hypothesis, function(x) as.numeric(unlist(strsplit(x, split="<"))))
    variances <- unique(unlist(inequalities))  

    samp <- matrix(NA, nrow=nsim, ncol=M)

    for (i in 1:length(variances)) {samp[, variances[i]] <- SSF[variances[i]]/rchisq(nsim, dfF[variances[i]])}  
    ineqs <- matrix(0, nrow=nsim, ncol=length(inequalities))
    for (i in 1:length(inequalities)) {for (j in 1:(length(inequalities[[i]])-1)) {ineqs[, i] <- ineqs[, i]+(samp[, which(variances==inequalities[[i]][j])]<samp[, which(variances==inequalities[[i]][j+1])])}}
    PF <- sum(apply(ineqs, 1, max)==M-1)/nsim

    for (i in 1:length(variances)) {samp[, variances[i]] <- SSFb[variances[i]]/rchisq(nsim, dfFb[variances[i]])}
    ineqs <- matrix(0, nrow=nsim, ncol=length(inequalities))                                   
    for (i in 1:length(inequalities)) {for (j in 1:(length(inequalities[[i]])-1)) {ineqs[, i] <- ineqs[, i]+(samp[, which(variances==inequalities[[i]][j])]<samp[, which(variances==inequalities[[i]][j+1])])}}
    PFb <- sum(apply(ineqs, 1, max)==M-1)/nsim

    for (i in 1:length(variances)) {samp[, variances[i]] <- SSB[variances[i]]/rchisq(nsim, dfB[variances[i]])}
    ineqs <- matrix(0, nrow=nsim, ncol=length(inequalities))                                    
    for (i in 1:length(inequalities)) {for (j in 1:(length(inequalities[[i]])-1)) {ineqs[, i] <- ineqs[, i]+(samp[, which(variances==inequalities[[i]][j])]<samp[, which(variances==inequalities[[i]][j+1])])}}
    PB <- sum(apply(ineqs, 1, max)==M-1)/nsim

    PBb <- length(inequalities)/factorial(M)

    logmbx <- logmbx+c(log((1-PF)/(1-PFb)), log((1-PB)/(1-PBb)), log((1-PF)/(1-PBb)))
    #if (0%in%c(1-PF, 1-PFb, 1-PB, 1-PBb)) {warning("Problem in numerical appoximation of probabilities: One or more probabilities are 0. Try using larger nsim.")}

  } else if (nearorder) {

    inequalities <- as.numeric(unlist(strsplit(hypothesis, split="<")))
    inequalities <- rep(list(inequalities), J-1)
    for (j in 1:(J-1)) {inequalities[[j]][c(j+1, j)] <- inequalities[[j]][c(j, j+1)]}
    variances <- unique(unlist(inequalities))                                                              

    samp <- matrix(NA, nrow=nsim, ncol=M)

    for (i in 1:length(variances)) {samp[, variances[i]] <- SSF[variances[i]]/rchisq(nsim, dfF[variances[i]])}  
    ineqs <- matrix(0, nrow=nsim, ncol=length(inequalities))
    for (i in 1:length(inequalities)) {for (j in 1:(length(inequalities[[i]])-1)) {ineqs[, i] <- ineqs[, i]+(samp[, inequalities[[i]][j]]<samp[, inequalities[[i]][j+1]])}}
    PF <- sum(apply(ineqs, 1, max)==M-1)/nsim

    for (i in 1:length(variances)) {samp[, variances[i]] <- SSFb[variances[i]]/rchisq(nsim, dfFb[variances[i]])}
    ineqs <- matrix(0, nrow=nsim, ncol=length(inequalities))                                   
    for (i in 1:length(inequalities)) {for (j in 1:(length(inequalities[[i]])-1)) {ineqs[, i] <- ineqs[, i]+(samp[, inequalities[[i]][j]]<samp[, inequalities[[i]][j+1]])}}  
    PFb <- sum(apply(ineqs, 1, max)==M-1)/nsim

    for (i in 1:length(variances)) {samp[, variances[i]] <- SSB[variances[i]]/rchisq(nsim, dfB[variances[i]])}
    ineqs <- matrix(0, nrow=nsim, ncol=length(inequalities))                                    
    for (i in 1:length(inequalities)) {for (j in 1:(length(inequalities[[i]])-1)) {ineqs[, i] <- ineqs[, i]+(samp[, inequalities[[i]][j]]<samp[, inequalities[[i]][j+1]])}}  
    PB <- sum(apply(ineqs, 1, max)==M-1)/nsim

    PBb <- (M-1)/factorial(M)

    logmbx <- logmbx+c(log(PF/PFb), log(PB/PBb), log(PF/PBb))
    #if (0%in%c(PF, PFb, PB, PBb)) {warning("Problem in numerical appoximation of probabilities: One or more probabilities are 0. Try using larger nsim.")}

  } else if (grepl("<", hypothesis)) {
  
    constraints <- unlist(strsplit(gsub(" ", "", hypothesis), split="[0123456789]"))
    constraints <- constraints[!constraints==""]
    inequalities <- vector("character", J+length(constraints))     
    if (constraints[1]=="(") {
      inequalities[c(T, F)] <- constraints
      inequalities[c(F, T)] <- as.character(1:J)
    } else {
      inequalities[c(T, F)] <- as.character(1:J)
      inequalities[c(F, T)] <- constraints
    }
    inequalities <- paste(inequalities, sep="", collapse="")
    
    equalities <- gsub("[()]", "", inequalities)      
    equalities <- unlist(strsplit(equalities, split="[,<]"))
    equalities <- lapply(as.list(unique(equalities)), function(x) as.numeric(unlist(strsplit(x, split="="))))  

    for (i in 1:length(equalities)) {inequalities <- sub(paste(equalities[[i]], collapse="="), as.character(i), inequalities)}     
    if (grepl("[(]", inequalities)) {
      brackets <- c(unlist(gregexpr(pattern="[(]", inequalities)), unlist(gregexpr(pattern="[)]", inequalities)))
      brackets <- matrix(brackets, nrow=length(brackets)/2, ncol=2)
      commas <- unlist(apply(brackets, 1, function(x) unlist(gregexpr(pattern=",", substr(inequalities, start=x[1], stop=x[2])))+(x[1]-1)))   
      for (i in 1:length(commas)) {substring(inequalities, commas[i]) <- ";"}       
    }              
    inequalities <- unlist(strsplit(inequalities, split=","))
    inequalities <- as.list(inequalities[grepl("<", inequalities)])                                                                
    for (i in 1:length(inequalities)) {                                
      if (grepl("[(]", inequalities[[i]])) {
        combinations <- unlist(strsplit(inequalities[[i]], split="<"))
        combinations <- gsub("[()]", "", combinations)
        combinations <- strsplit(combinations, split=";")
        inequalities[[i]] <- apply(as.matrix(expand.grid(combinations)), 1, function(x) paste(x, collapse="<"))
      }
    }
    inequalities <- unlist(inequalities)
    ninequalities <- length(unlist(gregexpr(pattern="<", inequalities)))         
    inequalities <- lapply(strsplit(unlist(inequalities), split="<"), as.numeric)                                       
    variances <- unique(unlist(inequalities))

    samp <- matrix(NA, nrow=nsim, ncol=M)                                            

    for (i in 1:length(variances)) {samp[, variances[i]] <- SSF[variances[i]]/rchisq(nsim, dfF[variances[i]])}
    ineqs <- rep(0, nsim)                                      
    for (i in 1:length(inequalities)) {for (j in 1:(length(inequalities[[i]])-1)) {ineqs <- ineqs+(samp[, inequalities[[i]][j]]<samp[, inequalities[[i]][j+1]])}}
    PF <- sum(ineqs==ninequalities)/nsim

    for (i in 1:length(variances)) {samp[, variances[i]] <- SSFb[variances[i]]/rchisq(nsim, dfFb[variances[i]])}
    ineqs <- rep(0, nsim)                                      
    for (i in 1:length(inequalities)) {for (j in 1:(length(inequalities[[i]])-1)) {ineqs <- ineqs+(samp[, inequalities[[i]][j]]<samp[, inequalities[[i]][j+1]])}}
    PFb <- sum(ineqs==ninequalities)/nsim

    for (i in 1:length(variances)) {samp[, variances[i]] <- SSB[variances[i]]/rchisq(nsim, dfB[variances[i]])}
    ineqs <- rep(0, nsim)                                      
    for (i in 1:length(inequalities)) {for (j in 1:(length(inequalities[[i]])-1)) {ineqs <- ineqs+(samp[, inequalities[[i]][j]]<samp[, inequalities[[i]][j+1]])}}
    PB <- sum(ineqs==ninequalities)/nsim

    if (length(inequalities)==1) {
      PBb <- 1/factorial(length(inequalities[[1]]))
    } else {
      intersections <- sapply(inequalities, paste, collapse=",")
      intersections <- combn(intersections, 2)
      intersections <- apply(intersections, 2, function(x) intersect(as.numeric(unlist(strsplit(x[1], split=","))), as.numeric(unlist(strsplit(x[2], split=",")))))     # as.numeric() not necessary 
      if (length(intersections)==0) {
        PBb <- 1/prod(factorial(sapply(inequalities, length)))
      } else {
        for (i in 1:length(variances)) {samp[, variances[i]] <- runif(nsim)}
        ineqs <- rep(0, nsim)
        for (i in 1:length(inequalities)) {for (j in 1:(length(inequalities[[i]])-1)) {ineqs <- ineqs+(samp[, inequalities[[i]][j]]<samp[, inequalities[[i]][j+1]])}}
        PBb <- sum(ineqs==ninequalities)/nsim
      }
    }

    logmbx <- logmbx+c(log(PF/PFb), log(PB/PBb), log(PF/PBb))
    #if (0%in%c(PF, PFb, PB, PBb)) {warning("Problem in numerical appoximation of probabilities: One or more probabilities are 0. Try using larger nsim.")}

  }

  return(logmbx)

}



########## log marginal likelihoods:
log.marginal.likelihoods <- function(s2, n, b="default", hypotheses, nsim=1e5) {   
  lml <- matrix(NA, nrow=length(hypotheses), ncol=3, dimnames=list(paste("H", as.character(1:length(hypotheses)), sep=""), c("FBF", "BBF", "aFBF")))      
  for (h in 1:length(hypotheses)) {lml[h, ] <- log.marginal.likelihood(s2, n, b, hypotheses[h], nsim)}
  return(lml)
}



########## bayes factors:
bayes.factors <- function(lml, log.BF=F) {
  FBF  <- unname(sapply(lml[, 1], function(x) lml[, 1]-x, simplify=T, USE.NAMES=F))
  BBF  <- unname(sapply(lml[, 2], function(x) lml[, 2]-x, simplify=T, USE.NAMES=F))
  aFBF <- unname(sapply(lml[, 3], function(x) lml[, 3]-x, simplify=T, USE.NAMES=F))
  BF <- list(FBF, BBF, aFBF)
  names(BF) <- c("FBF", "BBF", "aFBF")  
  if (log.BF==F) {BF <- lapply(BF, exp)}  
  return(BF)
}



########## posterior probabilities:
posterior.probabilities <- function(lml, prior.probabilities=rep(1/nrow(lml), nrow(lml))) {
  if (!identical(sum(prior.probabilities), 1)) {stop("The prior probabilities do not sum to 1.")}
  PP <- matrix(NA, nrow=nrow(lml), ncol=3, dimnames=list(rownames(lml), c("FBF", "BBF", "aFBF")))
  for (i in 1:3) {PP[, i] <- prior.probabilities*sapply(lml[, i], function(x) 1/sum(exp(lml[, i]-x)*prior.probabilities), simplify=T, USE.NAMES=F)}
  return(PP)
}

  

########################################################################## TEST
s2 <- c(1, 1.5, 2.25)
n <- c(80, 90, 100)                                       

hypothesis <- "near 1<2<3"

hypothesis <- "not 1<2<3"
hypothesis <- "not 3<1<2"
hypothesis <- "not (1<2<3 or 3<2<1)"
hypothesis <- "not (3<2<1 or 1<2<3)"
hypothesis <- "not (3 < 2 < 1 or 1 < 2 < 3)"
hypothesis <- "not (1<2<3 or 3<2<1 or 1<3<2)"
hypothesis <- "not (3<2<1 or 1<3<2 or 1<2<3)"

n <- c(200, 200, 200)
s2 <- c(1, 1.5, 2.25)
hypothesis <- "3<2<1"
hypothesis <- "1=2=3"
hypotheses <- c("3<2<1", "1=2=3")
nsim <- 1e5
b <- "default"

log.marginal.likelihood(s2, n, b, hypothesis, nsim)
lml <- log.marginal.likelihoods(s2, n, b, hypotheses, nsim)
bayes.factors(lml)
posterior.probabilities(lml)