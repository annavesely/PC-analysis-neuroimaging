
require(adaFilter)



# --------------------------------------------
# UTILS
# --------------------------------------------

# unconditional PC p-value through Fisher combination
# x: vector of p-values
# gamma: granularity
# output: PC p-value

poolPvals <- function(x, gamma){
  ns <- sum(!is.na(x))
  
  if(ns < gamma){return(1)}
  
  z <- x[gamma : ns]
  comb <- -2*sum(log(z))
  out <- pchisq(comb, df=length(z)*2, lower.tail=FALSE)
  return(out)
}




# --------------------------------------------
# ADAFILTER
# WANG, GUI, SU, SABATTI, OWEN (2020)
# --------------------------------------------

# adaFilter
# pmat0: matrix of pvalues (m rows = voxels, s columns = subjects), not sorted
# pmat: matrix of pvalues (m rows = voxels, s columns = subjects);
# within each row (voxel), the s p-values are already sorted in ascending order
# alpha: significance level
# output: vector of length m, keeping track of rejected voxels

adafilter <- function(pmat0, pmat, alpha=0.05){
  
  m <- nrow(pmat)
  s <- ncol(pmat)
  
  X <- rep(0,m) # for each feature, keeps track of the greatest gamma
  
  pv <- apply(pmat, 1, function(x) poolPvals(x, 1))
  tmp <- p.adjust(pv, method="BH")
  sel <- which(tmp <= alpha)
  X[sel] <- 1
  
  for(gamma in (2:s)){
    rej <- adaFilter(pmat0, r=gamma, type.I.err="FDR", alpha=alpha)$decision
    rej <- as.logical(rej)
    X[rej] <- gamma
  }
  
  return(X)
}




# --------------------------------------------
# BENJAMINI, HELLER, YEKUTIELI (2009)
# --------------------------------------------

# Benjamini, Heller, Yekutieli (2009)
# pmat: matrix of pvalues (m rows = voxels, s columns = subjects);
# within each row (voxel), the s p-values are already sorted in ascending order
# alpha: significance level
# output: vector of length m, keeping track of rejected voxels

BHY <- function(pmat, alpha=0.05){
  
  m <- nrow(pmat)
  s <- ncol(pmat)
  
  X <- rep(0,m) # for each feature, keeps track of the greatest gamma
  
  pv <- apply(pmat, 1, function(x) poolPvals(x, 1))
  tmp <- p.adjust(pv, method="BH")
  sel <- which(tmp <= alpha)
  m1 <- length(sel)
  thr <- alpha * m1 / m
  
  for(j in sel){
    gamma <- 0
    
    while(gamma < s){
      
      gamma <- gamma + 1
      
      pv <- poolPvals(pmat[j,], gamma)
      if(pv > thr){
        gamma <- gamma - 1
        break
      }
    }
    X[j] <- gamma
  }
  
  return(X)
}




# --------------------------------------------
# COFILTER
# DICKHAUS, HELLER, HOANG, RINOTT (2024)
# --------------------------------------------

# base for Cofilter
# pv: vector of PC p-values for different voxels
# gamma: granularity
# tau: threshold
# alpha: significance level
# output: vector of length m, keeping track of rejected voxels

baseCofilter <- function(pv, gamma, tau, alpha){
  m <- length(pv)
  rej <- rep(FALSE, m)
  sel <- which(pv < tau)
  m1 <- length(sel)
  
  if(m1 > 0){
    pvCond <- pv[sel]/tau
    tmp <- p.adjust(pvCond, method="BH")
    rej[sel] <- (tmp <= alpha) 
  }
  return(rej)
}



# Cofilter for single gamma (adaptive tau)
# pmat: matrix of pvalues (m rows = voxels, s columns = subjects);
# within each row (voxel), the s p-values are already sorted in ascending order
# gamma: granularity
# tauSet: vector of threshold values
# alpha: significance level
# output: vector of length m, keeping track of rejected voxels

cofilter_single_gamma <- function(pmat, gamma, tauSet=seq(from=0.01, to=1, by=0.01), alpha=0.05){
  m <- nrow(pmat)
  pv <- apply(pmat, 1, function(x) poolPvals(x, gamma))
  bestRej <- rep(FALSE, m)
  
  for(tau in tauSet){
    tmp <- baseCofilter(pv, gamma, tau, alpha)
    if(sum(tmp) >= sum(bestRej)){
      bestRej <- tmp
    }
  }
  return(bestRej)
}




# Cofilter (adaptive tau)
# pmat: matrix of pvalues (m rows = voxels, s columns = subjects);
# within each row (voxel), the s p-values are already sorted in ascending order
# gamma: granularity
# tauSet: vector of threshold values
# alpha: significance level
# output: vector of length m, keeping track of rejected voxels

cofilter <- function(pmat, tauSet=seq(from=0.01, to=1, by=0.01), alpha=0.05){
  
  m <- nrow(pmat)
  s <- ncol(pmat)
  
  X <- rep(0,m) # for each feature, keeps track of the greatest gamma

  for(gamma in (1:s)){
    rej <- cofilter_single_gamma(pmat, gamma, tauSet, alpha)
    X[rej] <- gamma
  }
  
  return(X)
}




# --------------------------------------------
# ANALYSIS OF SIMULATED DATA
# --------------------------------------------

# comparison between results and ground truth
# X: vector of results (for each voxel, number of detected subjects)
# delta: vector of ground truth (for each voxel, number of subjects for which the voxel is truly active)
# output:
# vector containing the FDP and the power

compareTruth <- function(X, delta, s){
  
  X[X < 2] <- 0 # we are interested only in gamma >= 2
  R <- sum(X > 0) # discoveries
  Q <- sum(X > delta) # false discoveries
  FDP <- ifelse(R==0, 0, Q/R)
  power <- sapply(2:s, function(x) sum(X >= x & delta >= x) / max(sum(delta >= x),1))
  
  out <- list("FDP"=FDP, "power"=power)
  return(out)
}




# internal function to run methods

run_method <- function(method_fun, delta, s){
  tmp <- method_fun()
  out <- compareTruth(tmp, delta, s)
  return(out)
}




# analysis of a single simulation
# pmat0: matrix of pvalues (m rows = voxels, s columns = subjects), not sorted
# delta: vector, keeping track of number of subjects for which each voxel is truly active
# tauSet: vector of threshold values
# alpha: significance level
# output:
# df: dataframe containing, for each method, the FDP
# power: dataframe containing, for each method and gamma, the power


singleAnalysis <- function(pmat0, delta, tauSet=seq(from=0.01, to=1, by=0.01), alpha=0.05){
  
  pmat <- t(apply(pmat0, 1, function(x) sort(x, na.last=TRUE)))
  s <- ncol(pmat)
  
  df <- data.frame(
    method = c("adaFilter", "BHY", "coFilter"),
    FDP = 0
  )
  
  power <- data.frame(
    gamma = 2:s
  )
  
  
  tmp <- run_method(function() adafilter(pmat0, pmat, alpha), delta, s)
  df$FDP[1] <- tmp$FDP
  power$adaFilter <- tmp$power
  
  tmp <- run_method(function() BHY(pmat, alpha), delta, s)
  df$FDP[2] <- tmp$FDP
  power$BHY <- tmp$power
  
  tmp <- run_method(function() cofilter(pmat, tauSet, alpha), delta, s)
  df$FDP[3] <- tmp$FDP
  power$coFilter <- tmp$power
  
  out <- list("df"=df, "power"=power)
  return(out)
}

