
# source("main_functions.R")



# simulation of a single study
# active: boolean vector keeping track if each of m voxels is active in the study
# rho: equicorrelation level between voxels
# n: number of observations
# alpha: significance level
# pw: power
# output: vector of p-values for each voxel

singleStudy <- function(active, rho=0, n=50, alpha=0.05, pw=0.95){
  
  m <- length(active)
  d <- power.t.test(power=pw, n=n, sig.level=alpha, type="one.sample", alternative="two.sided")$delta
  
  eps <- sqrt(1 - rho) * matrix(rnorm(m * n), ncol=m) + sqrt(rho)*matrix(rep(rnorm(n), m), ncol=m)
  
  mu <- rep(0,m)
  mu[active] <- d
  X <- matrix(rep(mu, each=n), ncol=m) + eps
  p <- apply(X,2,function(x) t.test(x, alternative="two.sided")$p.value)
  return(p)
}





# simulation of multiple studies
# m: number of voxels
# s: number of studies
# delta: boolean vector, keeping track of number of subjects for which each voxel is truly active
# rho: equicorrelation level between voxels
# n: number of observations
# alpha: significance level
# pw: power
# output: matrix of pvalues (m rows = voxels, s columns = subjects)

multiStudies <- function(m=1000, s=10, delta, rho=0, n=50, alpha=0.05, pw=0.95, seed=NULL){
  
  if(!is.null(seed)){set.seed(seed)}
  
  A <- matrix(FALSE, nrow=m, ncol=s)
  for(i in seq_along(delta)){
    sel <- sample(seq(s),delta[i])
    A[i,sel] <- TRUE
  }
  
  pmat0 <- apply(A, 2, function(x) singleStudy(x, rho, n, alpha, pw))
  return(pmat0)
}





# generate delta
# m: number of voxels
# s: number of studies
# deltaPar
# output: boolean vector, keeping track of number of subjects for which each voxel is truly active


getFixedDelta <- function(m=1000, s=10, deltaPar=1.5){
  v <- (0:s)
  pr <- deltaPar^(s-v)
  pr <- pr/sum(pr)
  ti <- round(m*pr)
  delta <- rep(v, times=round(m*pr))
  
  tmp <- length(delta)-m
  if(tmp>0){
    delta <- delta[-seq(tmp)]
  }else if(tmp<0){
    delta <- c(rep(0,-tmp),delta)
  }
  return(delta)
}






# multiple simulations in the same setting
# m: number of voxels
# s: number of studies
# delta: boolean vector, keeping track of number of subjects for which each voxel is truly active
# rho: equicorrelation level between voxels
# n: number of observations
# alpha: significance level
# pw: power
# nSim: number of repetitions
# tau: threshold for proposal 1
# tauSet: vector of threshold values for proposal 2
# output: dataframe containing param (rho), FDP, power, method

simsEqui <- function(m=1000, s=10, deltaPar=1, rho=0, n=50, alpha=0.05, pw=0.95, nSim=10,
                     tauSet=seq(from=0.01, to=1, by=0.01)){
  
  delta <- getFixedDelta(m, s, deltaPar)
  
  methods <- c("adaFilter", "BHY", "coFilter")
  
  df <- data.frame()
  power <- data.frame()
  
  
  for(i in seq(nSim)){
    if((i %% 10)==0){print(i)}
    pmat0 <- multiStudies(m, s, delta, rho, n, alpha, pw, seed=i)
    tmp <- singleAnalysis(pmat0, delta, tauSet, alpha)
    
    df <- rbind(df, tmp$df)
    
    tmp$power$sim <- i
    power <- rbind(power, tmp$power)
  }
  
  df$rho <- rho
  df$deltaPar <- deltaPar
  df$pw <- pw
  
  power$rho <- rho
  power$deltaPar <- deltaPar
  power$pw <- pw
  
  out <- list("df"=df, "power"=power)
  return(out)
}




# --------------------------------------------
# EXAMPLE
# --------------------------------------------

# res <- simsEqui(m=10000, s=20, deltaPar=1.5, rho=0.3, pw=0.95, nSim=10)


