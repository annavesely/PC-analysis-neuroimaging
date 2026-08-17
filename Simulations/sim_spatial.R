
require(neuRosim)

# source("main_functions.R")



# simulation of a single study
# active: boolean vector keeping track if each of m voxels is active in the study
# dim: dimensions of the 3D image
# rho: spatial correlation level
# n: number of observations
# alpha: significance level
# pw: power
# output: vector of p-values for each voxel

singleStudy <- function(active, dim=c(20,20,20), rho=0, n=50, alpha=0.05, pw=0.95){
  
  m <- prod(dim)
  d <- power.t.test(power=pw, n=n, sig.level=alpha, type="one.sample", alternative="two.sided")$delta
  
  mu <- rep(0,m)
  mu[active] <- d
  
  eps <- spatialnoise(dim=dim, sigma=1, nscan=n, method="corr", type="gaussian", rho=rho, verbose=FALSE)
  eps <- matrix(aperm(eps, c(4,1,2,3)), nrow=n, ncol=m)
  
  X <- eps + matrix(rep(mu, each=n), ncol=m)
  p <- apply(X, 2, function(x) t.test(x, alternative="two.sided")$p.value)
  
  return(p)
}




# generate spatial replicability structure
# dim: dimensions of the 3D image
# repLev: replicability levels of the 4 spheres
# output: vector of length m containing the number of subjects
# in which each voxel is truly active

getSpatialDelta <- function(dim=c(20,20,20), repLev=c(2,5,7,10)){
  
  m <- prod(dim)
  delta <- rep(0,m)
  
  # 3D grid
  grid <- array(seq_len(m), dim=dim)
  
  # 4 spheres, placed symmetrically in the volume
  q1 <- ceiling(dim/4)
  q3 <- ceiling(3*dim/4)
  
  centers <- rbind(
    c(q1[1], q1[2], q1[3]),
    c(q1[1], q1[2], q3[3]),
    c(q3[1], q3[2], q1[3]),
    c(q3[1], q3[2], q3[3])
  )
  
  radius <- floor(min(dim)*3/20)
  
  for(k in seq_along(repLev)){
    
    center <- centers[k,]
    active <- array(FALSE, dim=dim)
    
    for(i in seq_len(dim[1])){
      for(j in seq_len(dim[2])){
        for(l in seq_len(dim[3])){
          active[i,j,l] <- sum((c(i,j,l)-center)^2) <= radius^2
        }
      }
    }
    
    delta[as.vector(active)] <- repLev[k]
  }
  
  return(delta)
}





# generate active voxels for each subject from the spatial replicability structure
# delta: vector containing the number of subjects in which each voxel is active
# s: number of subjects
# output: logical matrix indicating active voxels for each subject

getSpatialActivity <- function(delta, s=10){
  
  m <- length(delta)
  A <- matrix(FALSE, nrow=m, ncol=s)
  
  for(i in seq_len(m)){
    if(delta[i] > 0){
      A[i, sample(seq(s), delta[i])] <- TRUE
    }
  }
  
  return(A)
}





# simulation of multiple studies
# dim: dimensions of the 3D image
# s: number of subjects
# delta: vector containing the true number of active subjects for each voxel
# rho: spatial correlation level
# n: number of observations
# alpha: significance level
# pw: power
# seed: random seed
# output: matrix of p-values (m rows = voxels, s columns = subjects)

multiStudiesSpatial <- function(dim=c(20,20,20), s=10, delta, rho=0, n=50, alpha=0.05, pw=0.95, seed=NULL){
  
  if(!is.null(seed)){set.seed(seed)}
  
  A <- getSpatialActivity(delta, s)
  
  pmat0 <- apply(A, 2, function(x) singleStudy(x, dim, rho, n, alpha, pw))
  
  return(pmat0)
}





# multiple simulations in the same setting
# dim: dimensions of the 3D image
# s: number of subjects
# rho: spatial correlation level
# n: number of observations
# alpha: significance level
# pw: power
# nSim: number of repetitions
# tauSet: threshold values for CoFilter
# output: dataframe containing rho, FDP, power, method

simsSpatial <- function(dim=c(20,20,20), repLev=c(2,5,7,10), s=10, rho=0, n=50, alpha=0.05, pw=0.95, nSim=10, tauSet=seq(from=0.01, to=1, by=0.01)){
  
  delta <- getSpatialDelta(dim, repLev)
  
  df <- data.frame()
  power <- data.frame()
  
  for(i in seq(nSim)){
    if((i %% 10)==0){print(i)}
    
    pmat0 <- multiStudiesSpatial(dim, s, delta, rho, n, alpha, pw, seed=i)
    tmp <- singleAnalysis(pmat0, delta, tauSet, alpha)
    
    df <- rbind(df, tmp$df)
    
    tmp$power$sim <- i
    power <- rbind(power, tmp$power)
  }
  
  df$rho <- rho
  df$pw <- pw
  
  power$rho <- rho
  power$pw <- pw
  
  out <- list("df"=df, "power"=power)
  return(out)
}




# --------------------------------------------
# EXAMPLE
# --------------------------------------------

# res <- simsSpatial(dim=c(20,20,20), repLev=c(2,5,7,10), s=10, rho=0.3, pw=0.95, nSim=10)




# --------------------------------------------
# PLOT OF SIGNAL CONFIGURATION
# --------------------------------------------

# source("plot_utils_spatial.R")

rep <- getSpatialDelta(dim=c(20,20,20), repLev=c(2,5,7,10))
rep <- array(rep, dim=c(20,20,20))
visualizeBrain(rep, x=5, y=5, z=5, maxG=10)
visualizeBrain(rep, x=15, y=15, z=15, maxG=10)

# landscape 10 x 4

