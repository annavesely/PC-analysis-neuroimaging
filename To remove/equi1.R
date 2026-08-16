dirname <- "/home/vesely/rep/"

source(paste0(dirname,"main_functions.R"))
source(paste0(dirname,"equicorr.R"))

m <- 10000
s <- 20
pw <- 0.95
nSim <- 1000

df <- data.frame()
power <- data.frame()

for(deltaPar in c(0.5, 1, 1.5)){
  for(rho in c(0,0.3,0.6,0.9)){
    res <- simsEqui(m=m, s=s, deltaPar=deltaPar, rho=rho, pw=pw, nSim=nSim)
    df <- rbind(df, res$df)
    power <- rbind(power, res$power)
  }
}


fname <- paste0(dirname,"res_equi_", pw*100, ".RData")
save(df, power, dim, s, repLev, file=fname)

