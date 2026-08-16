dirname <- "/home/vesely/rep/"

source(paste0(dirname,"main_functions.R"))
source(paste0(dirname,"spatial.R"))

dim <- c(20,20,20)
s <- 10
repLev <- c(2,5,7,10)
pw <- 0.8
nSim <- 1000

df <- data.frame()
power <- data.frame()

for(rho in c(0,0.3,0.6,0.9)){
  res <- simsSpatial(dim=dim, repLev=repLev, s=s, rho=rho, pw=pw, nSim=nSim)
  df <- rbind(df, res$df)
  power <- rbind(power, res$power)
}


fname <- paste0(dirname,"res_spatial_", pw*100, ".RData")
save(df, power, dim, s, repLev, file=fname)

