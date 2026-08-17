
# source("main_functions.R")

require(RNifti)

setwd("~")

mask <- RNifti::readNifti("~/neurovault_2447/mask.nii.gz")
mask_v <- as.vector(mask)

n <- 10
m <- sum(mask)

# Create matrix of t statistics, with n rows = subjects
# and m columns = voxels

copes <- matrix(0, nrow=n, ncol=m)
sids <- 1:n

for (i in seq(n)){
  fname <- paste0("~/neurovault_2447/tmap_", sids[i], ".nii.gz")
  copes[i,] <- as.vector(RNifti::readNifti(fname))[mask_v==1]
}


x <- as.vector(copes)

# Compute p-values approximating the t distribution with a standard normal
p <- 2*pnorm(abs(x), lower.tail=FALSE)

pmat0 <- matrix(p, ncol=ncol(copes))
pmat0 <- t(pmat0)


pmat <- t(apply(pmat0, 1, function(x) sort(x, na.last=TRUE)))

time_bhy <- system.time(res_bhy <- BHY(pmat, alpha=0.05))[3]

time_co <- system.time(res_co <- cofilter(pmat, tauSet=seq(from=0.01, to=1, by=0.01), alpha=0.05))[3]

time_ada <- system.time(res_ada <- adafilter(pmat0, pmat, alpha=0.05))[3]



getbrain <- function(res){
  res[res==1] <- 0
  out <- mask_v
  out[mask_v==1] <- res
  out[mask_v==0] <- NA
  out <- array(out, dim=dim(mask))
  return(out)
}


times <- list("BHY"=time_bhy, "CoFilter"=time_co, "adaFilter"=time_ada)
results <- list("BHY"=res_bhy, "CoFilter"=res_co, "adaFilter"=res_ada)
brains <- list("BHY"=getbrain(res_bhy), "CoFilter"=getbrain(res_co), "adaFilter"=getbrain(res_ada))


df_rej <- data.frame(
  gamma = 0:10,
  BHY = as.vector(table(factor(results$BHY, levels=0:10))),
  adaFilter = as.vector(table(factor(results$adaFilter, levels=0:10))),
  CoFilter = as.vector(table(factor(results$CoFilter, levels=0:10)))
)


df_power <- df_rej

df_power$BHY <- rev(cumsum(rev(df_rej$BHY)))
df_power$adaFilter <- rev(cumsum(rev(df_rej$adaFilter)))
df_power$CoFilter <- rev(cumsum(rev(df_rej$CoFilter)))

df_power <- df_power[df_power$gamma >= 2,]


save(times, results, brains, df_rej, df_power, mask, file="res_fmri.RData")









