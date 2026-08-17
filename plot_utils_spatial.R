
require(tidyr)
require(dplyr)
require(ggplot2)
require(ggpubr)
require(paletteer)
require(pals)


# visualize a 2D slice of a 3D brain map
# dataMap: 3D array
# slice: axis along which the slice is taken ("x", "y", or "z")
# x, y, z: coordinates
# maxG: maximum number of subjects, used when visualizing replicability levels
# output: plot

visualizeBrain0 <- function(dataMap, slice, x, y, z, maxG=NULL){
  
  if(slice=="x"){
    d <- dataMap[x,,]
    otherInd <- c(y,z)
    sliceInd <- x
  }else if(slice=="y"){
    d <- dataMap[,y,]
    otherInd <- c(x,z)
    sliceInd <- y
  }else{
    d <- dataMap[,,z]
    otherInd <- c(x,y)
    sliceInd <- z
  }
  
  remAxes <- setdiff(c("x","y","z"), slice)
  
  df <- expand.grid(a=seq(dim(d)[1]), b=seq(dim(d)[2]))
  df$subjects <- as.vector(d)
  
  if(!is.null(maxG)){
    df$subjects <- factor(df$subjects, levels=c(0:maxG))
    ncolor <- maxG+1
    legTitle <- "subjects"
  }else{
    ncolor <- 256
    legTitle <- "p-value"
  }
  
  linecol <- "darkorange2"
  
  p <- ggplot(df, aes(x=a, y=b, fill=subjects)) +
    geom_tile() +
    labs(title=paste0(slice, " = ", sliceInd), fill=legTitle, x=remAxes[1], y=remAxes[2]) +
    geom_vline(xintercept=otherInd[1], col=linecol) +
    geom_hline(yintercept=otherInd[2], col=linecol) +
    theme_minimal(base_size = 15) +
    theme(plot.title = element_text(hjust=0.5))
  
  cols <- paletteer::paletteer_c("pals::kovesi.linear_blue_5_95_c73", n=ncolor, direction=1)
  
  if(!is.null(maxG)){
    cols <- c("black", cols[-1])
    p <- p + scale_fill_manual(values=cols, na.value="white", breaks=(0:maxG))
  }else{
    p <- p + scale_fill_gradient(low=cols[200], high=cols[1], na.value="white")
  }
  
  p
}




# visualize 3 orthogonal slices of a 3D brain map
# dataMap: 3D array
# x, y, z: coordinates
# maxG: maximum number of subjects, used when visualizing replicability levels
# output: plot

visualizeBrain <- function(dataMap, x=32, y=32, z=16, maxG=NULL){
  
  p1 <- visualizeBrain0(dataMap, slice="x", x=x, y=y, z=z, maxG=maxG)
  p2 <- visualizeBrain0(dataMap, slice="y", x=x, y=y, z=z, maxG=maxG)
  p3 <- visualizeBrain0(dataMap, slice="z", x=x, y=y, z=z, maxG=maxG)
  
  ggarrange(p1, p2, p3, nrow=1, common.legend = TRUE, legend="right")
}



