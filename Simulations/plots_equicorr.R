
require(ggplot2)
require(tidyr)

load("res_equicorr.RData")



plotFDP <- function(df, rho=0.3){
  
  df <- df[df$rho == rho,]
  
  df$delta_label <- factor(paste0("c == ", df$deltaPar))
  df$pw_label <- factor(paste0("eta == ", df$pw))
  
  p <- ggplot(data=df, aes(x=method, y=FDP, fill=method)) +
    geom_boxplot() +
    geom_hline(yintercept=alpha, linetype="dashed") +
    labs(title = "", y = "FDP", x = "") +
    facet_grid(delta_label ~ pw_label, labeller=labeller(delta_label=label_parsed, pw_label=label_parsed)) +
    theme_bw(base_size = 15) +
    theme(legend.position = "none")
  
  print(p)
  
  out <- aggregate(df$FDP, list(df$method, df$pw, df$deltaPar), FUN=function(x) round(mean(x),4))
  names(out) <- c("method","pw","deltaPar","FDP")
  return(out)
}




plotPower <- function(power, rho=0.3){
  
  power <- power[power$rho==rho,]
  power$rho <- NULL
  
  df <- aggregate(power, by=list(power$gamma, power$pw, power$deltaPar), FUN=mean)
  df <- df[,c("pw","deltaPar","gamma","adaFilter","BHY","CoFilter")]
  
  df <- pivot_longer(
    df,
    cols = -c(pw,deltaPar,gamma),
    names_to = "method",
    values_to = "power"
  )
  
  df$method <- factor(df$method)
  df$delta_label <- factor(paste0("c == ", df$deltaPar))
  df$pw_label <- factor(paste0("eta == ", df$pw))
  
  p <- ggplot(data=df, aes(x=gamma, y=power, color=method)) +
    geom_line(aes(linetype=method)) + geom_point(aes(shape=method)) +
    labs(title = "", y = expression(beta(gamma)), x = expression(gamma)) +
    facet_grid(delta_label ~ pw_label, labeller=labeller(delta_label=label_parsed, pw_label=label_parsed)) +
    theme_bw(base_size = 15)
  
  print(p)
  
  out <- aggregate(df$power, list(df$method, df$pw, df$deltaPar), FUN=function(x) round(mean(x),4))
  names(out) <- c("method","pw","deltaPar","power")
  return(out)
}



plotFDP(df, rho=0.9)
plotPower(power, rho=0.9)


# landscape 10 x 6

