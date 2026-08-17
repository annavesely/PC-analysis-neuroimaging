

# load("res_fmri.RData")
# source("plot_utils_spatial.R")


# --------------------------------------------
# BRAIN SLICES
# --------------------------------------------

visualizeBrain(brains$adaFilter, maxG=10, x=24, y=32, z=24)
visualizeBrain(brains$BHY, maxG=10, x=24, y=32, z=24)
visualizeBrain(brains$CoFilter, maxG=10, x=24, y=32, z=24)

# landscape 10 x 4




# --------------------------------------------
# NUMBER OF REJECTIONS BY GAMMA
# --------------------------------------------

df <- pivot_longer(
  df_power,
  cols = -c(gamma),
  names_to = "method",
  values_to = "power"
)

df$method <- factor(df$method)

ggplot(data=df, aes(x=gamma, y=power, color=method)) +
  geom_line(aes(linetype=method)) + geom_point(aes(shape=method)) +
  labs(title = "", y = expression(r(gamma)), x = expression(gamma)) +
  theme_bw(base_size = 15)

