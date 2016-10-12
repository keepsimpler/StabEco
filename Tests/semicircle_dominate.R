if (! exists('semicircle_dominate_autonomy') || ! exists('semicircle_dominate_press'))
  load(file = 'Data/semicircle_dominate.RData')
timesteps = 10000

semicircle_dominate_autonomy_one <- semicircle_dominate_autonomy[[5]]
p1 <- ggplot_timeseries(semicircle_dominate_autonomy_one$out.critical, title = '', xlab = 'Time', ylab = 'Abundance', size = 0.1)
p1
plot(semicircle_dominate_autonomy_one$J_eigenvectors[,1], semicircle_dominate_autonomy_one$out.critical[6000,-1])
df_eigenvector_nstars <- data.frame(eigenvector = semicircle_dominate_autonomy_one$J_eigenvectors[,1], nstars = semicircle_dominate_autonomy_one$out.critical[timesteps,-1])
p2 <- ggplot(data = df_eigenvector_nstars, aes(x = eigenvector, y = nstars)) +
  geom_point() +
  geom_smooth() +  # method = 'lm'
  xlab('Elements of leading eigenvector') +
  ylab('Abundance at equilibrium') +  # at equilibrium
  theme_bw()
p2
