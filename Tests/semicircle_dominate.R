if (! exists('semicircle_dominate')) {
  load(file = 'Data/semicircle_dominate.RData')
}
timesteps = semicircle_dominate[[1]]$coeffs$timesteps # 10000

i = 5 # (1, 0.9) (2, 1) (3, 1.2) (4, 1.5) (5, 1.8)
semicircle_dominate_one <- semicircle_dominate[[i]]
semicircle_dominate_one$description

# find the critical point according to the sd among abundances
timeseries <- semicircle_dominate_one$timeseries
splittings <- which(apply(timeseries[,-1], 1, sd) > 1e-8)
if (length(splittings) == 0) {
  critical_point <- max(which(apply(timeseries[,-1], 1, mean) > 1e-10))
} else {
  critical_point <- min(splittings)
}
r_critical <- timeseries[critical_point, 1]
x_critical <- mean(timeseries[critical_point, -1])
r_critical

p1 <- ggplot_timeseries(semicircle_dominate_one$timeseries, xlab = 'r', ylab = 'Abundance', size = 0.4)
p1 <- p1 +
  geom_point(aes(x = r_critical, y = x_critical), size = 2, shape = 19, colour = 'black')
p1

p2 <- ggplot_timeseries(semicircle_dominate_one$out.critical[1:10000,], xlab = 'Time', size = 0.1) #, ylab = 'Abundance' (1, 3000) (2, 6000)  (3, 6000) (4, 10000)
p2

filename <- paste('Output/semicircle_dominate_c', semicircle_dominate[[i]]$coeffs$c_multiplicator, '.png', sep = '')
png(filename, w = 2000, h = 1600, res = 300)
#pdf("Output/semicircle_dominate_c12.pdf", width = 10, height = 8)
require(grid)
#A viewport taking up a fraction of the plot area
vp <- viewport(width = 0.6, height = 0.59, x = 0.96, y = 0.1, just = c("right", "bottom"))
print(p1)
print(p2, vp = vp)
if (i == 4) { # c_multiplitor = 1.5
  vp <- viewport(width = 0.23, height = 0.23, x = 0.71, y = 0.2, just = c("right", "bottom"))
  print(p3, vp = vp)
}
dev.off()

plot(semicircle_dominate_one$J_eigenvectors[,1], semicircle_dominate_one$out.critical[6000,-1])
df_eigenvector_nstars <- data.frame(eigenvector = semicircle_dominate_one$J_eigenvectors[,1], nstars = semicircle_dominate_one$out.critical[timesteps,-1])
p3 <- ggplot(data = df_eigenvector_nstars, aes(x = eigenvector, y = nstars)) +
  geom_point(size = 0.8) +
  geom_smooth(method = 'auto', size = 0.5) +  # method = 'lm' 'loess'
  xlab('Leading eigenvector') +
  ylab(NULL) +  # 'Abundance at equilibrium'
  scale_x_continuous(breaks=c(-0.2,0,0.2)) +
  scale_y_continuous(breaks=c(0, 0.3, 0.6)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text= element_text(size = 6),axis.title=element_text(size=6))
# panel.background = element_blank(), axis.line = element_line(colour = "gray"))
p3

