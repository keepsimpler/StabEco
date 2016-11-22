# explore Holling type II function
# y = m x / (1 + h m x)
holling_type_2 <- function(m, h) {
  y_max = 1 / h
  x_min = - 1 / (m * h)
  
}
# generate scale-free graphs
n1 = 50
km = 4
real_semicircle_bipartite_regular <- function(n1, km) {
  tmp <- sapply(1:200, function(i) {
    biGraph <- BiGraph$new(type = 'bipartite_regular', n1 = n1, k = km)
    get_m_semicircle_real(biGraph$get_graph())
  })
  mean(tmp)
}


### print Figure 2
filename <- paste('semicircle_dominate_', '3', '.png', sep = '')
#png(filename, w = 89, h = 89 / 1.3, units = 'mm', type = 'cairo', res = 300)
width = 120 # mm
png(filename, w = width, h = width / 1.3, units = 'mm', type = 'cairo', res = 300)
#pdf(filename, width = 5, height = 5 / 1.3)
print(p3)
dev.off()

### print Figure 3
i = 1
filename <- paste('sde_transitions_', i, '.png', sep = '')
width = 120 #mm
heights = c(1.5, 2.5, 2.5, 2.5)
png(filename, w = width, h = width / heights[1], units = 'mm', type = 'cairo', res = 300)
#pdf(filename, width = 5, height = 5 / 1.3)
plot_names = c('p1', 'p_Vc', 'p_Vc_deriv', 'p_syn')
print(get(plot_names[i]))  #
dev.off()

### print Figure 4
p_colorbar <- ggplot(data = graphs_xstars, aes_string(x = 'r', y = 'xstars_mean', group = 'variance_avg', color = 'variance_avg')) +
  scale_color_viridis(name = expression(H(G[m])), guide = guide_colorbar(direction = "horizontal", position = 'bottom') ) +
  theme(legend.position = "bottom", legend.box = "horizontal") 

i = 4
filename <- paste('hetero_abundance_', i, '.png', sep = '')
width = 89 #mm
heights = c(0.7, 0.7, 0.7, 1.2)
png(filename, w = width, h = width / heights[1], units = 'mm', type = 'cairo', res = 300)
#pdf(filename, width = 5, height = 5 / 1.3)
plot_names = c('p_degrees_xstars', 'p_xstars_mean', 'p_xstars_mean_2', 'p_colorbar')
print(get(plot_names[i]))  #
dev.off()

### print Figure 5
tmp = graphs_xstars[graphs_xstars$rho %in% c(2) , c('r', 'variance_avg', 'xstars_min')] # & graphs_xstars$graphs_index != 1 & graphs_xstars$J_lambda1 <= 0  & graphs_xstars$xstars_min > 0  & graphs_xstars$r >= 0
tmp$variance_avg = round(tmp$variance_avg, 5)
#tmp[tmp$J_lambda1 > 0, ]$J_lambda1 = 0 #
tmp[abs(tmp$r) < 1e-10, ]$r = 0 # 

variance_avg2 = seq(from = min(tmp$variance_avg), to = max(tmp$variance_avg), length.out = 200)
variance_avg2 = round(variance_avg2, 2)
get_category <- function(avector, a) {
  if (length(which(avector < a)) > 0)
    cur_less = max(which(avector < a))
  else
    cur_less = 1
  if (length(which(avector > a)) > 0)
    cur_larger = min(which(avector > a))
  else
    cur_larger = length(avector)
  if (abs(a - avector[cur_less]) < abs(a - avector[cur_larger]))
    cur = cur_less
  else
    cur = cur_larger
  avector[cur]
}
get_category(variance_avg2, 150)

tmp$id = 1:nrow(tmp)
tmp2 <- ddply(tmp, .variables = c('id'), .fun = function(one) {
  variance = get_category(variance_avg2, one$variance_avg)
})
tmp$variance = tmp2$V1
tmp2 = tmp
tmp2$variance_avg = tmp$variance

tmp2 <- ddply(tmp, .variables = c('r'), .fun = function(dat) {
  data.frame(
    variance_avg = dat$variance_avg,
    xstars_min_loess = predict(loess(xstars_min ~ variance_avg, dat, span = 0.5))  #0.1
    )
  })
rs = sort(unique(tmp2$r))
#graphs_index = unique(tmp$graphs_index)
variance_avg = sort(unique(tmp2$variance_avg))
#length(rs) * length(graphs_index) == nrow(tmp)
library(reshape2)
tmp3 <- acast(data = tmp2, r ~ variance_avg, mean, value.var = 'xstars_min')
tmp3[which(tmp3 < 0)] = 0
tmp4 <- sapply(1:ncol(tmp3), function(j) {
  #last_extinct = first(which(tmp3[, j] == 0))
  first_extinct = last(which(tmp3[, j] == 0)) + 1
  first_extinct = rs[first_extinct]
  first_extinct
})
tmp5 <- predict(loess(tmp4 ~ variance_avg, span = 0.5)) #, na.action = NULL
tmp5 = c(rep(NA, length(variance_avg) - length(tmp5)), tmp5)
tmp4 = rbind(tmp4, tmp5)
# remove NaN with 0
# tmp4 <- sapply(1:ncol(tmp3), function(j) {
#   nans = which(is.na(tmp3[, j]) | is.nan(tmp3[, j]))
#   nans_zero <- nans[nans != 1:length(nans)]
#   if (length(nans_zero) > 0)
#     tmp3[nans_zero, j] = 0
#   tmp3[, j]
# })

colnames(tmp3)
rownames(tmp3)
library(RColorBrewer) #for brewer.pal()
filled.contour(rs,variance_avg,tmp3,nlevels=5,col=brewer.pal(5,"YlOrRd"))

write.table(tmp3, file = 'xstars_min_20.csv', sep = ',', row.names = FALSE, col.names = FALSE)
write.table(rs, file = 'xstars_min_rs_20.csv', sep = ',', row.names = FALSE, col.names = FALSE)
write.table(variance_avg, file = 'xstars_min_variance_20.csv', sep = ',', row.names = FALSE, col.names = FALSE)
write.table(tmp4, file = 'first_last_extinct_20.csv', sep = ',', row.names = FALSE, col.names = FALSE)

library(pracma) # use meshgrid function
tmp3 = meshgrid(rs, variance_avg)
rs = tmp3$X
variance_avg = tmp3$Y


ggplot(tmp, aes(x = variance_avg, y = r, group = 'xstars_min', color = 'xstars_min')) + # , group = 'xstars_min', color = 'xstars_min'
  geom_point() +
#  geom_tile(aes(fill = xstars_min)) + #, interpolate = TRUE
#  stat_contour(geom="polygon",aes(fill = ..level..), size = 1) +
  scale_colour_gradient(low = "magenta", high = "green")

  labs( x = expression(H(G[m])), y = 'r') +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



tmp2 = graphs_xstars_first_extinct[graphs_xstars_first_extinct$rho == 2, ]
p2 <- ggplot() +
  geom_line(data = tmp2, aes_string(x = 'variance_avg', y = 'r', group = 'factor(rho)', color = 'factor(rho)'), size = 0.2, color = 'grey', alpha = 0.7) +
  geom_smooth(data = tmp2, aes_string(x = 'variance_avg', y = 'r', group = 'factor(rho)', color = 'factor(rho)'), method = 'loess', se = FALSE) +
  #  geom_smooth(data = graphs_xstars_last_extinct, aes_string(x = 'variance_avg', y = 'r', group = 'factor(rho)', color = 'factor(rho)'), method = 'loess', se = FALSE) +
  labs(x = expression(H(G[m])), y = 'First Extinction') +
  scale_color_discrete(guide = 'none') +
#  facet_grid(rho ~ ., scales = "free") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

library(gtable)
ggplot_dual_axis(p1, p2, "y")
