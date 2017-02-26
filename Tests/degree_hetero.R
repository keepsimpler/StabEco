
## simulate heterogeneity
library(plyr)
library(reshape2)
library(igraph)
n1 = 50 #n2 = n1 n = n1 + n2 kc = n1 - 1 
km = 6
s = 1
h =  0.1 #c(0.05, 0.1, 0.2, 0.4) have no effect
delta = 0. # 0 0.5
semicircle_real = real_semicircle_bipartite_regular(n1 = n1, km = km)
#semicircle_real = get_m_semicircle_real(inc_to_adj(graphs[[1]]))
# graphs
#graphs_all = get_graphs_all(n1, km, ntry = 5)
#graphs_sample = sample(1:length(graphs_all), 500)
#graphs <- graphs_all[graphs_sample]
alpha_min = 0
alpha_max = 4
by = 0.04
#length.out = 100
alphas = seq(from = alpha_min, to = alpha_max, by = by)
gammas = 1/ alphas + 1 #seq(from = 2, to = 10, by = 0.1) # exponent 

#graphs_all = list()
graphs_xstars <- ddply(data.frame(idx = 1:10), .variables = c('idx'), function(i) {
  graphs <- get_graphs_sf(n1, km, alpha_min = alpha_min, alpha_max = alpha_max, by = by, ntry = 500)
  #graphs_all = c(graphs_all, graphs)
  #sapply(graphs, sum)
  graphs_hetero = get_graphs_hetero(graphs)
  graphs_hetero$gamma = gammas
  graphs_hetero$alpha = alphas
  #graphs_degrees = get_graphs_degrees(graphs)
  
  coeffs = get_coeffs(graphs = graphs, graphs.start = 1, graphs.by = 1, n1 = n1, km = km, s = s, h = h, delta = delta, semicircle_real = semicircle_real, rhos = c(0.5, 1., 2), alphas = c(3), rmax = 1, r.stepwise = 0.01)  # rhos = c(0.5, 1., 2)  #rmax = 0.5
  nrow(coeffs)
  
  # simulate press
  graphs_xstars <- test_lv2_press(coeffs, graphs) #, perturb_num = 10
  graphs_xstars <- merge(graphs_xstars, graphs_hetero, by = 'graphs_index')
  graphs_xstars
})

library(viridis)
ggplot_lv2_press(graphs_xstars[graphs_xstars$rho %in% c(0.5,1,2), ], yvar = 'persistence', ylabel = 'Persistence', xlim = c(NA, 0.5))
p_xstars_mean_2 <- ggplot_lv2_press(graphs_xstars[graphs_xstars$rho %in% c(0.5,1,2), ], yvar = 'xstars_mean', ylabel = 'Mean abundance')
ggplot_lv2_press(graphs_xstars, yvar = 'xstars_min', ylabel = 'abundance min')
ggplot_lv2_press(graphs_xstars, yvar = 'M_tilde_lambda1', ylabel = 'M_tilde', first_extinction = FALSE, ylim = c(NA,2))
ggplot_lv2_press(graphs_xstars, yvar = 'Jshadow_lambda1', ylabel = 'Jshadow', ylim=c(NA,2))

# print feasibility(first extinction) and resilience(largest eigenvalue of Jacobian)
ggplot_lv2_resilience_feasibility(graphs_xstars[graphs_xstars$rho %in% c(0.5,1,2), ])

library(dplyr)
#graphs_xstars_first_extinct <- graphs_xstars %>% group_by(rho, alpha, graphs_index) %>% filter(J_lambda1 < 0) %>% filter(row_number() == n()) 
graphs_xstars_last_extinct <- graphs_xstars %>% group_by(rho, alpha, graphs_index) %>% filter(persistence == 0) %>% filter(row_number() == 1) 

tmp <- graphs_xstars[graphs_xstars$rho == 1, ]
tmp[tmp$Jshadow_lambda1 > 0.5,]$Jshadow_lambda1 = 0.5
tmp2 <- graphs_xstars_first_extinct[graphs_xstars_first_extinct$rho == 1, ]
ggplot(tmp, aes(x = r, y = Jshadow_lambda1, group = variance_avg, color = variance_avg)) + 
  geom_line(size = 0.2) +
  geom_point(data = tmp2, aes(x = r, y = Jshadow_lambda1), color = 'black', size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = 'r', y = expression(lambda[1](widetilde(J)))) +
  ylim(NA, 0.5) +
  #  xlim(NA, 0) +
  scale_color_viridis(name = expression(H(G[m]))) +
  #  facet_grid(rho ~ ., scales = "free") +
  theme_bw()

# find some critical points
#params[[length(params) + 1]] <- list(n1 = n1, km = km, s = s, c = c, m = ms, h = hs, rho = rho, rmin=rmin, rmin_hetero = min(graphs_xstars$r), xstart = graphs_xstars[graphs_xstars$graphs_index == 1, 'xstars_mean'][1], xstart_hetero = graphs_xstars[graphs_xstars$graphs_index == graphs_num, 'xstars_mean'][1], Jshadow_zero = max(which(graphs_xstars[graphs_xstars$graphs_index == graphs_num, 'Jshadow_lambda1'] < 0)), steps_hetero = length(graphs_xstars[graphs_xstars$graphs_index == graphs_num, 'Jshadow_lambda1']))

# plot(graphs_hetero$graphs_index, graphs_hetero$variance_avg)
# plot(graphs_hetero$entropy_avg)
# plot(graphs_hetero$lambda1)
# plot(graphs_hetero$lambda2)
# plot(graphs_hetero$gap)
# ggplot(data=graphs_hetero, aes(x = entropy_avg, y = lambda1)) +
#   geom_line() +
#   geom_smooth(method = 'lm')


r = 0.5 # 0 0.5 1
graphs_xstars_auto_r05 <- ddply(data.frame(idx = 1:20), .variables = c('idx'), function(one) {
  graphs <- get_graphs_sf(n1, km, alpha_min = alpha_min, alpha_max = alpha_max, by = by, ntry = 500)
  #graphs_all = c(graphs_all, graphs)
  #sapply(graphs, sum)
  graphs_hetero = get_graphs_hetero(graphs)
  graphs_hetero$gamma = gammas
  graphs_hetero$alpha = alphas
  graphs_degrees = get_graphs_degrees(graphs)
  
  coeffs = get_coeffs(graphs = graphs, graphs.start = 1, graphs.by = 1, n1 = n1, km = km, s = s, h = h, delta = delta, semicircle_real = semicircle_real, rhos = c(0.5, 1, 2, 4), alphas = c(3), rmax = r, r.stepwise = 0.01)  # rhos = c(0.5, 1., 2)
  nrow(coeffs)
  
  # simulate auto
  graphs_xstars_auto <- test_lv2(coeffs, graphs, flag = 'auto')
  #graphs_xstars_auto <- graphs_xstars_auto[, !duplicated(colnames(graphs_xstars_auto))]
  graphs_xstars_auto <- merge(graphs_xstars_auto, graphs_hetero, by = 'graphs_index')
  graphs_xstars_auto <- merge(graphs_xstars_auto, graphs_degrees, by = 'graphs_index')
  graphs_xstars_auto

  #rs = data.frame(r = seq(from = -0.8, to = -1.5, by = - 0.01))
  #ddply(rs, .variables = c('r'), function(one) {
  #})
})

print_lv2_auto(graphs_xstars_auto_r0)
print_lv2_auto(graphs_xstars_auto_r05)
