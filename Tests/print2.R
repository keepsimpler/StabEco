## print figures from output of press experiments
## input: graphs_xstars
## output: figures
{
  library(viridis)
  unique(graphs_xstars$idx)
  unique(graphs_xstars$rho)
  
  mean_fun <- function(persistences) {
    if (length(persistences) < idx_num)
      persistences = c(persistences, rep(0, idx_num - length(persistences)))
    mean(persistences)
  }
  idx_num = length(unique(graphs_xstars$idx))
  
  ## Figure: the mean abundance
  tmp <- aggregate(cbind(xstars_mean) ~ rho + alpha.y + r, graphs_xstars, mean_fun) #[graphs_xstars_rho4$persistence >= 0, ]
  tmp$rho_label = paste('rho', '==', as.character(tmp$rho), sep = '')
  tmp2 <- aggregate(cbind(xstars_mean) ~ rho + alpha.y + r, graphs_xstars, var) 
  tmp2[is.na(tmp2$xstars_mean),]$xstars_mean = 0
  tmp3 = tmp2$sd / tmp$xstars_mean
  quantile(tmp3, probs = seq(0,1,0.95), na.rm = TRUE)
  
  press_xstars_mean <- ggplot(tmp, aes_string(x = 'r', y = 'xstars_mean', group = 'alpha.y', color = 'alpha.y')) + 
    geom_line(size = 0.2) +
    labs( x = 'r', y = 'Mean Abundance') +
    scale_color_viridis(name = expression(beta), guide = 'none') +
    facet_grid(rho_label ~ ., scales = "free", labeller = label_parsed) +
    theme_bw()
  
  ## Figure: persistence
  tmp <- aggregate(cbind(persistence) ~ rho + alpha.y + r, graphs_xstars, mean_fun) #[graphs_xstars_rho4$persistence >= 0, ]

  tmp2 = tmp[tmp$persistence >= 60 , ] #& tmp$r <= 0
  library(data.table)
  tmp2 = data.table(tmp2)
  # add a column [r_min] equal to the minimal value of [r]
  tmp2[, r_min := min(r), by = list(rho, alpha.y)]
  # assign the minimal persistence for the minimal value of [r]
  tmp2[r == r_min, persistence := 60,]
  
  tmp2$rho_label = paste('rho', '==', as.character(tmp2$rho), sep = '')
  
  press_persistence_1 <- ggplot(tmp2, aes_string(x = 'r', y = 'persistence', group = 'alpha.y', color = 'alpha.y')) + 
    geom_line(size = 0.2) +
    labs( x = 'r', y = 'Persistence (%)') +
#    scale_y_continuous(labels=percent) +
    scale_color_viridis(name = expression(beta), guide = 'colorbar') + #none
    ylim(c(60, 100)) +
    facet_grid(rho_label ~ ., scales = "free", labeller = label_parsed) +
    theme_bw()

  tmp3 = tmp2[(tmp2$rho == 4 & tmp2$r < 0) | tmp2$rho %in% c(0.5, 1,2), ]
  press_persistence_2 <- ggplot(tmp3, aes_string(x = 'r', y = 'persistence', group = 'alpha.y', color = 'alpha.y')) + 
    geom_line(size = 0.2) +
    labs( x = 'r', y = 'Persistence') +
    scale_color_viridis(name = expression(beta), guide = 'none') +
    ylim(c(60, 100)) +
    facet_wrap( ~ rho_label, ncol = 1, scales = "free", labeller = label_parsed) +
    theme_bw()

  ## Figure: minimal abundances (xstars_min)
  tmp <- aggregate(cbind(xstars_min, M_lambda1, M_tilde_lambda1,Jshadow_lambda1 ) ~ rho + alpha.y + r, graphs_xstars, mean) # mean mean_fun
  tmp$rho_label = paste('rho', '==', as.character(tmp$rho), sep = '')
  tmp2 <- aggregate(cbind(xstars_min) ~ rho + alpha.y + r, graphs_xstars, var) # mean mean_fun
  tmp2[is.na(tmp2$xstars_min),]$xstars_min = 0
  tmp2$sd = sqrt(tmp2$xstars_min)
  tmp3 = tmp2$sd / tmp$xstars_min
  tmp3[which(is.na(tmp3))] = 0
  quantile(tmp3[tmp2$rho == 4], probs = seq(0,1,0.8))
  
  tmp2 = tmp[tmp$rho == 1, ]  # 0.5 1 2 4

  require(dplyr)
  #tmp2 <- graphs_xstars %>% group_by(rho, alpha.y) %>% filter(persistence == n) %>% filter(row_number() == n) 
  tmp3 <- tmp2 %>% group_by(rho, alpha.y) %>% filter(xstars_min > 0) %>% filter(row_number() == 1) 
  #tmp3 <- graphs_xstars %>% group_by(rho, alpha.y) %>% filter(persistence == 0) %>% filter(row_number() == 1) 
  tmp4 <- tmp2 %>% group_by(rho, alpha.y) %>% filter(xstars_min == 0) %>% filter(row_number() == 1) 
  
  p_xstars_min <- ggplot() + 
    geom_raster(data = tmp2, aes(x = alpha.y, y = r, fill = xstars_min), interpolate = TRUE) +
#    geom_contour(data = tmp, aes(x = alpha.y, y = r, z = xstars_min), size = 0.2, color = 'white',linemitre = 1) +
    geom_line(data = tmp3, aes(x = alpha.y, y = r), size = 0.2, color = 'darkgrey', alpha = 0.5) +
    geom_smooth(data = tmp3, aes(x = alpha.y, y = r), size = 1, method = 'loess', span = 0.5, color = 'green', se = FALSE) +
    geom_line(data = tmp4, aes(x = alpha.y, y = r), size = 1.5, color = 'black') +
    scale_x_continuous(sec.axis = sec_axis(trans = ~., name = expression(gamma), labels = c('Inf', 2, 1.5, 1.33, 1.25))) +
    labs( x = expression(beta), y = 'r') +
    scale_fill_gradient(name = 'ARS', limits=c(0, max(tmp2$xstars_min)), low = "blue", high = "red") +
#    facet_grid(rho_label ~ ., scales = "free", labeller = label_parsed) +
    theme_bw() 
  
  p_M_lambda1 <- ggplot() + 
    geom_raster(data = tmp2, aes(x = alpha.y, y = r, fill = M_lambda1), interpolate = TRUE) +
    scale_x_continuous(sec.axis = sec_axis(trans = ~., name = expression(gamma), labels = c('Inf', 2, 1.5, 1.33, 1.25))) +
    labs( x = expression(beta), y = 'r') +
    scale_fill_gradient(name = expression(lambda[1](G[m])), limits=c(min(tmp2$M_lambda1), max(tmp2$M_lambda1)), low = "blue", high = "red") +
    theme_bw() 
  
  p_M_tilde_lambda1 <- ggplot() + 
    geom_raster(data = tmp2, aes(x = alpha.y, y = r, fill = M_tilde_lambda1), interpolate = TRUE) +
    scale_x_continuous(sec.axis = sec_axis(trans = ~., name = expression(gamma), labels = c('Inf', 2, 1.5, 1.33, 1.25))) +
    labs( x = expression(beta), y = 'r') +
    scale_fill_gradient(name = expression(lambda[1](widetilde(M))), limits=c(min(tmp2$M_tilde_lambda1), max(tmp2$M_tilde_lambda1)), low = "blue", high = "red") + # max(tmp2$M_tilde_lambda1) 1
    theme_bw() 
  
  p_Jshadow_lambda1 <- ggplot() + 
    geom_raster(data = tmp2, aes(x = alpha.y, y = r, fill = Jshadow_lambda1), interpolate = TRUE) +
    scale_x_continuous(sec.axis = sec_axis(trans = ~., name = expression(gamma), labels = c('Inf', 2, 1.5, 1.33, 1.25))) +
    labs( x = expression(beta), y = 'r') +
    scale_fill_gradient(name = expression(lambda[1](widetilde(J))), limits=c(min(tmp2$Jshadow_lambda1), max(tmp2$Jshadow_lambda1)), low = "blue", high = "red") + 
    theme_bw() 

  tmp = graphs_xstars[graphs_xstars$J_lambda1<= 0, ]
  tmp2 <- aggregate(cbind(xstars_min, J_lambda1 ) ~ rho + alpha.y + r, tmp, mean) # mean mean_fun
  tmp2$rho_label = paste('rho', '==', as.character(tmp2$rho), sep = '')
  
  p_xstars_min_resilience <- 
    ggplot(tmp2, aes(x = xstars_min, y = J_lambda1, group = alpha.y, color = alpha.y)) + 
    geom_line(size = 0.2) +
    labs( x = 'ARS', y = expression(lambda[1](J))) +
    scale_color_viridis(name = expression(r), guide = 'none') +
    facet_grid(rho_label ~ ., scales = "free", labeller = label_parsed) +
    theme_bw()
  
  
}

## print figures from output of autonomous experiments
## input: graphs_xstars_auto
## output: figures
{
  print(paste('n1 = ', n1, ', km = ', km, sep = ''))
  print(paste('r = ', unique(graphs_xstars_auto$r), ', km = ', unique(graphs_xstars_auto$km)))
  
  hetero <- 'alpha.y' # 'entropy_avg' variance_avg
  m_tilde_names <- paste('m_tilde', 1:(2*n1), sep = '')
  xstars_names <- paste('xstars', 1:(2*n1), sep = '')
  degrees_names <- paste('V', 1:(2*n1), sep = '')
  require(reshape2)
  tmp <- melt(data = graphs_xstars_auto, id.vars = c(hetero, 'r', 'rho'), measure.vars = xstars_names, variable.name = 'node', value.name = 'abundance')
  tmp2 <- melt(data = graphs_xstars_auto, id.vars = c(hetero, 'r', 'rho'), measure.vars = m_tilde_names, variable.name = 'node', value.name = 'm_tilde')
  tmp3 <- melt(data = graphs_xstars_auto, id.vars = c(hetero, 'r', 'rho'), measure.vars = degrees_names, variable.name = 'node', value.name = 'degree')
  tmp <- cbind(tmp, m_tilde = tmp2$m_tilde, degree = tmp3$degree)

  tmp2 <- aggregate(cbind(abundance) ~ degree + alpha.y + rho, tmp, mean)
  tmp2$rho_label <- paste('rho', '==', as.character(tmp2$rho), sep = '')   
  tmp3 <- aggregate(cbind(abundance) ~ degree + alpha.y + rho, tmp, var)
  tmp3[is.na(tmp3$abundance),]$abundance = 0
  tmp4 = tmp3$abundance / tmp2$abundance
  quantile(tmp4, probs = seq(0,1,0.95))
  
  require(viridis)
  p_degrees_xstars <- ggplot() +
    geom_line(data = tmp2, aes_string(x = 'degree', y = 'abundance', group = hetero, color = hetero), size = 0.2, alpha = 0.5) +
    #geom_smooth(data = tmp, aes_string(x = 'degree', y = 'abundance', group = hetero, color = hetero), method = 'loess', se = FALSE, size = 0.5) +
    labs(x = expression(k[m]), y = 'Abundances of individual species') +
    scale_color_viridis(name = expression(beta), guide = 'none' ) +
    facet_grid(rho_label ~ ., scales = "free", labeller = label_parsed) +
    theme_bw()
  
  tmp2 <- aggregate(cbind(xstars_mean) ~ rho + alpha.y + gamma, graphs_xstars_auto, mean)
  #rho = 0.5
  #tmp3 = tmp2[tmp2$rho == rho, c('alpha.y', 'xstars_mean')]
  #write.table(tmp3, file = paste('means_', rho, '.csv', sep = ''), sep = ',', row.names = FALSE, col.names = TRUE)
  
  tmp3 <- aggregate(cbind(xstars_mean) ~ rho + alpha.y + gamma, graphs_xstars_auto, var) # error
  tmp3$sd = sqrt(tmp3$xstars_mean)
  tmp4 = tmp3$sd / tmp2$xstars_mean
  
  tmp2$rho_label <- paste('rho', '==', as.character(tmp2$rho), sep = '')   
  cols <- gg_color_hue(3); names(cols) <- c('color1', 'color2', 'color3');
  p_xstars_mean <- ggplot(data = tmp2, aes(x = alpha.y)) +
    #  geom_line(aes(y = xstars_max, colour = 'color1'), size = 1.) +
    geom_line(aes(y = xstars_mean), colour = cols[1], size = 1.) +
    #  geom_line(aes(y = xstars_min, colour = 'color3'), size = 1.) +
    labs(x = expression(beta), y = 'Mean abundance') +
    facet_grid(rho_label ~ ., scales = "free", labeller = label_parsed) +
    theme_bw()
  
  rho = 4 # 0.5 1 2 4
  tmp2 = tmp[tmp$rho == rho, ]
  
  tmp3 <- aggregate(cbind(abundance) ~ degree + alpha.y + r, tmp2, mean)
  tmp4 <- aggregate(cbind(abundance) ~ degree + alpha.y + r, tmp2, length) #length sum
  tmp3$count = tmp4$abundance
  #tmp5 <- aggregate(cbind(abundance) ~ alpha.y + r, tmp2, sum) # length
  write.table(tmp3, file = paste('means_counts_', rho, '.csv', sep = ''), sep = ',', row.names = FALSE, col.names = TRUE)

  ggplot() + 
    geom_raster(data = tmp3, aes(x = alpha.y, y = degree, fill = abundance), interpolate = TRUE) +
    scale_x_continuous(sec.axis = sec_axis(trans = ~., name = expression(gamma), labels = c('Inf', 2, 1.5, 1.33, 1.25))) +
    labs( x = expression(beta), y = expression(k[m])) +
    scale_fill_gradient(name = 'Abundance', limits=c(min(tmp3$abundance), max(tmp3$abundance)), low = "blue", high = "red") +
    theme_bw() 
  
  tmp3 <- acast(data = tmp2, alpha.y ~ degree, mean, value.var = 'abundance')
  library(zoo)
  tmp7 <- sapply(1:dim(tmp3)[1], function(row) {
    tmp4 = tmp3[row, ]
    tmp5 = zoo(tmp4)
    tmp6 = na.approx(tmp5)
    c(as.numeric(tmp6), rep(NaN, dim(tmp3)[2] - length(tmp6)))
  })
  write.table(tmp7, file = paste('xstars_', rho, '.csv', sep = ''), sep = ',', row.names = FALSE, col.names = FALSE)
  

  
  
}