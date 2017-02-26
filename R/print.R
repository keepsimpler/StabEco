#' @title print results of consistent transitions and splitting transitions
#' @param meanfiledMutual1 results of consistent transitions
#' @param meanfiledMutual2 results of splitting transitions
#' @examples 
#' load("Data/meanfieldMutual_splitting_c0006.RData")
#' meanfieldMutual1 = meanfieldMutual
#' load("Data/meanfieldMutual_splitting_c0008.RData")
#' meanfieldMutual2 = meanfieldMutual
#' print_splitting(meanfieldMutual1, meanfieldMutual2)
print_splitting <- function(meanfieldMutual1, meanfieldMutual2) {
  #
  simLV2Press = meanfieldMutual2$simLV2Press
  out = simLV2Press$out_press
  out_xstars = sapply(out, function(one) {
    one$xstars
  })
  out_xstars = t(out_xstars)
  out_length <- dim(out_xstars)[1]
  rs <- seq(from = simLV2Press$rmax, by = - simLV2Press$r.delta.mu, length.out = out_length)
  timeseries2 <- cbind(rs, out_xstars)
  # find the critical point
  splittings <- which(apply(timeseries2[,-1], 1, sd) > 1e-5)
  if (length(splittings) == 0) {
    critical_point <- max(which(apply(timeseries2[,-1], 1, mean) > 1e-10))
  } else {
    critical_point <- min(splittings)
  }
  r_critical_2 <- timeseries2[critical_point, 1]
  x_critical_2 <- mean(unlist(timeseries2[critical_point, -1]))
  #  ggplot_timeseries(timeseries, xlab = 'r', ylab = 'Abundance', size = 0.4)
  timeseries_length <- dim(timeseries2)[1]
  variable_num <- dim(timeseries2)[2]
  timeseries2 <- data.frame(timeseries2)
  colnames(timeseries2)[1] <- 'time'  # the name of the index column
  require(reshape2)
  # reshape to a long format
  timeseries2.long <- melt(data = timeseries2, id.vars = 'time', variable.name = "variable", value.name = "value")
  
  rs <- seq(from = -0.15, to = meanfieldMutual2$rmin * 1.01, length.out = 100) #by = - simLV2Press$r.delta.mu
  nstars <- sapply(rs, function(r) unlist(get_xstars(r, meanfieldMutual2$s, meanfieldMutual2$c, meanfieldMutual2$kc, meanfieldMutual2$m, meanfieldMutual2$km, meanfieldMutual2$h)))
  nstars <- t(nstars)
  nstars <- cbind(rs = rs, nstars)
  nstars <- data.frame(nstars)
  nstars[is.nan(nstars$X1),]$X1 = 0
  nstars[is.nan(nstars$X2),]$X2 = 0
  r_critical_1 = nstars[max(which(nstars$X1 > 0)), ]$rs
  x_critical_1 = nstars[max(which(nstars$X1 > 0)), ]$X1
  
    p1 <- ggplot() + #, lty = variable
    geom_line(data = timeseries2.long, aes(x = time,y = value, group = variable), size = 0.2, colour = 'blue') +
    geom_line(data = nstars, aes(x = rs, y = X1), size = 0.2, color = 'black', alpha = 0.5) +
    geom_point(aes(x = r_critical_2, y = x_critical_2), size = 2, shape = 19, colour = 'black') +
      geom_point(aes(x = r_critical_1, y = x_critical_1), size = 2, shape = 19, colour = 'black') +
      theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = ifelse(test = FALSE,"right", "none")) +
    labs(x = 'r', y = 'Abundacne')
  
  # print 
  simLV2 = meanfieldMutual2$simLV2
  out2 = simLV2$out[1:5000, ]
  #ggplot_timeseries(out, size = 0.2, xlab = 'Time', ylab = 'Abundance')
  #  ggplot_timeseries(timeseries, xlab = 'r', ylab = 'Abundance', size = 0.4)
  timeseries_length <- dim(out2)[1]
  variable_num <- dim(out2)[2]
  out2 <- data.frame(out2)
  colnames(out2)[1] <- 'time'  # the name of the index column
  require(reshape2)
  # reshape to a long format
  out2.long <- melt(data = out2, id.vars = 'time', variable.name = "variable", value.name = "value")
  
  simLV2 = meanfieldMutual1$simLV2
  out1 = simLV2$out[1:5000, ]
  #ggplot_timeseries(out, size = 0.2, xlab = 'Time', ylab = 'Abundance')
  #  ggplot_timeseries(timeseries, xlab = 'r', ylab = 'Abundance', size = 0.4)
  timeseries_length <- dim(out1)[1]
  variable_num <- dim(out1)[2]
  out1 <- data.frame(out1)
  colnames(out1)[1] <- 'time'  # the name of the index column
  require(reshape2)
  # reshape to a long format
  out1.long <- melt(data = out1, id.vars = 'time', variable.name = "variable", value.name = "value")

  p2 <- ggplot() + #, lty = variable
    geom_line(data = out2.long, aes(x = time,y = value, group = variable), size = 0.1, colour = 'blue') +
    geom_line(data = out1.long, aes(x = time,y = value, group = variable), size = 0.1, colour = 'black', alpha = 0.5) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = ifelse(test = FALSE,"right", "none")) +
    labs(x = 'Time', y = 'Abundacne')
  
  simLV2 = meanfieldMutual2$simLV2
  r = unique(simLV2$params$r)
  x1 = get_xstars(r, meanfieldMutual2$s, meanfieldMutual2$c, meanfieldMutual2$kc, meanfieldMutual2$m, meanfieldMutual2$km, meanfieldMutual2$h)$X1
  J = get_jacobian_from_params_xstars(simLV2$params, rep(x1, meanfieldMutual2$n))
  eigenvector1 = eigen(J$J)$vectors[,1]

  df_eigenvector_nstars <- data.frame(eigenvector = eigenvector1, nstars = simLV2$xstars)
  p3 <- ggplot(data = df_eigenvector_nstars, aes(x = eigenvector, y = nstars)) +
    geom_point(size = 0.8) +
    geom_smooth(method = 'auto', size = 0.8) +  # method = 'lm' 'loess' , se = FALSE
    xlab('Leading eigenvector') +
    ylab(NULL) +  # 'Abundance at equilibrium'
    scale_x_continuous(breaks=c(-0.2,0,0.2)) +
    scale_y_continuous(breaks=c(0, 0.3, 0.6)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # , axis.text= element_text(size = 6),axis.title=element_text(size=6)
  
  list(p1 = p1, p2 = p2, p3 = p3)
}
print_MeanFieldMutual_lv2_press <- function(meanfieldMutual) {
  if (! (class(meanfieldMutual)[1] == 'MeanFieldMutual' )) #&& classes[2] == 'R6'
    stop('Should be an instance of MeanFieldMutual R6 class')
  ## print the output of Stochastic simulation
  if (meanfieldMutual$state_lv2_press == TRUE) {
    simLV2Press = meanfieldMutual$simLV2Press
    out = simLV2Press$out_press
    out_xstars = sapply(out, function(one) {
      one$xstars
    })
    out_xstars = t(out_xstars)
    #matplot(out_xstars, type = 'l', lwd = 1., xlab = 'Time', ylab = 'Abundance')
    rs <- seq(from = simLV2Press$rmax, to = meanfieldMutual$rmin * 1.01, length.out = 500) #by = - simLV2Press$r.delta.mu
    nstars <- sapply(rs, function(r) unlist(get_xstars(r, meanfieldMutual$s, meanfieldMutual$c, meanfieldMutual$kc, meanfieldMutual$m, meanfieldMutual$km, meanfieldMutual$h)))
    nstars <- t(nstars)
    nstars <- cbind(rs = rs, nstars)
    nstars <- data.frame(nstars)
    out_length <- dim(out_xstars)[1]
    rs2 <- seq(from = simLV2Press$rmax, by = - simLV2Press$r.delta.mu, length.out = out_length)
    timeseries <- cbind(rs2, out_xstars)
    p1 <- ggplot_timeseries(timeseries, xlab = 'r', ylab = 'Abundance', size = 0.2)
    p1 <- p1 +
      geom_line(data = nstars, aes(x = rs, y = X1), color = 'black') +
      geom_line(data = nstars, aes(x = rs, y = X2), color = 'black', linetype = 2)
    return(p1)
  }
}

print_sde_transitions_pre <- function(meanfieldMutual, transition_type = c('consistent', 'splitting')) {
  transition_type = match.arg(transition_type)
  simSLV2 <- meanfieldMutual$simSLV2
  delta <- simSLV2$delta
  slv2.data <- simSLV2$get_out()
  slv2.data <- sapply(slv2.data, function(one) as.numeric(one))
  # any(abs(slv2.data) > 1e5) # check if the simulation is valid
  # matplot(slv2.data, type = 'l')
  timeseries <- preproc_trim_negative(slv2.data, trim_or_replace = 'replace')
  rs <- seq(from = simSLV2$r, to = simSLV2$rmin, length.out = simSLV2$steps + 1)
  timeseries <- cbind(rs, timeseries)
  timeseries <- data.frame(timeseries)
  colnames(timeseries)[1] <- 'time'  # the name of the index column
  require(reshape2)
  # reshape to a long format
  timeseries.long <- melt(data = timeseries, id.vars = 'time', variable.name = "variable", value.name = "value")
  timeseries.long$transition_type = transition_type
  timeseries.long$measure_type = ' '
  
  nstars <- sapply(rs, function(r) unlist(get_xstars(r, meanfieldMutual$s, meanfieldMutual$c, meanfieldMutual$kc, meanfieldMutual$m, meanfieldMutual$km, meanfieldMutual$h)))
  nstars <- t(nstars)
  nstars <- cbind(rs = rs, nstars)
  nstars <- data.frame(nstars)
  nstars$transition_type = transition_type
  nstars$measure_type = ' '
  
  # plot variance/covariance of the transient simulation
  slv2.data.trim <- preproc_trim_negative(slv2.data, trim_or_replace = 'trim')
  slv2.covariances <- run_stats_multivariate(slv2.data.trim, winsize = 50, type = 'cov')
  window_size <- round(dim(slv2.data.trim)[1] * 50/100) # the absolute window size
  window_num <- dim(slv2.data.trim)[1] - window_size + 1 # the number of windows
  slv2.covariances.Vc <- slv2.covariances$Vc
  slv2.covariances.Vc.pad <- c(rep(NaN, window_size - 1), slv2.covariances.Vc, rep(NaN, dim(slv2.data)[1] - dim(slv2.data.trim)[1]))
  slv2.covariances.Vs <- slv2.covariances$Vs
  slv2.covariances.Vs.pad <- c(rep(NaN, window_size - 1), slv2.covariances.Vs, rep(NaN, dim(slv2.data)[1] - dim(slv2.data.trim)[1]))
  slv2.covariances.dataframe <- data.frame(rs = rs, Vc = slv2.covariances.Vc.pad, Vs = slv2.covariances.Vs.pad)
  slv2.covariances.dataframe$transition_type = transition_type
  slv2.covariances.dataframe$measure_type = 'Vc'
  #slv2.covariances.dataframe = rbind(slv2.covariances.dataframe, tmp)
  
  slv2.syn.dataframe = slv2.covariances.dataframe
  slv2.syn.dataframe$syn = slv2.syn.dataframe$Vc / slv2.syn.dataframe$Vs
  slv2.syn.dataframe$transition_type = transition_type
  slv2.syn.dataframe$measure_type = 'Syn'
  
  slv2.covariances.Vc.deriv <-
    numeric_derivative(slv2.covariances.Vc, delta)
  slv2.covariances.Vc.deriv.pad <- c(rep(NaN, window_size - 1 + 1), slv2.covariances.Vc.deriv, rep(NaN, dim(slv2.data)[1] - dim(slv2.data.trim)[1] + 1))
  slv2.covariances.Vs.deriv <-
    numeric_derivative(slv2.covariances.Vs, delta)
  slv2.covariances.Vs.deriv.pad <- c(rep(NaN, window_size - 1 + 1), slv2.covariances.Vs.deriv, rep(NaN, dim(slv2.data)[1] - dim(slv2.data.trim)[1] + 1))
  slv2.covariances.deriv.dataframe <- data.frame(rs = rs, Vc = slv2.covariances.Vc.deriv.pad, Vs = slv2.covariances.Vs.deriv.pad)
  slv2.covariances.deriv.dataframe$transition_type = transition_type
  slv2.covariances.deriv.dataframe$measure_type = 'Deriv. of Vc'
  
  list(timeseries.long = timeseries.long, nstars = nstars, slv2.covariances.dataframe = slv2.covariances.dataframe, slv2.syn.dataframe = slv2.syn.dataframe, slv2.covariances.deriv.dataframe = slv2.covariances.deriv.dataframe)
}

#' @title print results of sde(transient stochastic dynamic) for consistent and splitting transitions
#' @param meanfieldMutual1 results of consistent transitions
#' @param meanfieldMutual2 results of splitting transitions
#' @examples 
#' load("Data/meanfieldMutual_sde_c002.RData")
#' meanfieldMutual1 = meanfieldMutual
#' load("~/Code/StabEco/Data/meanfieldMutual_sde_c008.RData")
#' meanfieldMutual2 = meanfieldMutual
#' print_sde_transitions_2(meanfieldMutual1, meanfieldMutual2)
print_sde_transitions_2 <- function(meanfieldMutual1, meanfieldMutual2) {
  out1 = print_sde_transitions_pre(meanfieldMutual1, transition_type = 'consistent')
  out2 = print_sde_transitions_pre(meanfieldMutual2, transition_type = 'splitting')
  timeseries.long.consistent = out1$timeseries.long
  timeseries.long.splitting = out2$timeseries.long
  nstars = out1$nstars  # out2$nstars
  slv2.covariances.consistent = out1$slv2.covariances.dataframe
  slv2.covariances.splitting = out2$slv2.covariances.dataframe
  slv2.covariances.deriv.consistent = out1$slv2.covariances.deriv.dataframe
  slv2.covariances.deriv.splitting = out2$slv2.covariances.deriv.dataframe
  slv2.syn.consistent = out1$slv2.syn.dataframe
  slv2.syn.splitting = out2$slv2.syn.dataframe
  
  p1 <- ggplot() + #, lty = variable
    geom_line(data = timeseries.long.consistent, aes(x = time,y = value, group = variable), size = 0.1, color = 'darkgray', alpha = 1) +
    geom_line(data = timeseries.long.splitting, aes(x = time,y = value, group = variable), size = 0.1, color = 'blue', alpha = 1) +
    geom_line(data = nstars, aes(x = rs, y = X1), color = 'black') +
    geom_line(data = nstars, aes(x = rs, y = X2), color = 'black', linetype = 2) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = 'none',
          axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    theme(panel.margin.x = unit(0., "lines"), panel.margin.y = unit(0., "lines")) + 
    labs(x = 'r', y = 'Abundances')

  max.splitting = max(slv2.covariances.splitting$Vc, na.rm = TRUE) # maximal value for splitting transition
  min.splitting = min(slv2.covariances.splitting$Vc, na.rm = TRUE) # maximal value for splitting transition
  min.consistent = min(slv2.covariances.consistent$Vc, na.rm = TRUE) # maximal value for splitting transition
  max.splitting.new = min.consistent * 0.8
  slv2.covariances.splitting$new_Vc = min.splitting + (slv2.covariances.splitting$Vc - min.splitting) * (max.splitting.new - min.splitting) /  (max.splitting - min.splitting)
  p_Vc <- ggplot() + #, lty = variable
    geom_line(data = slv2.covariances.consistent, aes(x = rs, y = Vc), colour = 'darkgray') +
    geom_line(data = slv2.covariances.splitting, aes(x = rs, y = new_Vc), colour = 'blue') +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = 'none',
          axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    theme(panel.margin.x = unit(0., "lines"), panel.margin.y = unit(0., "lines")) + 
    scale_y_continuous(name = expression(V[c]), breaks=c(max.splitting.new, min.consistent, 3, 5), labels = c(round(max.splitting,1), round(min.consistent,1), 3, 5)) +
    labs(x = 'r')
  
  slv2.covariances.deriv.splitting$new_Vc = slv2.covariances.deriv.splitting$Vc * 5
  p_Vc_deriv <- ggplot() + #, lty = variable
    geom_line(data = slv2.covariances.deriv.consistent, aes(x = rs, y = Vc), colour = 'darkgray') +
    geom_line(data = slv2.covariances.deriv.splitting, aes(x = rs, y = new_Vc), colour = 'blue') +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = 'none',
          axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    theme(panel.margin.x = unit(0., "lines"), panel.margin.y = unit(0., "lines")) + 
    labs(x = 'r', y = expression(partialdiff(V[c])/partialdiff(r)))
  
  max.splitting = max(slv2.syn.splitting$syn, na.rm = TRUE) # maximal value for splitting transition
  min.splitting = min(slv2.syn.splitting$syn, na.rm = TRUE) # maximal value for splitting transition
  min.consistent = min(slv2.syn.consistent$syn, na.rm = TRUE) # maximal value for splitting transition
  max.consistent = max(slv2.syn.consistent$syn, na.rm = TRUE) # maximal value for splitting transition
  min.consistent.new = max.splitting * 1.2
  slv2.syn.consistent$new_syn = max.consistent - (max.consistent - slv2.syn.consistent$syn) * (max.consistent - min.consistent.new) /  (max.consistent - min.consistent)
  p_syn <- ggplot() + #, lty = variable
    geom_line(data = slv2.syn.consistent, aes(x = rs, y = new_syn), colour = 'darkgray') +
    geom_line(data = slv2.syn.splitting, aes(x = rs, y = syn), colour = 'blue') +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = 'none') +
    theme(panel.margin.x = unit(0., "lines"), panel.margin.y = unit(0., "lines")) + 
    scale_y_continuous(name = expression(eta), breaks=c(0.25, max.splitting, min.consistent.new, 0.85), labels = c(0.25, round(max.splitting,2), round(min.consistent,2), 0.85)) +
    labs(x = 'r')
  
  list(p1 = p1, p_Vc = p_Vc, p_Vc_deriv = p_Vc_deriv, p_syn = p_syn)
}

#' @title print results of sde(transient stochastic dynamic) for consistent and splitting transitions
#' @param meanfieldMutual1 results of consistent transitions
#' @param meanfieldMutual2 results of splitting transitions
print_sde_transitions <- function(meanfieldMutual1, meanfieldMutual2) {
  out1 = print_sde_transitions_pre(meanfieldMutual1, transition_type = 'consistent')
  out2 = print_sde_transitions_pre(meanfieldMutual2, transition_type = 'splitting')
  timeseries.long = rbind(out1$timeseries.long, out2$timeseries.long)
  nstars = rbind(out1$nstars, out2$nstars)
  slv2.covariances.dataframe = rbind(out1$slv2.covariances.dataframe, out2$slv2.covariances.dataframe)
  slv2.syn.dataframe = rbind(out1$slv2.syn.dataframe, out2$slv2.syn.dataframe)
  max.splitting = max(slv2.syn.dataframe[slv2.syn.dataframe$transition_type=='splitting',]$syn, na.rm = TRUE) # maximal value for splitting transition
  min.consistent = min(slv2.syn.dataframe[slv2.syn.dataframe$transition_type=='consistent',]$syn, na.rm = TRUE) # maximal value for splitting transition
  #gap = min.consistent - max.splitting
  slv2.syn.dataframe$new_syn = 0
  slv2.syn.dataframe[slv2.syn.dataframe$transition_type=='splitting',]$new_syn = slv2.syn.dataframe[slv2.syn.dataframe$transition_type=='splitting',]$syn + (0.5 - max.splitting)
  slv2.syn.dataframe[slv2.syn.dataframe$transition_type=='consistent',]$new_syn = slv2.syn.dataframe[slv2.syn.dataframe$transition_type=='consistent',]$syn - (min.consistent - 0.5)
  
  
  slv2.covariances.deriv.dataframe = rbind(out1$slv2.covariances.deriv.dataframe, out2$slv2.covariances.deriv.dataframe)
  max.splitting = max(slv2.covariances.deriv.dataframe[slv2.covariances.deriv.dataframe$transition_type=='splitting',]$Vc, na.rm = TRUE) # maximal value for splitting transition
  max.consistent = max(slv2.covariances.deriv.dataframe[slv2.covariances.deriv.dataframe$transition_type=='consistent',]$Vc, na.rm = TRUE) # maximal value for splitting transition
  max.new = 4 * max.splitting
  slv2.covariances.deriv.dataframe$new_Vc <- ifelse(slv2.covariances.deriv.dataframe$transition_type == 'splitting', slv2.covariances.deriv.dataframe$Vc, ifelse(slv2.covariances.deriv.dataframe$Vc <= max.consistent, (slv2.covariances.deriv.dataframe$Vc) * (max.new) / (max.consistent), NA)) # 2 * max.splitting + (slv2.covariances.deriv.dataframe$Vc - 2 * max.splitting) * (max.new - 2 * max.splitting) / (max.consistent - 2 * max.splitting)
  
  p <- ggplot() + #, lty = variable
    geom_line(data = timeseries.long, aes(x = time,y = value, colour = variable), size = 0.1) +
    geom_line(data = nstars, aes(x = rs, y = X1), color = 'black') +
    geom_line(data = nstars, aes(x = rs, y = X2), color = 'black', linetype = 2) +
    geom_line(data = slv2.covariances.dataframe, aes(x = rs, y = Vc), colour = gg_color_hue(3)[1]) +
    geom_line(data = slv2.syn.dataframe, aes(x = rs, y = new_syn), colour = gg_color_hue(3)[2]) +
    geom_line(data = slv2.covariances.deriv.dataframe, aes(x = rs, y = new_Vc), colour = gg_color_hue(3)[3]) +
    facet_grid(measure_type ~ transition_type, scales = 'free_y') + # , space = 'free_y', heights = 2:1
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = ifelse(test = FALSE,"right", "none")) +
    theme(panel.margin.x = unit(0., "lines"), panel.margin.y = unit(0., "lines")) + 
    labs(x = 'r', y = 'Abundance')
  p
  
  gt = ggplotGrob(p)
  gt$heights[4] = unit(2, "null")
  grid.newpage()
  grid.draw(gt)
  gt
  
  gt$grobs[[4]]$children$axis$grobs[[1]]$children[[1]]$label
}

#' @title visualize the output of MeanFieldMutual object
#' @param meanfieldMutual a MeanFieldMutual object
#' @note can use a generic function \code{print}
print_MeanFieldMutual <- function(meanfieldMutual, remove.least = FALSE) {
  #classes <- class(meanfieldMutual)
  if (! (class(meanfieldMutual)[1] == 'MeanFieldMutual' )) #&& classes[2] == 'R6'
    stop('Should be an instance of MeanFieldMutual R6 class')
  ## print the output of Stochastic simulation
  if (meanfieldMutual$state_slv2 == TRUE) {
    simSLV2 <- meanfieldMutual$simSLV2
    delta <- simSLV2$delta
    slv2.data <- simSLV2$get_out()
    slv2.data <- sapply(slv2.data, function(one) as.numeric(one))
    # any(abs(slv2.data) > 1e5) # check if the simulation is valid
    # matplot(slv2.data, type = 'l')
    timeseries <- preproc_trim_negative(slv2.data, trim_or_replace = 'replace')
    if(remove.least == TRUE) {
      # remove the speceis with the least abundance
      timeseries <- timeseries[, - which(timeseries[100, ] == min(timeseries[100, ]))]
    }
    #matplot(timeseries, type = 'l')
    rs <- seq(from = simSLV2$r, to = simSLV2$rmin, length.out = simSLV2$steps + 1)
    timeseries <- cbind(rs, timeseries)
    p1 <- ggplot_timeseries(timeseries, xlab = 'r', ylab = 'Abundance', size = 0.1)
    nstars <- sapply(rs, function(r) unlist(get_xstars(r, meanfieldMutual$s, meanfieldMutual$c, meanfieldMutual$kc, meanfieldMutual$m, meanfieldMutual$km, meanfieldMutual$h)))
    nstars <- t(nstars)
    nstars <- cbind(rs = rs, nstars)
    nstars <- data.frame(nstars)
    p1 <- p1 +
      geom_line(data = nstars, aes(x = rs, y = X1), color = 'black') +
      geom_line(data = nstars, aes(x = rs, y = X2), color = 'black', linetype = 2)
  }
  
  # plot variance/covariance of the transient simulation
  slv2.data.trim <- preproc_trim_negative(slv2.data, trim_or_replace = 'trim')
  if(remove.least == TRUE) {
    # remove the speceis with the least abundance
    slv2.data.trim <- slv2.data.trim[, - which(slv2.data.trim[100,] == min(slv2.data.trim[100,]))]
  }
  slv2.covariances <- run_stats_multivariate(slv2.data.trim, winsize = 50, type = 'cov')
  window_size <- round(dim(slv2.data.trim)[1] * 50/100) # the absolute window size
  window_num <- dim(slv2.data.trim)[1] - window_size + 1 # the number of windows
  slv2.covariances.Vc <- slv2.covariances$Vc
  slv2.covariances.Vc.pad <- c(rep(NaN, window_size - 1), slv2.covariances.Vc, rep(NaN, dim(slv2.data)[1] - dim(slv2.data.trim)[1]))
  slv2.covariances.Vs <- slv2.covariances$Vs
  slv2.covariances.Vs.pad <- c(rep(NaN, window_size - 1), slv2.covariances.Vs, rep(NaN, dim(slv2.data)[1] - dim(slv2.data.trim)[1]))
  slv2.covariances.dataframe <- data.frame(rs = rs, Vc = slv2.covariances.Vc.pad, Vs = slv2.covariances.Vs.pad)
  
  p_Vc <- ggplot() +
    geom_line(data = slv2.covariances.dataframe, aes(x = rs, y = Vc), colour = gg_color_hue(3)[1]) +
    theme_bw() +
    labs(title = '', x = 'r', y = expression(V^c))
  p_syn <- ggplot() +
    geom_line(data = slv2.covariances.dataframe, aes(x = rs, y = Vc / Vs), colour = gg_color_hue(3)[3]) +
    theme_bw() +
    labs(title = '', x = 'r', y = expression(eta))
  
  slv2.covariances.Vc.deriv <-
    numeric_derivative(slv2.covariances.Vc, delta)
  slv2.covariances.Vc.deriv.pad <- c(rep(NaN, window_size - 1 + 1), slv2.covariances.Vc.deriv, rep(NaN, dim(slv2.data)[1] - dim(slv2.data.trim)[1] + 1))
  slv2.covariances.Vs.deriv <-
    numeric_derivative(slv2.covariances.Vs, delta)
  slv2.covariances.Vs.deriv.pad <- c(rep(NaN, window_size - 1 + 1), slv2.covariances.Vs.deriv, rep(NaN, dim(slv2.data)[1] - dim(slv2.data.trim)[1] + 1))
  slv2.covariances.deriv.dataframe <- data.frame(rs = rs, Vc = slv2.covariances.Vc.deriv.pad, Vs = slv2.covariances.Vs.deriv.pad)

  
  p_Vc_deriv <- ggplot() +
    geom_line(data = slv2.covariances.deriv.dataframe, aes(x = rs, y = Vc), colour = gg_color_hue(3)[1]) +
    theme_bw() +
    labs(title = '', x = 'r', y = expression(paste('Derivative of ', V^c)))

  {
    # p_Vs <- ggplot() +
    #   geom_line(data = slv2.covariances.dataframe, aes(x = rs, y = Vs), colour = gg_color_hue(3)[2]) +
    #   theme_bw() +
    #   labs(title = '', x = 'r', y = expression(V^s))
    
    #   p_Vs_derive <- ggplot() +
  #     geom_line(data = slv2.covariances.deriv.dataframe, aes(x = rs, y = Vs), colour = gg_color_hue(3)[2]) +
  #     theme_bw() +
  #     labs(title = '', x = 'r', y = expression(paste('Derivative of ', V^s)))
  #   
  #   dot_semicircle <- sapply(rs, function(r) unlist(get_dot_semicircle(r, meanfieldMutual$s, meanfieldMutual$c, meanfieldMutual$kc, meanfieldMutual$m, meanfieldMutual$km, meanfieldMutual$h, meanfieldMutual$n, meanfieldMutual$graphm$get_graph())))
  #   dot_semicircle <- t(dot_semicircle)
  #   dot_semicircle <- cbind(r = rs, dot_semicircle)
  #   dot_semicircle <- data.frame(dot_semicircle)
  #   # plot dot, semicircle and gap of the shadow Jacobian
  #   cols <- gg_color_hue(3); names(cols) <- c('color1', 'color2', 'color3');
  #   p_Jshadow_dot_semicircle <- ggplot(data = dot_semicircle, aes(x = r)) +
  #     geom_line(aes(y = Jshadow_dot, colour = 'color1')) +
  #     geom_line(aes(y = Jshadow_semicircle_real, colour = 'color2')) +
  #     geom_line(aes(y = Jshadow_dot - Jshadow_semicircle_real, colour = 'color3')) +
  #     scale_color_manual(name = '',
  #                        values = cols,
  #                        labels = c(expression(lambda[d]),
  #                                   expression(lambda[s]),
  #                                   expression(tilde(Delta))
  #                        )
  #     ) +
  #     theme_bw() +
  #     labs(title = '', x = 'r', y = '')  +
  #     theme(legend.position = c(.8, .6))
  #   
  #   p_J_dot_semicircle <- ggplot(data = dot_semicircle, aes(x = r)) +
  #     geom_line(aes(y = J_dot, colour = 'color1')) +
  #     geom_line(aes(y = J_semicircle_real, colour = 'color2')) +
  #     geom_line(aes(y = J_dot - J_semicircle_real, colour = 'color3')) +
  #     scale_color_manual(name = '',
  #                        values = cols,
  #                        labels = c(expression(lambda[d]),
  #                                   expression(lambda[s]),
  #                                   expression(tilde(Delta))
  #                        )
  #     ) +
  #     theme_bw() +
  #     labs(title = '', x = 'r', y = '')  +
  #     theme(legend.position = c(.8, .6))
}
  
  list(p1 = p1, p_Vc = p_Vc, p_syn = p_syn,
       p_Vc_deriv = p_Vc_deriv) # , p_Vs = p_Vs, p_Vs_derive = p_Vs_derive, p_Jshadow_dot_semicircle = p_Jshadow_dot_semicircle, p_J_dot_semicircle = p_J_dot_semicircle
}

#' @title ggplot a multivariate timeseries
#' @param timeseries a numeric matrix representing a multivariate timeseries, the first colnum is the index data
#' @param title, is_legend
ggplot_timeseries <- function(timeseries, title = NULL, xlab = NULL, ylab = NULL, is_legend = FALSE, size = 0.5) {
  timeseries_length <- dim(timeseries)[1]
  variable_num <- dim(timeseries)[2]
  timeseries <- data.frame(timeseries)
  colnames(timeseries)[1] <- 'time'  # the name of the index column
  require(reshape2)
  # reshape to a long format
  timeseries.long <- melt(data = timeseries, id.vars = 'time', variable.name = "variable", value.name = "value")
  p <- ggplot() + #, lty = variable
    geom_line(data = timeseries.long, aes(x = time,y = value, colour = variable), size = size) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = ifelse(test = is_legend,"right", "none")) +
    labs(title = title, x = xlab, y = ylab)
  p
}

#' @title plot the output of ode simulation
#' @param the output of ode simulation
#' \describe{
#'   \item{out}{output of one ODE simulation, including the trajectory of values of state variables}
#'   \item{nstar}{the values of state variables in equilibrium}
#'   \item{Phi}{the Jacobian matrix in equilibrium}
#'   \item{params}{parameters assigned to the model}
#'   \item{species.survived}{a vector of species that survived}
#' }
plot_ode_output_1 <- function(out) {
  out_xstars = sapply(out, function(one) {
    one$xstars
  })
  out_xstars = t(out_xstars)
  matplot(out_xstars, type = 'l', lwd = 1.8, xlab = 'Time', ylab = 'Abundance')
}

plot_ode_output_2 <- function(out, step) {
  ode.out <- out[[step]]$out
  matplot(ode.out[, -1], type = 'l', lwd = 1.8, xlab = 'Time', ylab = 'Abundances of species')
}

#' @title print output of pressure simulation
#' @param data, output of pressure simulation
#' @param yvar, 'xstars_mean', 'xstars_min', 'M_lambda1', 'M_tilde_lambda1', 'Jshadow_lambda1', 'persistence'
#' @examples 
#' ggplot_lv2_press(graphs_xstars[graphs_xstars$rho %in% c(1,2,4), ], yvar = 'persistence', ylabel = 'persistence')
#' ggplot_lv2_press(graphs_xstars, yvar = 'xstars_mean', ylabel = 'abundance')
ggplot_lv2_press <- function(data, yvar, ylabel, xlim = NULL, ylim = NULL, first_extinction = FALSE, data_first_extinction = NULL) {
  data$rho_label = paste('rho', '==', as.character(data$rho), sep = '')
  p1 <- ggplot(data, aes_string(x = 'r', y = yvar, group = 'variance_avg', color = 'variance_avg')) + # M_tilde_lambda1 J_lambda1 nstars_min nstars_mean
    geom_line(size = 0.2) +
    labs( x = 'r', y = ylabel) +
    scale_color_viridis(name = expression(H(G[m])), guide = 'none') +
    facet_grid(rho_label ~ ., scales = "free", labeller = label_parsed) +
    theme_bw()
  if (! is.null(xlim)) 
    p1 <- p1 + xlim(xlim)
  if (! is.null(ylim)) 
    p1 <- p1 + ylim(ylim)
  if (first_extinction == TRUE) {
    p1 +
      geom_point(data = data_first_extinction, aes_string(x = 'r', y = yval), color = 'black', size = 1)
  }
  p1
}

print_lv2_press <- function(graphs_xstars) {
  graphs_xstars_agg <- aggregate(cbind(xstars_min, xstars_mean, xstars_sd, xstars_max, persistence, M_lambda1, M_tilde_lambda1, Jshadow_lambda1, J_lambda1, variance_avg,  entropy_avg) ~ rho + alpha.y + r, graphs_xstars, mean)

  library(dplyr)
  #graphs_xstars_first_extinct <- graphs_xstars %>% group_by(rho, alpha.y) %>% filter(persistence == n) %>% filter(row_number() == n) 
  graphs_xstars_first_extinct <- graphs_xstars_agg %>% group_by(rho, alpha.y) %>% filter(xstars_min > 0) %>% filter(row_number() == 1) 
  graphs_xstars_last_extinct <- graphs_xstars_agg %>% group_by(rho, alpha.y) %>% filter(persistence == 0) %>% filter(row_number() == 1) 
  graphs_xstars_last_extinct <- graphs_xstars_agg %>% group_by(rho, alpha.y) %>% filter(xstars_min == 0) %>% filter(row_number() == 1) 
  

  # ggplot() + 
  #   geom_raster(data = tmp2, aes(x = alpha.y, y = r, fill = persistence), interpolate = TRUE) +
  #   scale_x_continuous(sec.axis = sec_axis(trans = ~., name = expression(gamma), labels = c('Inf', 2, 1.5, 1.33, 1.25))) +
  #   labs( x = expression(beta), y = 'r') +
  #   scale_fill_gradient(name = 'Persistence', limits=c(80, 99), low = "blue", high = "red") +
  #   theme_bw() 
  
  
  

  tmp = graphs_xstars_agg[graphs_xstars_agg$rho == 2, ]
  p_persistence <- 
    ggplot() + 
    geom_raster(data = tmp, aes(x = alpha.y, y = r, fill = persistence), interpolate = TRUE) +
    scale_x_continuous(sec.axis = sec_axis(trans = ~., name = expression(gamma), labels = c('Inf', 2, 1.5, 1.33, 1.25))) +
    labs( x = expression(beta), y = 'r') +
    scale_fill_gradient(name = expression(lambda[1](widetilde(J))), limits=c(min(tmp$persistence), max(tmp$persistence)-1), low = "blue", high = "red") + 
    theme_bw() 
  
  
}

#' @title print output of autonormous simulation
#' @examples 
#' print_lv2_auto(graphs_xstars[round(graphs_xstars$r,10) == 0 & graphs_xstars$rho %in% c(0.5, 1., 2), ], graphs_degrees)
#' print_lv2_auto(graphs_xstars[round(graphs_xstars$r,10) == 1, ], graphs_degrees)
print_lv2_auto <- function(graphs_xstars_auto) {
  stopifnot(nrow(graphs_xstars_auto) > 0)

  
  tmp4 = tmp[tmp$rho == 2, ]
  p_degrees_m_tilde <- ggplot() +
    geom_smooth(data = tmp4, aes_string(x = 'degree', y = 'm_tilde', group = hetero, color = hetero), size = 0.2, se = FALSE, method = 'loess', span = 0.5) +
    labs(x = expression(k[m]), y = expression(phi)) + # tilde(m)
    scale_color_viridis(name = expression(beta)) +
#    facet_grid(rho ~ ., scales = "free") +
    theme_bw()

  # p_xstars_m_tilde <- ggplot() +
  #   geom_point(data = tmp, aes_string(x = 'abundance', y = 'm_tilde', group = hetero, color = hetero), size = 1., alpha = 0.5) +
  #   labs(x = 'Abundance', y = expression(phi)) + # tilde(m)
  #   scale_color_viridis(name = expression(beta)) +
  #   facet_grid(rho ~ ., scales = "free") +
  #   theme_bw()
  
  # ggplot() +
  #   geom_point(data = tmp, aes_string(x = hetero, y = 'm_tilde * degree', group = 'log(degree)', color = 'log(degree)'), size = 1., alpha = 0.5) +
  #   labs(x = expression(H(G[m])), y = expression(phi*k[m])) + # tilde(m)
  #   scale_color_viridis(name = expression(k[m])) +
  #   facet_grid(rho ~ ., scales = "free") +
  #   theme_bw()
  
  p_xstars_min <- ggplot() +
    geom_line(data = tmp2, aes(x = alpha.y, y = xstars_min), colour = cols[1], size = 1.) +
    labs(x = expression(beta), y = 'ARS') +
    facet_grid(rho ~ ., scales = "free") +
    theme_bw()

  p_Jshadow_lambda1 <- ggplot() +
    geom_line(data = tmp2, aes(x = alpha.y, y = Jshadow_lambda1), colour = cols[1], size = 1.) +
    labs(x = expression(beta), y = expression(lambda[1](widetilde(J)))) +
    facet_grid(rho ~ ., scales = "free") +
    theme_bw()

  # ggplot() +
  #   geom_line(data = tmp2, aes(x = alpha.y, y = M_tilde_lambda1), colour = cols[1], size = 1.) + # J_lambda1 M_lambda1 M_tilde_lambda1
  #   labs(x = expression(beta), y = expression(lambda[1](widetilde(J)))) +
  #   facet_grid(rho ~ ., scales = "free") +
  #   theme_bw()
  
  list(p_degrees_xstars = p_degrees_xstars,  p_xstars_mean = p_xstars_mean, p_xstars_min = p_xstars_min, p_Jshadow_lambda1 = p_Jshadow_lambda1) #p_degrees_m_tilde = p_degrees_m_tilde, p_xstars_m_tilde = p_xstars_m_tilde,
}

#' @title print effect of heterogeneity on the feasiblility and resilience
#' @examples 
#' ggplot_lv2_resilience_feasibility(graphs_xstars[graphs_xstars$rho %in% c(0.5,1,2), ])
ggplot_lv2_resilience_feasibility <- function(graphs_xstars) {
  library(dplyr)
  graphs_xstars_first_extinct <- graphs_xstars %>% group_by(rho, graphs_index) %>% filter(persistence == n) %>% filter(row_number() == n()) 
  #graphs_xstars_first_extinct <- graphs_xstars %>% group_by(rho, graphs_index) %>% filter(J_lambda1 < 0) %>% filter(row_number() == n()) 
  graphs_xstars_first_extinct <- graphs_xstars %>% group_by(rho, graphs_index) %>% filter(xstars_min > 0) %>% filter(row_number() == n()) 
  graphs_xstars_last_extinct <- graphs_xstars %>% group_by(rho, graphs_index) %>% filter(persistence == 0) %>% filter(row_number() == 1) 
  
  p_first_extinct_0 <- ggplot() +
    geom_smooth(data = graphs_xstars_first_extinct, aes_string(x = 'alpha.y', y = 'r', group = 'factor(rho)', color = 'factor(rho)')) +
    theme_bw()
  p_last_extinct_0 <- ggplot() +
    geom_line(data = graphs_xstars_last_extinct, aes_string(x = 'variance_avg', y = 'r', group = 'factor(rho)', color = 'factor(rho)')) +
    theme_bw()
  
  p_xstars_min <- ggplot(graphs_xstars[graphs_xstars$J_lambda1 <= 0 & graphs_xstars$graphs_index != 1 & graphs_xstars$r > 0,], aes_string(x = 'variance_avg', y = 'xstars_min', group = 'r', color = 'r')) + # M_tilde_lambda1 J_lambda1 nstars_min nstars_mean
    geom_smooth(size = 0.2, se = FALSE, method = 'loess') +
    ylim(0, NA) +
    labs( x = expression(H(G[m])), y = 'ARS') +
    scale_color_viridis(name = expression(r), guide = 'none') +
    facet_grid(rho ~ ., scales = "free") +
    theme_bw()
  
  p_xstars_min_resilience <- ggplot(graphs_xstars[graphs_xstars$J_lambda1 <= 0 & graphs_xstars$graphs_index != 1,], aes_string(x = 'xstars_min', y = 'J_lambda1', group = hetero, color = hetero)) + 
    geom_line(size = 0.2) +
    labs( x = 'ARS', y = expression(lambda[1](J))) +
    scale_color_viridis(name = expression(r), guide = 'none') +
    facet_grid(rho ~ ., scales = "free") +
    theme_bw()
  
  p_resilience <- ggplot() +
    geom_smooth(data = graphs_xstars[graphs_xstars$J_lambda1 <= 0 & graphs_xstars$graphs_index != 1 & graphs_xstars$r > 0,], aes_string(x = 'variance_avg', y = 'J_lambda1', group = 'r', color = 'r'), size = 0.2, se = FALSE, method = 'loess') + 
    ylim(NA, 0) +
    labs( x = expression(H(G[m])), y = expression(lambda[1](J))) +
    scale_color_viridis(name = expression(r), guide = 'none') +
    facet_grid(rho ~ ., scales = "free") +
    theme_bw()
  
  graphs_xstars_first_extinct = graphs_xstars_first_extinct[graphs_xstars_first_extinct$graphs_index != 1, ]
  p_first_extinct <- ggplot() +
    geom_line(data = graphs_xstars_first_extinct, aes_string(x = 'alpha.y', y = 'r', group = 'factor(rho)', color = 'factor(rho)'), size = 0.2, color = 'grey', alpha = 0.5) +
    geom_smooth(data = graphs_xstars_first_extinct, aes_string(x = 'alpha.y', y = 'r', group = 'factor(rho)', color = 'factor(rho)'), method = 'loess', formula = x~y, se = FALSE) +
    #  geom_smooth(data = graphs_xstars_last_extinct, aes_string(x = 'variance_avg', y = 'r', group = 'factor(rho)', color = 'factor(rho)'), method = 'loess', se = FALSE) +
    labs(x = expression(H(G[m])), y = 'First Extinction') +
    scale_color_discrete(guide = 'none') +
    facet_grid(rho ~ ., scales = "free") +
    theme_bw()
  
  list(p_first_extinct_0 = p_first_extinct_0, p_last_extinct_0 = p_last_extinct_0, p_xstars_min = p_xstars_min, p_xstars_min_resilience = p_xstars_min_resilience, p_resilience = p_resilience, p_first_extinct = p_first_extinct)
}
