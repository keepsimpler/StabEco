#' @title visualize the output of MeanFieldMutual object
#' @param meanfieldMutual a MeanFieldMutual object
#' @note can use a generic function \code{print}
print_MeanFieldMutual <- function(meanfieldMutual, remove.least = FALSE) {
  classes <- class(meanfieldMutual)
  if (! (classes[1] == 'MeanFieldMutual' && classes[2] == 'R6'))
    stop('Should be an instance of MeanFieldMutual R6 class')
  ## print the output of Stochastic simulation
  if (meanfieldMutual$simulated_sde == TRUE) {
    refSim <- meanfieldMutual$refSim
    delta <- refSim$delta
    slv2.data <- refSim$get_out()
    slv2.data <- sapply(slv2.data, function(one) as.numeric(one))
    # any(abs(slv2.data) > 1e5) # check if the simulation is valid
    # matplot(slv2.data, type = 'l')
    timeseries <- preproc_trim_negative(slv2.data, trim_or_replace = 'replace')
    if(remove.least == TRUE) {
      # remove the speceis with the least abundance
      timeseries <- timeseries[, - which(timeseries[100, ] == min(timeseries[100, ]))]
    }
    #matplot(timeseries, type = 'l')
    rs <- seq(from = refSim$r, to = refSim$rmin, length.out = refSim$steps + 1)
    timeseries <- cbind(rs, timeseries)
    p1 <- ggplot_timeseries(timeseries, xlab = 'r', ylab = 'Abundance', size = 0.2)
    nstars <- sapply(rs, function(r) unlist(meanfieldMutual$get_xstars(r, meanfieldMutual$s, meanfieldMutual$c, meanfieldMutual$kc, meanfieldMutual$m, meanfieldMutual$km, meanfieldMutual$h)))
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
  p_Vs <- ggplot() +
    geom_line(data = slv2.covariances.dataframe, aes(x = rs, y = Vs), colour = gg_color_hue(3)[2]) +
    theme_bw() +
    labs(title = '', x = 'r', y = expression(V^s))
  p_asyn <- ggplot() +
    geom_line(data = slv2.covariances.dataframe, aes(x = rs, y = Vs / Vc), colour = gg_color_hue(3)[3]) +
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
  p_Vs_derive <- ggplot() +
    geom_line(data = slv2.covariances.deriv.dataframe, aes(x = rs, y = Vs), colour = gg_color_hue(3)[2]) +
    theme_bw() +
    labs(title = '', x = 'r', y = expression(paste('Derivative of ', V^s)))
  
  dot_semicircle <- sapply(rs, function(r) unlist(meanfieldMutual$get_dot_semicircle(r, meanfieldMutual$s, meanfieldMutual$c, meanfieldMutual$kc, meanfieldMutual$m, meanfieldMutual$km, meanfieldMutual$h, meanfieldMutual$n, meanfieldMutual$graphm$get_graph())))
  dot_semicircle <- t(dot_semicircle)
  dot_semicircle <- cbind(r = rs, dot_semicircle)
  dot_semicircle <- data.frame(dot_semicircle)
  # plot dot, semicircle and gap of the shadow Jacobian
  cols <- gg_color_hue(3); names(cols) <- c('color1', 'color2', 'color3');
  p_Jshadow_dot_semicircle <- ggplot(data = dot_semicircle, aes(x = r)) +
    geom_line(aes(y = Jshadow_dot, colour = 'color1')) +
    geom_line(aes(y = Jshadow_semicircle_real, colour = 'color2')) +
    geom_line(aes(y = Jshadow_dot - Jshadow_semicircle_real, colour = 'color3')) +
    scale_color_manual(name = '',
                       values = cols,
                       labels = c(expression(lambda[d]),
                                  expression(lambda[s]),
                                  expression(tilde(Delta))
                       )
    ) +
    theme_bw() +
    labs(title = '', x = 'r', y = '')  +
    theme(legend.position = c(.8, .6))
  
  p_J_dot_semicircle <- ggplot(data = dot_semicircle, aes(x = r)) +
    geom_line(aes(y = J_dot, colour = 'color1')) +
    geom_line(aes(y = J_semicircle_real, colour = 'color2')) +
    geom_line(aes(y = J_dot - J_semicircle_real, colour = 'color3')) +
    scale_color_manual(name = '',
                       values = cols,
                       labels = c(expression(lambda[d]),
                                  expression(lambda[s]),
                                  expression(tilde(Delta))
                       )
    ) +
    theme_bw() +
    labs(title = '', x = 'r', y = '')  +
    theme(legend.position = c(.8, .6))
  
  list(p1 = p1, p_Vc = p_Vc, p_Vs = p_Vs, p_asyn = p_asyn,
       p_Vc_deriv = p_Vc_deriv, p_Vs_derive = p_Vs_derive, p_Jshadow_dot_semicircle = p_Jshadow_dot_semicircle, p_J_dot_semicircle = p_J_dot_semicircle)
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

