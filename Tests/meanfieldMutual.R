## simulate LV2 for critical transitions
{
  graphm = BiGraph$new(type = 'bipartite_regular', n1 = 50, k = 5, directed = T, is_adj = T)
  # variance of inital values (xinit_sd) influence the occurrence of critical transitions
  meanfieldMutual <- MeanFieldMutual$new(n1 = 50, s = 1, c = 0.004, km = 5, m = 0.5, h = 0.5, xinit_sd = 0.2, graphm = graphm)
  meanfieldMutual$rho
  meanfieldMutual$alpha_real
  s_plus_kcc <- meanfieldMutual$s + meanfieldMutual$kc * meanfieldMutual$c
  kmm <- meanfieldMutual$km * meanfieldMutual$m
  kmm / s_plus_kcc == meanfieldMutual$rho
  cs <- c(0.004, 0.006, 0.008, 0.01, 0.012)
  ss <- s_plus_kcc - meanfieldMutual$kc * cs
  
  meanfieldMutual$update(s = ss[3], c = cs[3])
  meanfieldMutual$rho
  meanfieldMutual$alpha_real
  meanfieldMutual$Delta_real
  meanfieldMutual$rmin
  r_semicircle0 <- get_r_semicircle0(meanfieldMutual$s, meanfieldMutual$c, meanfieldMutual$kc, meanfieldMutual$m, meanfieldMutual$km, meanfieldMutual$h, meanfieldMutual$n, meanfieldMutual$graphm$get_graph())
  r_semicircle0
  r_tiny_perturb = -1e-5  # a tiny perturbation of r to ensure 
  if (meanfieldMutual$alpha_real < 1) {
    r = r_semicircle0 + r_tiny_perturb
  } else {
    r = meanfieldMutual$rmin + r_tiny_perturb
  }
  r
  meanfieldMutual$sim_lv2(r = r, xinit=NULL, steps = 10000, stepwise = 1, method = 'lsoda', jactype = 'fullusr', atol = 1e-12, rtol = 1e-12, extinct_threshold = 1e-10) #lsoda lsode
  ggplot_timeseries(meanfieldMutual$simLV2$out[1:10000,], size = 0.1, xlab = 'Time', ylab = 'Abundance')
  
  # simulate pressure
  # start from 5 steps in advance to the critical r value where semicircle equal to 0
  r.stepwise = abs(meanfieldMutual$rmin / 1500)
  r.start <- r_semicircle0 + 5 * r.stepwise
  meanfieldMutual$sim_lv2_press(perturb_type = 'growthrate_all', perturb_num = 150, r.delta.mu = r.stepwise, r.delta.sd = 0, xstars.sd = 0.2, is.out = FALSE, r = r.start, r.sd = 0, steps = 10000, stepwise = 1,  method = 'lsoda', jactype = 'fullusr', atol = 1e-12, rtol = 1e-12, extinct_threshold = 1e-10)
  out <- meanfieldMutual$simLV2Press$out_press
  plot_ode_output_1(out)

  #get_stability_from_params_xstars(simLV2$params, simLV2$xstars)
  #matplot(meanfieldMutual$simLV2$out[1:100,-1], type = 'l')
  
}

## simulate LV2 press for variance of intrinsic growth rates (r), self-regulation(s), mutualistic strength(m)
{
  graphm = BiGraph$new(type = 'bipartite_regular', n1 = 50, k = 5, directed = T, is_adj = T)
  coeffs_variance = get_coeffs_variance(graphm = graphm$get_graph(), n1 = 50, km = 5, s= 1, h = 0.5, rhos = c(3), alphas = c(0.5, 2), r.sd = c(0), s.sd = c(0), c.sd = c(0), m.sd = c(0.5), rmax = 0, r.stepwise = 0.005)  # r.sd = c(0.2, 0.5, 0.8 1) alphas = c(3) rhos = c(3) s.sd = c(0.5, 0.8)
  nrow(coeffs_variance)
  xstars_variance = test_lv2_press_variance(coeffs = coeffs_variance, graphm = graphm)
  #meanfieldMutual <- MeanFieldMutual$new(n1 = 50, s = 1, c = 0.002, km = 5, m = 1, h = 0.5, graphm = graphm) # m:0.5,1 c: 0.002, 0.004, 0.008, 0.01
  #plot_ode_output_1(meanfieldMutual$simLV2Press$out_press)
  
  tmp = xstars_variance[xstars_variance$s.sd==0.5 & xstars_variance$alpha == 0.5 & xstars_variance$r <= 0, ]
  tmp2 <- tmp[,c('r', xstars_names)]
  p1 <- ggplot_timeseries(tmp2, size = 0.2, xlab = 'r', ylab = 'Abundance')
  
  rs <- seq(from = max(tmp$r), to = unique(tmp$rmin) * 1.01, length.out = 500) #by = - simLV2Press$r.delta.mu
  nstars <- sapply(rs, function(r) unlist(get_xstars(r, s=unique(tmp$s), c=unique(tmp$c), kc= unique(tmp$kc), m = unique(tmp$m), km = unique(tmp$km), h = unique(tmp$h))))
  nstars <- t(nstars)
  nstars <- cbind(rs = rs, nstars)
  nstars <- data.frame(nstars)
  p1 +
    geom_line(data = nstars, aes(x = rs, y = X1), color = 'black') +
    geom_line(data = nstars, aes(x = rs, y = X2), color = 'black', linetype = 2)
  
  
}

## simulate LV2 for Initial values
{
  meanfieldMutual <- MeanFieldMutual$new(n1 = 50, s = 1, c = 0.004, km = 5, m = 1, h = 0.5)
  meanfieldMutual$rho
  meanfieldMutual$alpha_real
  r = meanfieldMutual$rmin * 0.9
  x <- get_xstars(r = r, meanfieldMutual$s, meanfieldMutual$c, meanfieldMutual$kc, meanfieldMutual$m, meanfieldMutual$km, meanfieldMutual$h)
  mu = x$X2 # x$X1   x$X2/2   x$X2
  sd = 0.2*x$X2  # x$X1 - x$X2   x$X2/2   0.2*x$X2
  xinit <- runif2(meanfieldMutual$n, mu, sd)
  meanfieldMutual$sim_lv2(r = r, xinit, steps = 200, stepwise = 1, method = 'lsoda', jactype = 'fullusr', atol = 1e-8, rtol = 1e-8, extinct_threshold = 1e-8)
  #matplot(meanfieldMutual$simLV2$out[1:100,-1], type = 'l')
  p <- ggplot_timeseries(meanfieldMutual$simLV2$out[1:100,], size = 0.2, xlab = 'Time', ylab = 'Abundance')
  p +
    geom_point(aes(x = 0, y = mu), size = 2, shape = 19, colour = 'black') +
    geom_point(aes(x = 0, y = mu + sd), size = 2, shape = 25, colour = 'black') +
    geom_point(aes(x = 0, y = mu - sd), size = 2, shape = 24, colour = 'black')
}

## simulate SLV2
{
  meanfieldMutual <- MeanFieldMutual$new(n1 = 10, s = 1, c = 0.02, km = 4, m = 0.8, h = 0.5)
  meanfieldMutual$sim_slv2(rmax = 0, steps = 1000, stepwise = 1, sigma = 0.02, extension = 1.1, xinit_sd = 0.1)
  print_MeanFieldMutual(meanfieldMutual)
  meanfieldMutual$get_params()
  meanfieldMutual$sim_slv2(rmax = 0, steps = 2000, stepwise = 1, sigma = 0.03, extension = 1.1, xinit_sd = 0.1)
  print_MeanFieldMutual(meanfieldMutual)
  meanfieldMutual$update(h = 0.4) # decrease [h] lead to unstable
  meanfieldMutual$update(h = 0.5, m = 0.8)
  meanfieldMutual$update(m = 0.6) # decrease [m] lead to earlier critical transitions
  meanfieldMutual$update(m = 0.4)
  
  meanfieldMutual$update(m = 0.8, km = 4, c = 0.02)
  meanfieldMutual$update(m = 0.4, km = 8) # increase [km], increase gap, and
  s_plus_kcc <- meanfieldMutual$s + meanfieldMutual$kc * meanfieldMutual$c
  kmm <- meanfieldMutual$km * meanfieldMutual$m
  kmm / s_plus_kcc == meanfieldMutual$rho
  cs <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08)
  ss <- s_plus_kcc - meanfieldMutual$kc * cs
  meanfieldMutual$update(s = ss[8], c = cs[8])
  meanfieldMutual$rho
  meanfieldMutual$alpha_real
  meanfieldMutual$Delta_real
}

