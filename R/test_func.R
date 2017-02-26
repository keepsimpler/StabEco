library(R6)
library(plyr)

### Test functions

get_coeffs_variance <- function(graphm, n1 = 50, km = 5, s = 1, h = 0.1, rhos = c(0.5, 1., 2, 3, 4), alphas = c(3), r.sd = c(0.1, 0.2), s.sd = c(0.5), c.sd = c(0.5), m.sd = c(0.5), rmax = 1, r.stepwise = 0.005) {
  semicircle_real = get_m_semicircle_real(graphm)
  stopifnot(sum(graphm) == km * n1 * 2)
  kc = n1 - 1
  alpha_min = kc * (km - semicircle_real) / ((kc + 1) * km)
  stopifnot(all(alphas > alpha_min))
  rhos = data.frame(rho = rhos) #  c(0.8, 1.0, 1.1, 1.5)
  alphas = data.frame(alpha = alphas)
  coeffs = merge(rhos, alphas)
  # should be > 1, 
  coeffs$Delta = coeffs$rho * coeffs$alpha
  coeffs$n1 = n1
  coeffs$kc = kc
  coeffs$n = 2 * n1
  coeffs$s = s
  coeffs$km = km
  coeffs$h = h
  coeffs$semi = semicircle_real
  #Delta_real <- ((km - semicircle_real) * m) / ((kc + 1) * c)
  #Delta * (kc + 1) * c * km / (km - semicircle_real) = rho * (s + kc * c) = rho * s + rho * kc * c
  coeffs$c = coeffs$rho * coeffs$s / ( coeffs$Delta * (coeffs$kc + 1) * coeffs$km / (coeffs$km - coeffs$semi) - coeffs$rho * coeffs$kc)
  coeffs$m = coeffs$Delta * (coeffs$kc + 1) * coeffs$c / (coeffs$km - coeffs$semi)
  coeffs$rmin = get_rmin(coeffs$s, coeffs$c, coeffs$kc, coeffs$m, coeffs$km, coeffs$h)
  if (nrow(coeffs[coeffs$rho < 1,]) > 0)
    coeffs[coeffs$rho < 1,]$rmin = 0
  coeffs$rmax = rmax
  coeffs$r.stepwise = r.stepwise
  coeffs$r.steps = (coeffs$rmax - coeffs$rmin) / coeffs$r.stepwise
  # merge with r.sd, s.sd, c.sd, m.sd
  coeffs = merge(coeffs, data.frame(r.sd = r.sd))
  coeffs = merge(coeffs, data.frame(s.sd = s.sd))
  coeffs = merge(coeffs, data.frame(c.sd = c.sd))
  coeffs = merge(coeffs, data.frame(m.sd = m.sd))
  coeffs$id <- 1:nrow(coeffs)
  #nrow(coeffs)
  coeffs
}

test_lv2_press_variance <- function(coeffs, graphm, extension = 10) {
  graphs_xstars <- ddply(coeffs, .variables = .(id), function(coeff) {
    print(coeff$id)
    n1 = coeff$n1
    s = coeff$s
    c = coeff$c
    kc = coeff$kc
    m = coeff$m
    km = coeff$km
    h = coeff$h
    id = coeff$id
    rho = coeff$rho
    Delta = coeff$Delta
    alpha = coeff$alpha
    rmax = coeff$rmax
    rmin = coeff$rmin
    r.stepwise = coeff$r.stepwise
    r.steps = coeff$r.steps
    r.sd = coeff$r.sd
    s.sd = coeff$s.sd
    c.sd = coeff$c.sd
    m.sd = coeff$m.sd

    meanfieldMutual <- MeanFieldMutual$new(n1 = n1, s = s, c = c, km = km, m = m, h = h, s.sd = s.sd, c.sd = c.sd, m.sd = m.sd, graphm = graphm)
    meanfieldMutual$sim_lv2_press(perturb_type = 'growthrate_all', perturb_num = r.steps * extension, r.delta.mu = r.stepwise, r.delta.sd = 0, xstars.sd = 0.1, is.out = FALSE, r = rmax, r.sd = r.sd, steps = 50, stepwise = 1,  method = 'lsoda', jactype = 'fullusr', atol = 1e-8, rtol = 1e-8, extinct_threshold = 1e-8)
    out <- meanfieldMutual$simLV2Press$out_press
    ldply(out, function(one) {
      xstars <- one$xstars
      params <- one$params
      c(s = s, c = c, kc = kc, m = m, km = km, h = h, rho = rho, alpha = alpha, Delta = Delta, rmin = rmin, rmax = rmax, r.steps = r.steps, r.sd = r.sd, s.sd = s.sd, c.sd = c.sd, m.sd = m.sd, r = params$r, get_stability_from_params_xstars(params, xstars))
    })
  })
}

#' @title Construct coefficients for simulations
#' @param 
get_coeffs <- function(graphs, graphs.start = 1, graphs.by = 1, n1 = 50, km = 5, s = 1, h = 0.1, delta = 0, semicircle_real = NULL, rhos = c(1), alphas = c(3), rmax = 1, r.stepwise = 0.005) {
  graphm = inc_to_adj(graphs[[1]])
  stopifnot(sum(graphm) == km * n1 * 2)
  if(is.null(semicircle_real)) 
    semicircle_real <- estimate_semicircle_bipartite_regular(n = 2 * n1, km = km)
  kc = n1 - 1
  alpha_min = kc * (km - semicircle_real) / ((kc + 1) * km)
  stopifnot(all(alphas > alpha_min))
  rhos = data.frame(rho = rhos) #  c(0.8, 1.0, 1.1, 1.5)
  alphas = data.frame(alpha = alphas)
  coeffs = merge(rhos, alphas)
  coeffs$Delta = coeffs$rho * coeffs$alpha
  # should be > 1, 
  #Deltas = data.frame(Delta = Deltas) #  = c(0.8, 1.0, 1.1, 1.5) seq(from = 0.5, to = 2, by = 0.4)
  coeffs$n1 = n1
  coeffs$kc = kc
  coeffs$n = 2 * n1
  coeffs$s = s
  coeffs$km = km
  coeffs$h = h
  coeffs$delta = delta
  coeffs$semi = semicircle_real
  #Delta_real <- ((km - semicircle_real) * m) / ((kc + 1) * c)
  #Delta * (kc + 1) * c * km / (km - semicircle_real) = rho * (s + kc * c) = rho * s + rho * kc * c
  coeffs$c = coeffs$rho * coeffs$s / ( coeffs$Delta * (coeffs$kc + 1) * coeffs$km / (coeffs$km - coeffs$semi) - coeffs$rho * coeffs$kc)
  coeffs$m = coeffs$Delta * (coeffs$kc + 1) * coeffs$c / (coeffs$km - coeffs$semi)
  coeffs$rmin = get_rmin(coeffs$s, coeffs$c, coeffs$kc, coeffs$m, coeffs$km, coeffs$h)
  if (nrow(coeffs[coeffs$rho < 1,]) > 0)
    coeffs[coeffs$rho < 1,]$rmin = 0
  coeffs$rmax = rmax
  coeffs$r.stepwise = r.stepwise
  coeffs$r.steps = (coeffs$rmax - coeffs$rmin) / coeffs$r.stepwise
  graphs_num = length(graphs)
  graphs_index = seq(from = graphs.start, to = graphs_num, by = graphs.by)
  coeffs = merge(coeffs, data.frame(graphs_index = graphs_index))
  coeffs$id <- 1:nrow(coeffs)
  #nrow(coeffs)
  coeffs
}

#' @title Test autonomous equilibrium of LV2 dynamics
#' @param coeffs \code{get_coeffs}
test_lv2 <- function(coeffs, graphs, xinit_sd = 0.5, flag = 'auto') {
  graphs_xstars_auto <- ddply(coeffs, .variables = .(id), function(coeff) {
    print(coeff$id)
    n1 = coeff$n1
    s = coeff$s
    c = coeff$c
    m = coeff$m
    km = coeff$km
    h = coeff$h
    id = coeff$id
    r = coeff$rmax
    rho = coeff$rho
    alpha = coeff$alpha
    Delta = coeff$Delta
    graphs_index <- coeff$graphs_index
    graphm <- graphs[[coeff$graphs_index]]
    
    meanfieldMutual <- MeanFieldMutual$new(n1 = n1, s = s, c = c, km = km, m = m, h = h, graphm = graphm, xinit_sd = xinit_sd)
    meanfieldMutual$sim_lv2(r = r, steps = 200, stepwise = 1, method = 'lsoda', jactype = 'fullusr', atol = 1e-8, rtol = 1e-8, extinct_threshold = 1e-5)
    xstars <- meanfieldMutual$simLV2$xstars
    params <- meanfieldMutual$simLV2$params
    c(n1 = n1, s = s, c = c, m = m, km = km, h = h, rho = rho, alpha = alpha, Delta = Delta, graphs_index = graphs_index, get_stability_from_params_xstars(params, xstars, flag = flag)) # r = r, 
  })
}

test_lv2_press <- function(coeffs, graphs, perturb_num = NULL,  extension = 10) {
  graphs_xstars <- ddply(coeffs, .variables = .(id), function(coeff) {
    print(coeff$id)
    n1 = coeff$n1
    s = coeff$s
    c = coeff$c
    m = coeff$m
    km = coeff$km
    h = coeff$h
    delta = coeff$delta
    id = coeff$id
    rho = coeff$rho
    Delta = coeff$Delta
    alpha = coeff$alpha
    rmax = coeff$rmax
    rmin = coeff$rmin
    r.stepwise = coeff$r.stepwise
    r.steps = coeff$r.steps
    graphs_index <- coeff$graphs_index
    graphm <- graphs[[coeff$graphs_index]]
    
    if (is.null(perturb_num))
      perturb_num = r.steps * extension
    meanfieldMutual <- MeanFieldMutual$new(n1 = n1, s = s, c = c, km = km, m = m, h = h, delta = delta, graphm = graphm)
    meanfieldMutual$sim_lv2_press(perturb_type = 'growthrate_all', perturb_num = perturb_num, r.delta.mu = r.stepwise, r.delta.sd = 0, xstars.sd = 0.5, is.out = FALSE, r = rmax, steps = 50, stepwise = 1,  method = 'lsoda', jactype = 'fullusr', atol = 1e-8, rtol = 1e-8, extinct_threshold = 1e-5)
    out <- meanfieldMutual$simLV2Press$out_press
    ldply(out, function(one) {
      xstars <- one$xstars
      params <- one$params
      c(s = s, c = c, m = m, h = h, rho = rho, alpha = alpha, Delta = Delta, rmin = rmin, rmax = rmax, r.steps = r.steps, graphs_index = graphs_index, get_stability_from_params_xstars(params, xstars, flag = 'press'))
    })
  })
  
}