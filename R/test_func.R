library(R6)
library(plyr)

### Test functions

#' @title Construct coefficients for simulations
#' @param 
get_coeffs <- function(graphs, graphs.by = 1, n1 = 50, km = 4, s = 1, h = 0.1, rhos = c(1), Deltas = c(1), rmax = 1, r.stepwise = 0.005) {
  graphm = inc_to_adj(graphs[[1]])
  stopifnot(sum(graphm) == km * n1 * 2)
  rhos = data.frame(rho = rhos) #  c(0.8, 1.0, 1.1, 1.5)
  # should be > 1, 
  Deltas = data.frame(Delta = Deltas) #  = c(0.8, 1.0, 1.1, 1.5) seq(from = 0.5, to = 2, by = 0.4)
  coeffs = merge(rhos, Deltas)
  coeffs$n1 = n1
  coeffs$kc = n1 - 1
  coeffs$n = 2 * n1
  coeffs$s = s
  coeffs$km = km
  coeffs$h = h
  semicircle_real <- get_m_semicircle_real(graphm)
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
  graphs_index = seq(from = 1, to = graphs_num, by = graphs.by)
  coeffs = merge(coeffs, data.frame(graphs_index = graphs_index))
  coeffs$id <- 1:nrow(coeffs)
  #nrow(coeffs)
  coeffs
}

#' @title Test autonomous equilibrium of LV2 dynamics
#' @param coeffs \code{get_coeffs}
test_lv2 <- function(coeffs, graphs) {
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
    Delta = coeff$Delta
    graphs_index <- coeff$graphs_index
    graphm <- graphs[[coeff$graphs_index]]
    
    meanfieldMutual <- MeanFieldMutual$new(n1 = n1, s = s, c = c, km = km, m = m, h = h, graphm = graphm, xinit_sd = 0.5)
    meanfieldMutual$sim_lv2(r = r, steps = 200, stepwise = 1, method = 'lsoda', jactype = 'fullusr', atol = 1e-8, rtol = 1e-8, extinct_threshold = 1e-8)
    xstars <- meanfieldMutual$simLV2$xstars
    params <- meanfieldMutual$simLV2$params
    c(id = id, n1 = n1, s = s, c = c, m = m, km = km, h = h, r = r, rho = rho, Delta = Delta, graphs_index = graphs_index, get_stability_from_params_xstars(params, xstars)) # 
  })
  
}