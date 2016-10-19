library(R6)

#' @title A Mean Field Model for Mutualistic ecological communities
#' @field n1 species number of group 1
#' @field n2 species number of group 2
#' @examples 
#' meanfieldMutual <- MeanFieldMutual$new(n1 = 20, s = 1, c = 0.01, km = 5, m = 1, h = 0.5)
#' meanfieldMutual$sim_slv2()
MeanFieldMutual <- R6Class('MeanFieldMutual',
  public = list(
    n1 = NULL,
    n2 = NULL,
    n = NULL,
    s = NULL,
    kc = NULL,
    c = NULL,
    km = NULL,
    m = NULL,
    h = NULL,
    graphc = NULL,
    graphm = NULL,
    rho = NULL,
    Delta_real = NULL,
    Delta_est = NULL,
    alpha_real = NULL,
    alpha_est = NULL,
    rmin = NULL,
    refSim = NULL,  # referenced simulation object
    simulated_sde = FALSE,  # if the SDE has been simulated
    get_params = function() {
      list(n1 = self$n1, n2 = self$n2, n = self$n, s = self$s, kc = self$kc, c= self$c,
           km = self$km, m = self$m, h = self$h, rho = self$rho, Delta = self$Delta_real,
           alpha = self$alpha_real, rmin = self$rmin, slv2_rmax = self$refSim$r,
           slv2_rmin = self$refSim$rmin, slv2_sigma = self$refSim$sigma,
           slv2_steps = self$refSim$steps, slv2_stepwise = self$refSim$stepwise)
    },
    update = function(s = NULL, c = NULL, km = NULL, m = NULL, h = NULL) {
      if (! is.null(s)) self$s = s
      if (! is.null(c)) self$c = c
      if (! is.null(km)) {
        if (km != self$km)
          self$km = km
          self$graphm <- Graph$new(type = 'bipartite_regular', n1 = self$n1, k = self$km, directed = T, is_adj = T)
      }
      if (! is.null(m)) self$m = m
      if (! is.null(h)) self$h = h
      self$rho <- private$get_rho(self$s, self$c, self$kc, self$m, self$km)
      self$Delta_real <- private$get_Delta_real(self$c, self$kc, self$m, self$km, self$graphm$get_graph())
      self$Delta_est <- private$get_Delta_est(self$c, self$kc, self$m, self$km, self$n)
      self$alpha_real <- self$Delta_real / self$rho
      self$alpha_est <- self$Delta_est / self$rho
      self$rmin <- private$get_rmin(self$s, self$c, self$kc, self$m, self$km, self$h)
    },
    initialize = function(n1 = 20, s = 1, c = 0.01, km = 5, m = 1, h = 0.5) {
      self$n1 <- n1
      self$n2 <- n1
      self$n <- self$n1 + self$n2
      self$s <- s
      self$kc <- self$n1 - 1
      self$c <- c
      self$km <- km
      self$m <- m
      self$h <- h
      self$graphc <- Graph$new(type = 'two_blocks_regular', n1 = self$n1, n2 = self$n2, k = self$kc)
      self$graphm <- Graph$new(type = 'bipartite_regular', n1 = self$n1, k = self$km, directed = T, is_adj = T)
      self$rho <- private$get_rho(self$s, self$c, self$kc, self$m, self$km)
      self$Delta_real <- private$get_Delta_real(self$c, self$kc, self$m, self$km, self$graphm$get_graph())
      self$Delta_est <- private$get_Delta_est(self$c, self$kc, self$m, self$km, self$n)
      self$alpha_real <- self$Delta_real / self$rho
      self$alpha_est <- self$Delta_est / self$rho
      self$rmin <- private$get_rmin(self$s, self$c, self$kc, self$m, self$km, self$h)
    },
    sim_slv2 = function(rmax = 0, steps = 1000, stepwise = 1, sigma = 0.02, extension = 1.1, xinit_sd = 0.1) {
      self$simulated_sde = FALSE # reset the state to not simulated
      # rmin = private$get_rmin(self$s, self$c, self$kc, self$m, self$km, self$h)
      rmin = self$rmin * extension
      terminal_time = steps * stepwise  # time interval
      delta = (rmax - rmin) / terminal_time # increment of parameter at each step
      x <- self$get_xstars(rmax, self$s, self$c, self$kc, self$m, self$km, self$h)
      x1 <- x$X1
      xinit = runif2(self$n, x1, x1 * xinit_sd) 
      
      self$refSim <- SimSLV2$new()
      self$refSim$r = rmax
      self$refSim$rmin = rmin
      self$refSim$sigma = sigma
      self$refSim$extension = extension
      self$refSim$xinit_sd = xinit_sd
      self$refSim$set_drift(n1 = self$n1, n2 = self$n2)
      self$refSim$set_diffusion(n = self$n)
      self$refSim$set_variables(n = self$n)
      self$refSim$set_model()
      self$refSim$set_times(steps = steps, stepwise = stepwise)
      self$refSim$set_init(xinit = xinit)
      self$refSim$set_params(n1 = self$n1, n2 = self$n2, r = rmax, delta = delta, s = self$s, c = self$c, h = self$h, M = self$m * self$graphm$get_graph(), sigma = sigma)
      self$refSim$simulate()
      self$simulated_sde = TRUE  # set the state to have been simulated
#      self$out = self$refSim$get_out()
    },
    # the two possible equilibrium abundances
    get_xstars = function(r, s, c, kc, m, km, h) {
      A = (s + kc * c) * h * km * m
      B = (s + kc * c) - km * m - r * h * km * m
      C = - r
      delta = B^2 - 4 * A * C
      if (delta >= 0) {
        X1 = (-B + sqrt(B^2 - 4 * A * C)) / (2 * A)
        X2 = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
        return(list(X1 = X1, X2 = X2))
      }
      else {
        warning('delta = B^2 - 4AC < 0')
        return(list(X1 = NaN, X2 = NaN))
      }
    },
    # the dot eigenvalue and the real and estimated semicircle eigenvalue of the Shadow Jacobian matrix and the Jacobian matrix
    get_dot_semicircle = function(r, s, c, kc, m, km, h, n, graphm) {
      m_semicircle_real <- private$get_m_semicircle_real(graphm)
      m_semicircle_est <- private$get_m_semicircle_est(n, km)
      m_tilde <- private$get_m_tilde(r, s, c, kc, m, km, h)
      x <- self$get_xstars(r, s, c, kc, m, km, h)$X1 # the stable equilibrium
      # the dot eigenvalue of shadow-Jacobian
      Jshadow_dot <- km * m_tilde - kc * c - s
      # the semicircle eigenvalue of shadow-Jacobian
      Jshadow_semicircle_real <- m_semicircle_real * m_tilde + c - s
      Jshadow_semicircle_est <- m_semicircle_est * m_tilde + c - s
      # the dot eigenvalue of Jacobian
      J_dot <- Jshadow_dot * x
      # the semicircle eigenvalue of Jacobian
      J_semicircle_real <- Jshadow_semicircle_real * x
      J_semicircle_est <- Jshadow_semicircle_est * x
      list(Jshadow_dot = Jshadow_dot, Jshadow_semicircle_real = Jshadow_semicircle_real, Jshadow_semicircle_est = Jshadow_semicircle_est, J_dot = J_dot, J_semicircle_real = J_semicircle_real, J_semicircle_est = J_semicircle_est)  # , Vc = Vc, Vs = Vs, asyn = asyn
}

  ),
  private = list(
    # the ratio between mutualistic strength and competitive strength
    get_rho = function(s, c, kc, m, km) {
      return(km * m / (s + kc * c))
    },
    # the real semicircle eigenvalue of the mutualistic adjacency matrix
    get_m_semicircle_real = function(graphm) {
      semicircle_real <- sort(eigen(graphm)$values, decreasing = TRUE)[2]
    },
    # the estimated semicircle eigenvalue of the mutualistic adjacency matrix
    get_m_semicircle_est = function(n, km) {
      semicircle_est <- estimate_semicircle_bipartite_regular(n, km)$est
    },
    # the real ratio between mutual spectral gap and competitive spectral gap
    get_Delta_real = function(c, kc, m, km, graphm) {
      semicircle_real <- private$get_m_semicircle_real(graphm)
      Delta_real <- ((km - semicircle_real) * m) / ((kc + 1) * c)
      Delta_real
    },
    # the estimated ratio between mutual spectral gap and competitive spectral gap
    get_Delta_est = function(c, kc, m, km, n) {
      semicircle_est <- private$get_m_semicircle_est(n, km)
      Delta_est <- ((km - semicircle_est) * m) / ((kc + 1) * c)
      Delta_est
    },
    # the dot eigenvalue
    get_dot = function(r, rho, h) {
      focus <- sqrt((rho + r * h * rho + 1)^2 - 4 * rho)  # B^2-4AC
      numerator <- (rho + r * h * rho - 1 + focus) * focus
      denominator <- h * rho * (rho + r * h * rho + 1 + focus)
      - numerator / denominator
    },
    # the minimal value of [r] where critical transitions happen
    get_rmin = function(s, c, kc, m, km, h) {
      return(-(sqrt(km * m) - sqrt(s + kc * c))^2 / (h * km * m))
    },
    # the effective mutualistic strength caused by h > 0
    get_m_tilde = function(r, s, c, kc, m, km, h) {
      x <- self$get_xstars(r, s, c, kc, m, km, h)$X1
      m_tilde <- m  / (1 + h * km * m * x)^2
      return(m_tilde)
    }
  ),
  cloneable = FALSE
  )