library(R6)

#' @title A Mean Field Model for Mutualistic systems 
#' @description inherit from class \code{MutualModel}, add two restrictions:
#' First, the number of species of two groups are same, n1 == n2
#' Second, the competition within groups is full, kc == n1-1
#' Third, all the variances are zero, s.sd = c.sd = m.sd = h.sd = 0
#' @examples 
#' meanfieldMutual <- MeanFieldMutual$new(n1 = 20, s = 1, c = 0.01, km = 5, m = 1, h = 0.5)
#' meanfieldMutual$sim_slv2()
MeanFieldMutual <- R6Class('MeanFieldMutual',
  inherit = MutualModel,
  public = list(
    rho = NULL,
    Delta_real = NULL,
    Delta_est = NULL,
    alpha_real = NULL,
    alpha_est = NULL,
    rmin = NULL,
    extension = 1.1, # extension of simulation, in order to ensure the happen of bifurcations
    xinit_sd = 0.1, # the sd of initial values
    get_params = function() {
      list(n1 = self$n1, n2 = self$n2, n = self$n, s = self$s, kc = self$kc, c= self$c,
           km = self$km, m = self$m, h = self$h, rho = self$rho, Delta = self$Delta_real,
           alpha = self$alpha_real, rmin = self$rmin, slv2_rmax = self$simSLV2$r,
           slv2_rmin = self$simSLV2$rmin, slv2_sigma = self$simSLV2$sigma,
           slv2_steps = self$simSLV2$steps, slv2_stepwise = self$simSLV2$stepwise)
    },
    update = function(s = NULL, c = NULL, km = NULL, m = NULL, h = NULL) {
      if (! is.null(s)) self$s = s
      if (! is.null(c)) self$c = c
      if (! is.null(km)) {
        if (km != self$km)
          self$km = km
          self$graphm <- BiGraph$new(type = 'bipartite_regular', n1 = self$n1, k = self$km, directed = T, is_adj = T)
      }
      if (! is.null(m)) self$m = m
      if (! is.null(h)) self$h = h
      private$update_inner()
    },
    initialize = function(n1 = 20, s = 1, c = 0.01, km = 5, m = 1, h = 0.5, delta = 0, xinit_sd = 0.1, graphc = NULL, graphm = NULL) {
      if (is.null(graphc)) {
        graphc <- TwoblocksGraph$new(type = 'two_blocks_regular', n1 = n1, n2 = n1, k = n1 - 1)
      }
      else {
        if (class(graphc)[1] == 'TwoblocksGraph') {
          stopifnot(graphc$n1 == n1 && graphc$n2 == n1 && graphc$k == n1 - 1)
        }
        else {
          stop('graphc shoud be TwoblocksGraph class')
          #graphc <- TwoblocksGraph$new(n1 = n1, k = km)$set_graph_inc(graph_inc = graphm)
        }
      }
      if (is.null(graphm)) {
        graphm <- BiGraph$new(type = 'bipartite_regular', n1 = n1, k = km, directed = T, is_adj = T)
      }
      else {
        if (class(graphm)[1] == 'BiGraph') {
          stopifnot(graphm$n1 == n1 && graphm$n2 == n1 && graphm$k == km)
        }
        else {
          tmp <- graphm
          graphm <- BiGraph$new(n1 = n1, k = km)
          graphm$set_graph_inc(graph_inc = tmp)
        }
      }
      self$xinit_sd = xinit_sd
      super$initialize(n1, n2 = n1, s, kc = n1 - 1, c, km, m, h, s.sd = 0, c.sd = 0, m.sd = 0, h.sd = 0, delta, graphc, graphm)
      private$update_inner()
    },
    sim_slv2 = function(rmax = 0, steps = 1000, stepwise = 1, sigma = 0.02, extension = 1.1, xinit_sd = 0.1) {
      # self$state_slv2 = FALSE # reset the state to not simulated
      rmin = self$rmin * extension
      terminal_time = steps * stepwise  # time interval
      delta = (rmax - rmin) / terminal_time # increment of parameter at each step
      x <- get_xstars(rmax, self$s, self$c, self$kc, self$m, self$km, self$h)
      x1 <- x$X1
      xinit = runif2(self$n, x1, x1 * xinit_sd) 
      
      super$sim_slv2(rmax, delta, steps, stepwise, sigma, xinit)
      
      self$extension = extension
      self$xinit_sd = xinit_sd
      #self$state_slv2 = TRUE  # set the state to have been simulated
      # self$out = self$simSLV2$get_out()
    },
    sim_lv2 = function(r, steps = 100, stepwise = 1, method = c('lsoda', 'lsode'),
                       jactype = c('fullint', 'fullusr'),
                       atol = 1e-8, rtol = 1e-8,
                       extinct_threshold = 1e-8) {
      x <- get_xstars(r, self$s, self$c, self$kc, self$m, self$km, self$h)
      x1 <- x$X1
      if (is.nan(x1) || is.na(x1) ||  # r < rmin
          x1 <= 0) # rho < 1
        xinit = runif2(self$n, 1, 1 * self$xinit_sd)
      else
        xinit = runif2(self$n, x1, x1 * self$xinit_sd) 
      super$sim_lv2(r, r.sd = 0, steps, stepwise, xinit, method, jactype, atol, rtol, extinct_threshold)
    },
    sim_lv2_press = function(
      # parameters for pressure
      perturb_type = c('growthrate_all'), perturb_num, r.delta.mu, r.delta.sd = 0, xstars.immigration = 0, xstars.sd = 0, is.out = FALSE,
      r, steps = 100, stepwise = 1, method = c('lsoda', 'lsode'),
                       jactype = c('fullint', 'fullusr'),
                       atol = 1e-8, rtol = 1e-8,
                       extinct_threshold = 1e-8) {
      x <- get_xstars(r, self$s, self$c, self$kc, self$m, self$km, self$h)
      x1 <- x$X1
      xinit = runif2(self$n, x1, x1 * self$xinit_sd) 
      super$sim_lv2_press(perturb_type, perturb_num, r.delta.mu, r.delta.sd, xstars.immigration, xstars.sd, is.out, r, r.sd = 0, steps, stepwise, xinit, method, jactype, atol, rtol, extinct_threshold)
    }
  ),
  private = list(
    # update inner parameters according to external parameters
    update_inner = function() {
      self$rho <- get_rho(self$s, self$c, self$kc, self$m, self$km)
      self$Delta_real <- get_Delta_real(self$c, self$kc, self$m, self$km, self$graphm$get_graph())
      self$Delta_est <- get_Delta_est(self$c, self$kc, self$m, self$km, self$n)
      self$alpha_real <- self$Delta_real / self$rho
      self$alpha_est <- self$Delta_est / self$rho
      self$rmin <- get_rmin(self$s, self$c, self$kc, self$m, self$km, self$h)
    }
  ),
  cloneable = FALSE
  )