library(R6)

#' @title The root class of Model for Mutualistic systems.
#' @description This class has restrictions. 
#' First, species have to be splitted to two groups, mutualistic interactions happen between groups while competitive interactions happen within groups, thus, the graph of mutualistic interactions is bipartite. 
#' Second, method \code{sim_slv2} using class \code{SimSLV2} to simulate SLV2 for mutualistic system, has limitation because of the limitation of class \code{SimSLV2}.
#' @field n1 species number of group 1
#' @field n2 species number of group 2
MutualModel <- R6Class('MutualModel',
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
    s.sd = 0,
    c.sd = 0,
    m.sd = 0,
    h.sd = 0,
    delta = 0,
    graphc = NULL,
    graphm = NULL,
    simSLV2 = NULL,  # referenced SLV2 simulation object
    simLV2 = NULL,   # referenced LV2 simulation object
    simLV2Press = NULL,   # referenced LV2 simulation object
    state_slv2 = FALSE,  # if the SLV2 has been simulated
    state_lv2 = FALSE,   # if the LV2 has been simulated
    state_lv2_press = FALSE,   # if the LV2 Press has been simulated
    initialize = function(n1, n2, s, kc, c, km, m, h, s.sd, c.sd, m.sd, h.sd, delta, graphc, graphm) {
      #cat('Initialize MutualModel object.\n')
      self$n1 <- n1
      self$n2 <- n2
      self$n <- self$n1 + self$n2
      self$s <- s
      self$kc <- kc # self$n1 - 1
      self$c <- c
      self$km <- km
      self$m <- m
      self$h <- h
      self$s.sd <- s.sd
      self$c.sd <- c.sd
      self$m.sd <- m.sd
      self$h.sd <- h.sd
      self$graphc <- graphc
      self$graphm <- graphm
    },
    sim_slv2 = function(rmax, delta, steps, stepwise, sigma, xinit) {
      cat('Simulate SLV2.\n')
      self$state_slv2 = FALSE # reset the state to not simulated
      self$simSLV2 <- SimSLV2$new()
      self$simSLV2$sim(
        n1 = self$n1, n2 = self$n2, r = rmax, delta = delta, s = self$s,
        c = self$c, h = self$h, M = self$m * self$graphm$get_graph(), 
        sigma = sigma, steps = steps, stepwise = stepwise, xinit = xinit)
      self$state_slv2 = TRUE  # set the state to have been simulated
      # self$out = self$simSLV2$get_out()
    },
    sim_lv2 = function(r, r.sd =0, steps = 100, stepwise = 1, xinit, method = c('lsoda', 'lsode'),
                       jactype = c('fullint', 'fullusr'),
                       atol = 1e-8, rtol = 1e-8,
                       extinct_threshold = 1e-6) {
      #cat('Simulate LV2.\n')
      self$state_lv2 = FALSE 
      coeffs <- list(n1 = self$n1, n2 = self$n2, r.row.mu = r, r.row.sd = r.sd, r.col.mu = r, r.col.sd = r.sd, s.mu = self$s, s.sd = self$s.sd, c.mu = self$c, c.sd = self$c.sd, m.mu = self$m, m.sd = self$m.sd, h.mu = self$h, h.sd = self$h.sd, delta = self$delta)
      self$simLV2 <- SimLV2Mutual$new(method, jactype, atol, rtol, extinct_threshold)
      self$simLV2$sim(steps, stepwise, xinit, coeffs, graphc = self$graphc$get_graph(), graphm = self$graphm$get_graph())
      self$state_lv2 = TRUE 
    },
    sim_lv2_press = function(
      # parameters for pressure
      perturb_type = c('growthrate_all'), perturb_num, r.delta.mu, r.delta.sd = 0, xstars.immigration = 0, xstars.sd = 0, is.out = FALSE, r, r.sd =0, steps = 100, stepwise = 1, xinit, method = c('lsoda', 'lsode'), jactype = c('fullint', 'fullusr'), atol = 1e-8, rtol = 1e-8, extinct_threshold = 1e-8) {
      cat('Simulate LV2 Press.\n')
      self$state_lv2_press = FALSE 
      coeffs <- list(n1 = self$n1, n2 = self$n2, r.row.mu = r, r.row.sd = r.sd, r.col.mu = r, r.col.sd = r.sd, s.mu = self$s, s.sd = self$s.sd, c.mu = self$c, c.sd = self$c.sd, m.mu = self$m, m.sd = self$m.sd, h.mu = self$h, h.sd = self$h.sd, delta = self$delta)
      self$simLV2Press <- SimLV2Mutual$new(method, jactype, atol, rtol, extinct_threshold)
      self$simLV2Press$sim_press(perturb_type,perturb_num, r.delta.mu, r.delta.sd, xstars.immigration, xstars.sd, is.out, steps, stepwise, xinit, coeffs, graphc = self$graphc$get_graph(), graphm = self$graphm$get_graph())
      self$state_lv2_press = TRUE 
    }
))