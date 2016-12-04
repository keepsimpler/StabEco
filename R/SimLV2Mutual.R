library(R6)
library(deSolve)

##########################Competition-Mutualism##############################
#' @title Lotka-Volterra (LV) Equations of Holling type II for a community mixed by Competition and Mutualism interactions
#' @param time, time steps of simulation
#' @param init, the initial state of the LV system, a vector
#' @param params, parameters passed to LV model, a list of:
#' \describe{
#'   \item{r}{a vector of the intrinsic growth rates of species}
#'   \item{C}{a matrix of intra-species and inter-species competitions}
#'   \item{M}{a matrix of mutualism interactions among species}
#'   \item{h}{the saturate coefficient, handling time of species feed}
#' }
#' @return the derivation
#' @import deSolve
model_lv2_cm <- function(time, init, params, ...) {
  N = init  # initial state
  with(params, {
    dN <- N * ( r - C %*% N + (M %*% N) / (1 + h * M %*% N) )
    list(c(dN))
  })
}

#' @title Jacobian function for model \code{model_lv2_cm}
jacfunc_lv2_cm <- function(time, init, params) {
  N = init
  with(params, {
    s <- diag(C)
    diag(C) <- 0
    fixi <- diag(c(r - 2 * s * N - C %*% N + (M %*% N) / (1 + h * M %*% N)))
    fixj_m <- diag( N / c( (1 + h * M %*% N)^2 ) ) %*% M
    fixj_c <- - diag(N) %*% C
    J <- fixi + fixj_m + fixj_c
    return(J)
  })
}

#' @title parameters for model \code{\link{model_lv2_cm}}
#' @description assign parameters for model \code{\link{model_lv2_cm}} according to a structural network and a couple of coefficients
#' @param coeffs a list of coefficients : list(n1 = n1, n2 = n2, r.row.mu = r, r.row.sd = r.sd, r.col.mu = r, r.col.sd = r.sd, s.mu = s, s.sd = s.sd, c.mu = c, c.sd = c.sd, m.mu = m, m.sd = m.sd, h.mu = h, h.sd = h.sd, delta = delta)
#' @param graphc, the adjacency matrix of competitive interactions
#' @param graphm, the adjacency matrix of mutualistic interactions
#' @return a list of parameters for model \code{\link{model_lv2_cm}}
params_lv2_cm <- function(coeffs, graphc, graphm) {
  with(coeffs, {
    n = n1 + n2
    C = graphc # the competition part
    C = runif2(n * n, c.mu, c.sd) * C
    diag(C) = runif2(n, s.mu, s.sd)
    
    M = graphm
    edges = sum(M > 0)  # the number of all mutualistic interactions(edges)
    degrees = rowSums(M) # the degrees of all species
    M[M > 0] = runif2(edges, m.mu, m.sd) # values of interspecies cooperation
    total_strength_old <- sum(M)
    M = M / degrees^delta  # trade-off of mutualistic strength and number
    # in order to keep the total strength constant
    # we need multiple according to delta
    total_strength_new <- sum(M)
    M = M * (total_strength_old / total_strength_new)
    
    h = runif2(n, h.mu, h.sd)
    r = c(runif2(n1, r.row.mu, r.row.sd), runif2(n2, r.col.mu, r.col.sd))
    self$params <- list(r = r, C = C, M = M, h = h)  # the [parms] of ode model
  })
}

#' @title Simulation for LV2 model
#' @examples 
#' simLV2Mutual <- SimLV2Mutual$new()
#' simLV2Mutual$sim(steps, stepwise, xinit, coeffs, graphc, graphm)
SimLV2Mutual <- R6Class('SimLV2Mutual',
  inherit = SimODE,
  public = list(
    perturb_type = NULL,
    perturb_func = NULL,
    perturb_num = NULL,
    r.delta.mu = NULL,
    r.delta.sd = 0,
    xstars.immigration = 0,
    xstars.sd = 0,
    is.out = FALSE,
    rmax = NULL,
    out_press = NULL, # output of pressure simulation
    
    initialize = function(method = c('lsoda', 'lsode'),
                          jactype = c('fullusr', 'fullint'),
                          atol = 1e-8, rtol = 1e-8,
                          extinct_threshold = 1e-6) {
      self$model <- model_lv2_cm
      self$jacfunc <- jacfunc_lv2_cm
      self$method <- match.arg(method)
      self$jactype <- match.arg(jactype)
      self$atol <- atol
      self$rtol <- rtol
      self$extinct_threshold <- extinct_threshold
    },
    set_params = params_lv2_cm,
    sim = function(steps, stepwise, xinit, coeffs, graphc, graphm) {
      self$set_times(steps = steps, stepwise = stepwise)
      self$set_init(xinit = xinit)
      self$set_params(coeffs = coeffs, graphc = graphc, graphm = graphm)
      super$sim()
    },
    # simulate pressure effect on [r]
    sim_press = function(
      # parameters for pressure
      perturb_type = c('growthrate_all'), perturb_num, r.delta.mu, r.delta.sd = 0, xstars.immigration = 0, xstars.sd = 0, is.out = FALSE,
      # parameters for autonormous
      steps, stepwise, xinit, coeffs, graphc, graphm) {
      self$perturb_type = match.arg(perturb_type)
      self$perturb_num = perturb_num
      self$perturb_func = perturb
      self$r.delta.mu = r.delta.mu
      self$r.delta.sd = r.delta.sd
      self$xstars.immigration = xstars.immigration
      self$xstars.sd = xstars.sd
      self$is.out = is.out
      self$rmax = coeffs$r.row.mu # maximal r 
      perturb_coeffs = list(r.delta.mu = r.delta.mu, r.delta.sd = r.delta.sd, xstars.immigration = xstars.immigration, xstars.sd = xstars.sd)
      
      self$set_times(steps = steps, stepwise = stepwise)
      self$set_init(xinit = xinit)
      self$set_params(coeffs = coeffs, graphc = graphc, graphm = graphm)
      
      out_press = list()
      for (i in 1:perturb_num) {
        super$sim()
        species.survived = which(self$xstars > 0)  # survived species
        flag = 0
        # if all species are extinct, will end the simulation
        if (length(species.survived) == 0) flag = 1
        J = jacfunc_lv2_cm(time = self$times, init = self$xstars, params = self$params)
        if (is.out) {
          ret = list(out = self$out, xstars = self$xstars, J = J, params = self$params, species.survived = species.survived, flag = flag)
        }
        else {
          ret = list(xstars = self$xstars, J = J, params = self$params, species.survived = species.survived, flag = flag)
        }
        out_press[[length(out_press) + 1]] = ret
        # if all species are extinct, end the simulation
        if (flag == 1)  # || flag == 2
          break;
        
        # perturbation that returns new parameters and initial values
        perturb.ret = perturb(self$params, self$xstars, perturb_coeffs, perturb_type)
        self$params = perturb.ret$params
        self$xinit = perturb.ret$xstars
      }
      self$out_press = out_press
      
    }
  ))

