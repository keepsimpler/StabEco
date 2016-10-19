library(R6)
library(yuima) # Simulation of Stochastic Differential Equaitons

#' @title Simulation for (Multivariate) Stochastic Process implemented by \code{Yuima} package inherited from Interface \code{SimObj}
#' @field of a stochastic process:
#' \describe{
#' \item{drift}{the drift vector}
#' \item{diffusion}{the diffusion matrix}
#' \item{variables}{name of state variables}
#' }
SimYuima <- R6Class('SimYuima',
  inherit = SimObj,
  public = list(
    drift = NULL,
    diffusion = NULL,
    variables = NULL,
    initialize = function() {
      cat('Initialize the SimYuima object.')
      self$refObj = setYuima()
      self$state = 0
    },
    set_drift = function() {
      cat('Set drift vector for SimYuima object')
    },
    set_diffusion = function() {
      cat('Set diffusion matrix for SimYuima object')
    },
    set_variables = function() {
      cat('Set names of variables for SimYuima object')
    },
    set_model = function() {
      cat('Set model for SimYuima object.')
      model = setModel(drift = self$drift, diffusion = self$diffusion, solve.variable = self$variables, state.variable = self$variables)
      self$refObj@model = model
      self$state = 1
    },
    set_times = function(steps, stepwise) {
      cat('Set timesteps for SimYuima object.')
      sampling = setSampling(Terminal = steps * stepwise, n = steps)
      self$refObj@sampling = sampling
      self$steps = steps
      self$stepwise = stepwise
      self$state = 2
    },
    simulate = function() {
      cat('simulate SimYuima object, produce the output.')
      self$refObj = simulate(self$refObj, true.parameter = self$params, xinit = self$xinit)
    },
    get_out = function() {
      cat('Get the output of SimYuima object')
      self$refObj@data@zoo.data
    }
  ))