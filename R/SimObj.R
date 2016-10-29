library(R6)

#' @title The root class(interface) for all simulation objects
#' @description A classical simulation includes several steps: 
#' 1. set model, such as Ordinary Differential Equaitons, Stochastic Differential Equation, or some kind of agent-based model.
#' 2. set timesteps, [steps] with [stepwise] increment.
#' 3. set initial values, for all the state variables.
#' 4. set parameters, which adhere to the model.
#' 5. start simulate, which will produce the output data.
#' 6. get output of the simulation, in order for further data analysis.
SimObj <- R6Class('SimObj', 
  public = list(
    set_model = function() {
      cat('Set model for SimObj object.\n')
    },
    set_times = function() {
      cat('Set timesteps for SimObj object.\n')
    },
    set_init = function(xinit) {
      cat('Set initial values for SimObj object.\n')
      self$xinit = xinit
    },
    set_params = function() {
      cat('Set parameters for SimObj object.\n')
    },
    sim = function() {
      cat('Simulate for SimObj object, produce the output.\n')
    },
    get_out = function() {
      cat('Get the output of SimObj object.\n')
    },
    initialize = function() {
      cat('Initialize SimObj object.\n')
    },
#  ),
#  private = list(
    # reference to other object that implement the real simulation
    refObj = NULL,
    # current state of the simulation object: is_setted_(model, times, init, params, out)
    state = NULL,
    steps = NULL,
    stepwise = NULL,
    xinit = NULL,
    params = NULL
  ))