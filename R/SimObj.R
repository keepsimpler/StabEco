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
      cat('Set model for the simulation.')
    },
    set_times = function() {
      cat('Set timesteps for the simulation.')
    },
    set_init = function() {
      cat('Set initial values for the simulation.')
    },
    set_params = function() {
      cat('Set parameters for the simulation.')
    },
    simulate = function() {
      cat('Actually simulate, produce the output.')
    },
    get_out = function() {
      cat('Get the output of simulation.')
    },
    initialize = function() {
      cat('Initialize the SimObj object.')
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