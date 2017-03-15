library(R6)

#' @title The root class(interface) for all simulation objects
#' @description A classical simulation includes several steps: 
#' \describe{
#' \item{1.}{set Model, such as Ordinary Differential Equaitons, Stochastic Differential Equation, or some kind of agent-based model.}
#' \item{2.}{set Timesteps, [steps] with [stepwise] increment.}
#' \item{3.}{set Initial values, for all the state variables.}
#' \item{4.}{set Parameters, which adhere to the Model.}
#' \item{5.}{start Simulation, which will produce the output data.}
#' \item{6.}{get Output of the Simulation, in order for further data analysis.}
#' }
#' @field model, times, xinit, params, out
#' @field status, current status of the simulation object, status can only change from 0 to 5:
#' \describe{
#' \item{0}{the simulation object is initialized.}
#' \item{1}{the Model is setted.}
#' \item{2}{the Times is setted.}
#' \item{3}{the Initial values of all state variables are setted.}
#' \item{4}{the Parameters of Model is setted.}
#' \item{5}{the simulation is correctly completed, and the Output is gotten.}
#' }
SimObj <- R6Class('SimObj', 
  public = list(
    initialize = function() {
      #self$status = 0
      cat('Initialize SimObj object...\n')
    },
    set_model = function() {
      # current status of the simulation object must be Initialized(0)
      #stopifnot(self$status == 0) 
      #self$status = 1
      cat('Set model for SimObj object...\n')
    },
    set_times = function() {
      # current status of the simulation object must be Model-setted(1)
      #stopifnot(self$status == 1) 
      #self$status = 2
      cat('Set timesteps for SimObj object...\n')
    },
    set_init = function(xinit) {
      # current status of the simulation object must be Times-setted(2)
      #stopifnot(self$status == 2) 
      #self$status = 3
      #cat('Set initial values for SimObj object...\n')
      self$xinit = xinit
    },
    set_params = function() {
      # current status of the simulation object must be InitialValues-setted(3)
      #stopifnot(self$status == 3) 
      #self$status = 4
      cat('Set parameters for SimObj object...\n')
    },
    sim = function() {
      # current status of the simulation object must be Params-setted(4)
      #stopifnot(self$status == 4) 
      #self$status = 5
      cat('Simulate for SimObj object, produce the output...\n')
    },
    get_out = function() {
      cat('Get the output of SimObj object.\n')
    },
    # reference to other object that implement the real simulation
    refObj = NULL,
    
    model = NULL,
    times = NULL,
    steps = NULL,
    stepwise = NULL,
    xinit = NULL,
    params = NULL,
    out = NULL,
    
    status = NULL
  ))