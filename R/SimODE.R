library(R6)
library(deSolve)

SimODE <- R6Class('SimODE',
  inherit = SimObj,
  public = list(
    model = NULL, # function of equations of ODE
    times = NULL,
    method = NULL,  # solve method
    jactype = NULL, # Jacobian type
    jacfunc = NULL, # function of compute Jacobian
    atol = 1e-8, # absolute error tolerance
    rtol = 1e-8, # relative error tolerance
    extinct_threshold = 1e-10, # species extinction threshold
    out = NULL,  # output of ODE simulation
    xstars = NULL, # equilibrium values of state variables
    set_times = function(steps, stepwise) {
      cat('Set timesteps for SimODE object.\n')
      self$steps = steps
      self$stepwise = stepwise
      self$times = seq(from = 0, by = stepwise, length.out = steps + 1)
      self$state = 2
    },
    sim = function() {
      #cat('simulate SimODE object, produce the output.\n')
      ode.out = deSolve::ode(self$xinit, self$times, func = self$model, self$params, method = self$method, jactype = self$jactype, jacfunc = self$jacfunc, atol = self$atol, rtol = self$rtol) #  method = "ode45", atol = 1e-13, rtol = 1e-13   # , rootfun = rootfun, method = "lsodar"
      xstars = as.numeric(ode.out[nrow(ode.out), 2:ncol(ode.out)]) # species biomass at equilibrium
      xstars[xstars < self$extinct_threshold] = 0  # species with biomass less than extinct threshold is considered to be extinct
      # if any species' abundance is NaN, that means the ODE dynamic is unstable,
      # the simulation will be ended
      if (any(is.nan(self$xstars))) {
        stop('The ODE dynamic is unstable.')
      }
      self$out = ode.out
      self$xstars = xstars
    }
    
  ))