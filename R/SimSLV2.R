library(R6)
library(yuima) # Simulation of Stochastic Differential Equaitons

#' @title Simulation for Stochastic LV2 model
#' @description inherit from class \code{SimYuima}
#' dX = X (r - C X + M X / (1 + h M X)) dt + Sigma dW
#' X (r - C X + M X / (1 + h M X)) is the drift vector,
#' which include:
#'  r - vector of intrinsic growth rates
#'  C - competitive interactions matrix, all self-regulation strength is [s], all inter-species competitive strength is [c]
#'  M - mutualistic interactions matrix
#'  h - handling time
#' Note: we only support mean field approximation
#' Note: we only support the bipartite network of mutualistic interactions
#' Sigma is the diffusion matrix, we assume it's a diagonal matrix
#' W is a vector of Wiener process
#' @examples 
#'       simSLV2 <- SimSLV2$new()
#'       simSLV2$set_drift(n1 = 1, n2 = 2)
#'       simSLV2$set_diffusion(n = 1+2)
#'       simSLV2$set_variables(n = 1+2)
#'       simSLV2$set_model()
#'       simSLV2$set_times(steps = 1000, stepwise = 1)
#'       simSLV2$set_init(xinit = c(1,1,1))
#'       simSLV2$set_params(n1 = 1, n2 = 2, r = 1, delta = 0.001, s = 1, c = 0.01, h = 0.5, M = 1 * matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0), ncol = 3), sigma = 0.01)
#'       simSLV2$simulate()
#'       out = simSLV2$get_out()

SimSLV2 <- R6Class('SimSlv2',
  inherit = SimYuima,
  public = list(
    set_drift = function(n1, n2) {
      cat('Set drift vector for SimSLV2 object')
      rows <- 1:(n1 + n2)
      drift <- sapply(rows, function(row) {
        if (row <= n1) {
          # competitive block (1:n1, 1:n1)
          compet <- 1:n1
          # mutualistic block ( 1:n1, (n1+1):(n1+n2) )
          mutual <- (n1+1):(n1+n2) 
        } else {
          # competitive block ( (n1+1):(n1+n2), (n1+1):(n1+n2) )
          compet <- (n1+1):(n1+n2) 
          # mutualistic block ( (n1+1):(n1+n2), 1:n1 )
          mutual <- 1:n1 
        }
        if (length(compet[- which(compet == row)]) > 0) {
          C <- paste('c*(', paste('x', compet[-which(compet == row)],
                                  sep = '', collapse = '+'), ')', sep = '')
        } else {
          C <- ""
        }
        M <- paste('m', row, mutual, '*x', mutual, sep = '', collapse = '+')
        if (C == "") {
          drift <- paste('x', row, '*(r - delta * t - s * x', row, ' + (', M,
                         ')/(1+h*(', M, ')', '))', sep = '')
        } else {
          drift <- paste('x', row, '*(r - delta * t - s * x', row, ' - ', C, 
                         ' + (', M, ')/(1+h*(', M, ')', '))', sep = '')
        }
      })
      self$drift = drift
    },
    set_diffusion = function(n) {
      cat('Set diffusion matrix for SimSLV2 object')
      diffusion <- apply(diag(1:n), c(1, 2), function(ij)
        if (ij != 0)
          'sigma'
        else '0'
      )
      self$diffusion = diffusion
    },
    set_variables = function(n) {
      cat('Set names of variables for SimSLV2 object')
      variables <- sapply(1:n, function(i) paste('x', i, sep = ''))
      self$variables = variables
    },
    set_params = function(n1, n2, r, delta, s, c, h, M, sigma) {
      self$r = r
      self$delta = delta
      self$sigma = sigma
      cat('Set parameters for SimSLV2 object.')
      params <- list(r = r, delta = delta, s = s, c = c, h = h)
      params_M <- as.list(M[1:n1, (n1+1):(n1+n2)]) # by column
      Mij <- outer(1:n1, (n1+1):(n1+n2), FUN = paste, sep='')
      names(params_M) <- sapply(Mij, function(ij) paste('m', ij, sep = '')) # by column
      params <- c(params, params_M)
      params_M <- as.list(M[(n1+1):(n1+n2), 1:n1]) # by column
      Mij <- outer((n1+1):(n1+n2), 1:n1, FUN = paste, sep='')
      names(params_M) <- sapply(Mij, function(ij) paste('m', ij, sep = '')) # by column
      params <- c(params, params_M)
      params <- c(params, list(sigma = sigma))
      self$params = params
    },
    sim = function(n1, n2, r, delta, s, c, h, M, sigma, steps, stepwise, xinit) {
      # delta == (r - rmin) / (steps * stepwise)
      self$rmin = r - delta * steps * stepwise
      self$set_drift(n1, n2)
      self$set_diffusion(n1 + n2)
      self$set_variables(n1 + n2)
      self$set_model()
      self$set_times(steps, stepwise)
      self$set_init(xinit)
      self$set_params(n1, n2, r, delta, s, c, h, M, sigma)
      super$sim()
    },
    r = 0,  # the start value of [r]
    delta = 0, # decrease of [r] at each step, == (r - rmin) / (steps * stepwise)
    sigma = 0.02,
    rmin = NULL # the end value of [r]
  ))