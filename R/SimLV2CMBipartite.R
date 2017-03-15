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
#' @param coeffs a list of coefficients : 
#' \describe{
#' \item{n1, n2}{node number of two groups}
#' \item{r.row.mu, r.row.sd}{mean and sd of intrinsic growth rates of group 1}
#' \item{r.col.mu, r.col.sd}{mean and sd of intrinsic growth rates of group 2}
#' \item{s.mu, s.sd}{mean and sd of self-regulation strength}
#' \item{c.mu, c.sd}{mean and sd of strength of competitive interactions}
#' \item{m.mu, m.sd}{mean and sd of strength of mutualistic interactions}
#' \item{h.mu, h.sd}{mean and sd of handling time}
#' \item{delta}{trade-off between strength and number of mutualistic interactions}
#' }
#' list(n1 = n1, n2 = n2, r.row.mu = r, r.row.sd = r.sd, r.col.mu = r, r.col.sd = r.sd, s.mu = s, s.sd = s.sd, c.mu = c, c.sd = c.sd, m.mu = m, m.sd = m.sd, h.mu = h, h.sd = h.sd, delta = delta)
#' @param graphc, the adjacency matrix of competitive interactions
#' @param graphm, the adjacency matrix of mutualistic interactions
#' @return a list of parameters for model \code{\link{model_lv2_cm}}
params_lv2_cm <- function(coeffs, graphc, graphm) {
  with(coeffs, {
    n = n1 + n2
    C = graphc # the competition part
    C = runif2(n * n, c.mu, c.sd) * C
    diag(C) = runif2(n, s.mu, s.sd)
    
    M = graphm  # the mutualistic part
    edges = sum(M > 0)  # the number of all mutualistic interactions(edges)
    degrees = rowSums(M) # the degrees of all species
    M[M > 0] = runif2(edges, m.mu, m.sd) # strength of interspecies cooperation
    total_strength_old <- sum(M)
    M = M / degrees^delta  # trade-off of mutualistic strength and number
    # in order to keep the total strength constant
    # we need multiple according to delta
    total_strength_new <- sum(M)
    M = M * (total_strength_old / total_strength_new)
    
    h = runif2(n, h.mu, h.sd)
    r = c(runif2(n1, r.row.mu, r.row.sd), runif2(n2, r.col.mu, r.col.sd))
    params <- list(r = r, C = C, M = M, h = h)  # the [parms] of ode model
    params
  })
}

#' @title Simulation for LV2 model
#' @examples 
#' SimLV2CMBipartite <- SimLV2CMBipartite$new()
#' SimLV2CMBipartite$sim(steps, stepwise, xinit, coeffs, graphc, graphm)
SimLV2CMBipartite <- R6Class('SimLV2CMBipartite',
  inherit = SimODE,
  public = list(
    coeffs = NULL,
    graphc = NULL,
    graphm = NULL,
    set_model = function() {
      #super$set_model()
      self$model <- model_lv2_cm
    },
    set_params = function(coeffs, graphc, graphm) {
      #super$set_params()
      self$coeffs = coeffs
      self$graphc = graphc
      self$graphm = graphm
      self$params = params_lv2_cm(coeffs, graphc, graphm)
    },
    simulate = function(steps, stepwise, xinit, coeffs, graphc, graphm) {
      self$set_model()
      self$set_times(steps = steps, stepwise = stepwise)
      self$set_init(xinit = xinit)
      self$set_params(coeffs = coeffs, graphc = graphc, graphm = graphm)
      super$sim()
    }
  ))

