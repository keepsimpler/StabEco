library(R6)

#' @title Perturbations effect on systems by its parameters [params] and state variables [xstars]
#' @param params parameters
#' @param xstars state variables(species abundances)
#' @param type perturbation types:
#' \describe{
#'   \item{growthrate_all}{increasing/decreasing the intrinsic growth rates of ALL species}
#'   \item{growthrate_part}{increasing/decreasing the intrinsic growth rates of a part of species}
#'   \item{mutualistic_strength}{increasing/decreasing strengths of mutualistic interactions}
#'   \item{competitive_strength}{increasing/decreasing strengths of competitive interactions}
#'   \item{primary_extinct}{remove one species}
#' }
#' @param coeffs list of coefficients of perturbation:
#' \describe{
#'   \item{r.delta.mu, r.delta.sd}{stepwise of inc/dec of intrinsic growth rate}
#'   \item{xstars.sd}{randomly perturb species abundances at the end of each step of ode simulation to simulate(explore) the Initial Value Problem}
#'   \item{xstars.immigration}{simulating a immigration factor to (re)establish extinct species}
#' }
perturb <- function(params, xstars, coeffs, type = c('growthrate_all')) {
  type <- match.arg(type)
  with(coeffs, {
    if (type == 'growthrate_all') {
      params$r = params$r - runif2(length(params$r), r.delta.mu, r.delta.sd)
    }
    if (xstars.immigration > 0)
      xstars[xstars == 0] <- xstars[xstars == 0] +
        runif2(length(xstars[xstars == 0]), xstars.immigration, 0.5 * xstars.immigration)
    # random perturb xstars
    if (xstars.sd > 0) {
      #set.seed(123) # fix the random number to check its influence
      xstars <- xstars * runif2(length(xstars), 1, xstars.sd)
    }
    return(list(params = params, xstars = xstars))
  })
}

