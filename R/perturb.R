library(R6)

#' @title Perturbations on model \code{model_lv2_cm} through by parameters [params]
#' @param params parameters
#' @param type perturbation types:
#' \describe{
#'   \item{growthrate_all}{increasing/decreasing the intrinsic growth rates of ALL species}
#'   \item{growthrate_part}{increasing/decreasing the intrinsic growth rates of a part of species}
#'   \item{mutualistic_strength}{increasing/decreasing strengths of mutualistic interactions}
#'   \item{competitive_strength}{increasing/decreasing strengths of competitive interactions}
#'   \item{primary_extinct}{remove one species}
#' }
#' @param coeffs list of coefficients of perturbations:
#' \describe{
#'   \item{r.delta.mu, r.delta.sd}{stepwise of inc/dec of intrinsic growth rate}
#'   \item{m.delta.mu, m.delta.sd}{stepwise of inc/dec of mutualistic strength}
#'   \item{extinct_species}{a (index) vector of the removed species}
#' }
perturb_lv2_cm <- function(params, type = c('growthrate_all', 'mutualistic_strength', 'primary_extinct'), coeffs) {
  type <- match.arg(type)
  with(coeffs, {
    if (type == 'growthrate_all') {
      params$r = params$r - runif2(length(params$r), r.delta.mu, r.delta.sd)
    }
    else if (type == 'mutualistic_strength') {
      edges = length(params$M[params$M > 0])
      params$M[params$M > 0] = params$M[params$M > 0] + runif(edges, min = m.delta.mu - m.delta.sd, max = m.delta.mu + m.delta.sd)
    }
    else if (type == 'primary_extinct') {
      params$r = params$r[- extinct_species]
      params$C = params$C[- extinct_species, - extinct_species]
      params$M = params$M[- extinct_species, - extinct_species]
      params$h = params$h[- extinct_species]
    }
    return(params)
  })
}

#' @title Perturbations on complex systems through by state variables at equilibrium [xstars]
#' @param xstars state variables at equilibrium
#' @param type perturbation types:
#' \describe{
#'   \item{immigration}{simulating a immigration factor to (re)establish extinct species}
#'   \item{random}{randomly perturb species abundances at equilibrium}
#' }
#' @param coeffs list of coefficients of perturbation:
#' \describe{
#'   \item{immigration.mu, immigration.sd}{mean and sd of the immigration}
#'   \item{xstars.sd}{sd of random perturbations}
#' }
perturb_xstars <- function( xstars, type = c('growthrate_all'), coeffs) {
  type <- match.arg(type)
  with(coeffs, {
    if (type == 'immigration')
      xstars[xstars == 0] <- xstars[xstars == 0] +
        runif2(length(xstars[xstars == 0]), immigration.mu, immigration.sd)
    # random perturb xstars
    if (type == 'random') {
      #set.seed(123) # fix the random number to check its influence
      xstars <- xstars * runif2(length(xstars), 1, xstars.sd)
    }
    return(xstars = xstars)
  })
}

