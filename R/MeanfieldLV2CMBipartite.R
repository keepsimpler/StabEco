library(R6)

#' @title The mean field model for LV2-CM-Bipartite, i.e. Lotka-Volterra (LV) Equations of Holling type II for a community with two groups with Competitive intra-group interactions and Mutualistic inter-group interactions
#' @note Restrictions:
#' \describe{
#' \item{1.}{all species have the same number of mutualistic and competitive interactions}
#' \item{2.}{all species have the same strength for intrinsic growth rates, self-regulations, mutualistic and competitive interactions}
#' \item{3.}{two groups have the same number of species, i.e. n1 == n2}
#' }
#'  

#' @title get equilibrium values of state variables
get_xstars = function(r, s, c, kc, m, km, h) {
  if (h == 0) {
    if (s + kc * c > km * m) {  # rho < 1
      X1 = r / (s + kc * c - km * m)
      return(list(X1 = X1, X2 = NaN))
    }
    else {  # rho > 1
      warning('h = 0, rho > 1, with infinite equilibrium value!')
      return(list(X1 = Inf, X2 = NaN))
    }
  }
  # h > 0
  A = (s + kc * c) * h * km * m
  B = (s + kc * c) - km * m - r * h * km * m
  C = - r
  delta = B^2 - 4 * A * C
  if (delta >= 0) {    # >= rmin
    X1 = (-B + sqrt(B^2 - 4 * A * C)) / (2 * A)
    X2 = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    return(list(X1 = X1, X2 = X2))
  }
  else {  # < rmin
    warning('delta = B^2 - 4AC < 0')
    return(list(X1 = NaN, X2 = NaN))
  }
}

#' @title estimate semicircle eigenvalue of a bipartite regular graph using alon-method, plus-twins-method, or an empirical method
#' @param n, node number
#' @param km, node degree
estimate_semicircle_bipartite_regular <- function(n, km) {
  threshold1 = (sqrt(2 * n - 3) + 1) / 2
  threshold2 = n * (n - 2) / (2 * n + 12)
  alon <- 2 * sqrt(km - 1)
  plustwins <- 2 * sqrt(km * (n - 2 * km) / (n - 2))
  #plusone <- plustwins - km / (n / 2 - 1)
  empirical <- n / 2 - km
  if (km >= 2 && km < threshold1) { # using alon-method
    est <- alon
  }
  else if (km >= threshold1 && km < threshold2) { # using plus-twins method
    est <- plustwins
  }
  else if (km >= threshold2 && km <= n / 2 - 1) {
    est <- empirical
  }
  list(est = est, threshold1 = threshold1, threshold2 = threshold2, alon = alon,
       plustwins = plustwins, empirical = empirical) #, plusone = plusone
}

#' @title get the estimated semicircle eigenvalue of the mutualistic adjacency matrix
get_m_semicircle_est = function(n, km) {
  semicircle_est <- estimate_semicircle_bipartite_regular(n, km)$est
}

#' @title get the real semicircle eigenvalue of the mutualistic adjacency matrix
get_m_semicircle_real = function(graphm) {
  graphm = graphm$get_graph()
  semicircle_real <- sort(eigen(graphm)$values, decreasing = TRUE)[2]
}

#' @title get the effective mutualistic strength caused by h > 0
get_m_tilde = function(r, s, c, kc, m, km, h) {
  x <- get_xstars(r, s, c, kc, m, km, h)$X1
  stopifnot(x >= 0)
  m_tilde <- m  / (1 + h * km * m * x)^2
  return(m_tilde)
}

#' @title get the ratio between mutualistic strength and competitive strength
get_rho = function(s, c, kc, m, km) {
  return(km * m / (s + kc * c))
}

#' @title get the dot eigenvalue of the Jacobian matrix
get_dot = function(r, rho, h) {
  focus <- sqrt((rho + r * h * rho + 1)^2 - 4 * rho)  # B^2-4AC
  numerator <- (rho + r * h * rho - 1 + focus) * focus
  denominator <- h * rho * (rho + r * h * rho + 1 + focus)
  - numerator / denominator
}

#' @title get the estimated semicircle eigenvalue of the Jacobian matrix
get_semicircle_est = function(n, r, s, c, kc, m, km, h) {
  m_semicircle_est <- get_m_semicircle_est(n, km)
  x <- get_xstars(r, s, c, kc, m, km, h)$X1 # the stable equilibrium
  stopifnot(x >= 0)
  m_tilde <- get_m_tilde(r, s, c, kc, m, km, h)
  Jshadow_semicircle_est <- m_semicircle_est * m_tilde + c - s
  J_semicircle_est <- Jshadow_semicircle_est * x
}

#' @title get the real semicircle eigenvalue of the Jacobian matrix
get_semicircle_real = function(n, r, s, c, kc, m, km, h, graphm) {
  m_semicircle_real <- get_m_semicircle_real(graphm)
  x <- get_xstars(r, s, c, kc, m, km, h)$X1 # the stable equilibrium
  stopifnot(x >= 0)
  m_tilde <- get_m_tilde(r, s, c, kc, m, km, h)
  Jshadow_semicircle_real <- m_semicircle_real * m_tilde + c - s
  J_semicircle_real <- Jshadow_semicircle_real * x
}

#' @title get the value of \code{r} where the dot eigenvalue equal 0 and critical transitions can happen
get_rmin = function(s, c, kc, m, km, h) {
  return(-(sqrt(km * m) - sqrt(s + kc * c))^2 / (h * km * m))
}

#' @title get the value of \code{r} where the semicircle eigenvalue equal 0 and critical transitions can happen
get_r_semicircle0 <- function(s, c, kc, m, km, h, n, graphm) {
  rho <- get_rho(s, c, kc, m, km)
  m_semicircle_real <- get_m_semicircle_real(graphm)
  semi_ratio_real <- (m_semicircle_real * m) / (s - c)
  r_real <- 1/h * ( 1/rho * (semi_ratio_real^0.5 - 1) +  semi_ratio_real^-0.5 - 1 )
  r_real
}

#' @title get the real ratio between mutual spectral gap and competitive spectral gap
get_Delta_real = function(c, kc, m, km, graphm) {
  m_semicircle_real <- get_m_semicircle_real(graphm)
  Delta_real <- ((km - m_semicircle_real) * m) / ((kc + 1) * c)
  Delta_real
}

#' @title  the estimated ratio between mutual spectral gap and competitive spectral gap
get_Delta_est = function(c, kc, m, km, n) {
  m_semicircle_est <- get_m_semicircle_est(n, km)
  Delta_est <- ((km - m_semicircle_est) * m) / ((kc + 1) * c)
  Delta_est
}


MeanfieldLV2CMBipartite <- R6Class('MeanfieldLV2CMBipartite',
    public = list(
      n1 = NULL,
      n2 = NULL,
      n = NULL,
      r = NULL,
      s = NULL,
      kc = NULL,
      c = NULL,
      km = NULL,
      m = NULL,
      h = NULL,
      graphc = NULL,
      graphm = NULL,
      initialize = function(n, r, s, kc = NULL, c, km, m, h, graphc, graphm) {
        stopifnot(graphc$type == 'two_blocks_regular' && graphc$n == n &&
                    graphc$n1 == n / 2 && graphc$n2 == n / 2 &&
                    graphm$type == 'bipartite_regular' && graphm$n == n &&
                    graphm$n1 == n / 2 && graphm$n2 == n / 2)
        self$n <- n
        self$n1 <- n / 2
        self$n2 <- n / 2
        self$s <- s
        if (is.null(kc))
          self$kc <- self$n1 - 1
        else
          self$kc <- kc
        self$c <- c
        self$km <- km
        self$m <- m
        self$h <- h
        self$graphc <- graphc
        self$graphm <- graphm
      },
      get_xstars = get_xstars,
      estimate_semicircle_bipartite_regular = estimate_semicircle_bipartite_regular
      
    ))