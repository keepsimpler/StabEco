#' @title estimate semicircle eigenvalue of a bipartite regular graph using alon-method, plus-twins-method, or an empirical method
#' @param n, node number
#' @param km, node degree
estimate_semicircle_bipartite_regular <- function(n, km) {
  threshold1 = (sqrt(2 * n - 3) + 1) / 2
  threshold2 = n * (n - 2) / (2 * n + 12)
  alon <- 2 * sqrt(km - 1)
  plustwins <- 2 * sqrt(km * (n - 2 * km) / (n - 2))
  plusone <- plustwins - km / (n / 2 - 1)
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
       plustwins = plustwins, plusone = plusone, empirical = empirical)
}

# the two possible equilibrium abundances
get_xstars = function(r, s, c, kc, m, km, h) {
  A = (s + kc * c) * h * km * m
  B = (s + kc * c) - km * m - r * h * km * m
  C = - r
  delta = B^2 - 4 * A * C
  if (delta >= 0) {
    X1 = (-B + sqrt(B^2 - 4 * A * C)) / (2 * A)
    X2 = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    return(list(X1 = X1, X2 = X2))
  }
  else {
    warning('delta = B^2 - 4AC < 0')
    return(list(X1 = NaN, X2 = NaN))
  }
}

# the ratio between mutualistic strength and competitive strength
get_rho = function(s, c, kc, m, km) {
  return(km * m / (s + kc * c))
}

# the real semicircle eigenvalue of the mutualistic adjacency matrix
get_m_semicircle_real = function(graphm) {
  semicircle_real <- sort(eigen(graphm)$values, decreasing = TRUE)[2]
}

# the estimated semicircle eigenvalue of the mutualistic adjacency matrix
get_m_semicircle_est = function(n, km) {
  semicircle_est <- estimate_semicircle_bipartite_regular(n, km)$est
}

# the real ratio between mutual spectral gap and competitive spectral gap
get_Delta_real = function(c, kc, m, km, graphm) {
  semicircle_real <- get_m_semicircle_real(graphm)
  Delta_real <- ((km - semicircle_real) * m) / ((kc + 1) * c)
  Delta_real
}

# the estimated ratio between mutual spectral gap and competitive spectral gap
get_Delta_est = function(c, kc, m, km, n) {
  semicircle_est <- get_m_semicircle_est(n, km)
  Delta_est <- ((km - semicircle_est) * m) / ((kc + 1) * c)
  Delta_est
}

# the dot eigenvalue
get_dot = function(r, rho, h) {
  focus <- sqrt((rho + r * h * rho + 1)^2 - 4 * rho)  # B^2-4AC
  numerator <- (rho + r * h * rho - 1 + focus) * focus
  denominator <- h * rho * (rho + r * h * rho + 1 + focus)
  - numerator / denominator
}

# the minimal value of [r] where critical transitions happen
get_rmin = function(s, c, kc, m, km, h) {
  return(-(sqrt(km * m) - sqrt(s + kc * c))^2 / (h * km * m))
}

# the effective mutualistic strength caused by h > 0
get_m_tilde = function(r, s, c, kc, m, km, h) {
  x <- get_xstars(r, s, c, kc, m, km, h)$X1
  m_tilde <- m  / (1 + h * km * m * x)^2
  return(m_tilde)
}

# the dot eigenvalue and the real and estimated semicircle eigenvalue of the Shadow Jacobian matrix and the Jacobian matrix
get_dot_semicircle = function(r, s, c, kc, m, km, h, n, graphm) {
  m_semicircle_real <- get_m_semicircle_real(graphm)
  m_semicircle_est <- get_m_semicircle_est(n, km)
  m_tilde <- get_m_tilde(r, s, c, kc, m, km, h)
  x <- get_xstars(r, s, c, kc, m, km, h)$X1 # the stable equilibrium
  # the dot eigenvalue of shadow-Jacobian
  Jshadow_dot <- km * m_tilde - kc * c - s
  # the semicircle eigenvalue of shadow-Jacobian
  Jshadow_semicircle_real <- m_semicircle_real * m_tilde + c - s
  Jshadow_semicircle_est <- m_semicircle_est * m_tilde + c - s
  # the dot eigenvalue of Jacobian
  J_dot <- Jshadow_dot * x
  # the semicircle eigenvalue of Jacobian
  J_semicircle_real <- Jshadow_semicircle_real * x
  J_semicircle_est <- Jshadow_semicircle_est * x
  list(Jshadow_dot = Jshadow_dot, Jshadow_semicircle_real = Jshadow_semicircle_real, Jshadow_semicircle_est = Jshadow_semicircle_est, J_dot = J_dot, J_semicircle_real = J_semicircle_real, J_semicircle_est = J_semicircle_est)  # , Vc = Vc, Vs = Vs, asyn = asyn
}

#' @title calculate Jocobian matrix in equilibrium of model \code{model_lv2_cm}
#' @param params, parameters of model
#' @param xstars, species densities in equilibrium
#' @return Jacobian matrix
get_jacobian_from_params_xstars <- function(params, xstars) {
  with(params, {
    # Competitive part
    Jc <- - C
    # Mutualistic part
    m_tilde <- 1 / (1 + h * M %*% xstars)^2
    Jm <- diag(c(m_tilde)) %*% M
    #Jm <-  diag(xstars / ((1  + h * rowSums(M %*% diag(xstars)))^2)) %*% M
    Jshadow <- Jc + Jm
    J <- diag(xstars) %*% (Jc + Jm)
    return(list(Jc = Jc, M = M, m_tilde = m_tilde, M_tilde = Jm, Jshadow = Jshadow, J = J))
  })
}

get_stability_from_params_xstars <- function(params, xstars) {
  n <- length(xstars)
  r <- unique(params$r)
  xstars_mean = mean(xstars)
  xstars_sd = sd(xstars)
  persistence = length(xstars[xstars > 0])
  Jacobians <- get_jacobian_from_params_xstars(params, xstars)
  M <- Jacobians$M
  M_tilde <- Jacobians$M_tilde
  Jshadow <- Jacobians$Jshadow
  J <- Jacobians$J
  m_tilde <- Jacobians$m_tilde
  
  eigenvalues <- Re(eigen(M)$values)
  M_eigenvalues <- sort(eigenvalues, decreasing = TRUE)
  M_lambda1 <- M_eigenvalues[1]
  M_lambda2 <- M_eigenvalues[2]
  
  eigenvalues <- Re(eigen(M_tilde)$values)
  M_tilde_eigenvalues <- sort(eigenvalues, decreasing = TRUE)
  M_tilde_lambda1 <- M_tilde_eigenvalues[1]
  M_tilde_lambda2 <- M_tilde_eigenvalues[2]
  M_tilde_dot <- sum(M_tilde) / n
  
  eigenvalues <- Re(eigen(Jshadow)$values)
  Jshadow_eigenvalues <- sort(eigenvalues, decreasing = TRUE)
  Jshadow_lambda1 <- Jshadow_eigenvalues[1]
  Jshadow_lambda2 <- Jshadow_eigenvalues[2]
  Jshadow_dot <- sum(Jshadow) / n
  
  eigenvalues <- Re(eigen(J)$values)
  J_eigenvalues <- sort(eigenvalues, decreasing = TRUE)
  J_lambda1 <- J_eigenvalues[1]
  J_lambda2 <- J_eigenvalues[2]
  J_dot <- sum(J) / n
  
  m_tilde_mean <- mean(m_tilde)
  m_tilde_sd <- sd(m_tilde)
  
  c(r = r, xstars = xstars, xstars_mean = xstars_mean, xstars_min = min(xstars), xstars_max = max(xstars), xstars_sd = xstars_sd, persistence = persistence, M_lambda1 = M_lambda1, M_lambda2 = M_lambda2, M_tilde_lambda1 = M_tilde_lambda1, M_tilde_lambda2 = M_tilde_lambda2, M_tilde_dot = M_tilde_dot, Jshadow_lambda1 = Jshadow_lambda1, Jshadow_lambda2 = Jshadow_lambda2, Jshadow_dot = Jshadow_dot, J_lambda1 = J_lambda1, J_lambda2 = J_lambda2, J_dot = J_dot, m_tilde = m_tilde, m_tilde_mean = m_tilde_mean, m_tilde_sd = m_tilde_sd) # id = id, m = m, h = h, graphs_index = graphs_index, 
}
