

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

get_stability_from_params_xstars <- function(params, xstars, flag = 'press') {
  n <- length(xstars)
  r <- mean(params$r)
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
  if (flag == 'auto') {
    ret <- c(r = r, xstars = xstars, xstars_mean = xstars_mean, xstars_min = min(xstars), xstars_max = max(xstars), xstars_sd = xstars_sd, persistence = persistence, M_lambda1 = M_lambda1, M_lambda2 = M_lambda2, M_tilde_lambda1 = M_tilde_lambda1, M_tilde_lambda2 = M_tilde_lambda2, M_tilde_dot = M_tilde_dot, Jshadow_lambda1 = Jshadow_lambda1, Jshadow_lambda2 = Jshadow_lambda2, Jshadow_dot = Jshadow_dot, J_lambda1 = J_lambda1, J_lambda2 = J_lambda2, J_dot = J_dot, m_tilde = m_tilde, m_tilde_mean = m_tilde_mean, m_tilde_sd = m_tilde_sd) # id = id, m = m, h = h, graphs_index = graphs_index, 
  }
  else if (flag == 'press') {
    ret <- c(r = r, xstars_mean = xstars_mean, xstars_min = min(xstars), xstars_max = max(xstars), xstars_sd = xstars_sd, persistence = persistence, M_lambda1 = M_lambda1, M_lambda2 = M_lambda2, M_tilde_lambda1 = M_tilde_lambda1, M_tilde_lambda2 = M_tilde_lambda2, M_tilde_dot = M_tilde_dot, Jshadow_lambda1 = Jshadow_lambda1, Jshadow_lambda2 = Jshadow_lambda2, Jshadow_dot = Jshadow_dot, J_lambda1 = J_lambda1, J_lambda2 = J_lambda2, J_dot = J_dot, m_tilde_mean = m_tilde_mean, m_tilde_sd = m_tilde_sd) 
  }
  ret
}
