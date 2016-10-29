#' @title numerical differentiation algorithm to estimate the derivative of a function based on numerical data
#' @references https://en.wikipedia.org/wiki/Numerical_differentiation
numeric_derivative <- function(data, h) {
  (data[-c(1,2)] - data[-c(length(data), length(data) - 1)]) / (2 * h)
}

#' @title get covariances within moving windows along a multivariate timeseries
#' @param timeseries a numeric matrix of multivariate timeseries
#' @param winsize window size expressed as percentage of the timeseries length (must be numeric between 0 and 100). Default is 50.
#' @param type c('cov', 'cor', 'relative')
#' @return a list of covariances measurements
run_stats_multivariate <- function(timeseries, winsize = 50,
                                   type = c('cov', 'cor', 'relative')) {
  type <- match.arg(type)
  timeseries_length <- dim(timeseries)[1]
  variable_num <- dim(timeseries)[2]
  window_size <- round(timeseries_length * winsize/100) # the absolute window size
  window_num <- timeseries_length - window_size + 1 # the number of windows
  covariances <- sapply(1:window_num, function(i) {
    data_of_window <- timeseries[i:(i + window_size - 1),]
    if (type == 'cov')
      covariance_matrix <- cov(data_of_window)
    else if (type == 'cor')
      covariance_matrix <- cor(data_of_window)
    else if (type == 'relative') {
      means <- apply(data_of_window, 2, mean)
      covariance_matrix <- cov(data_of_window)
      covariance_matrix <-
        diag(1/means) %*% covariance_matrix %*% diag(1/means)
    }
    Vc <- sum(covariance_matrix)
    Vs <- sum(sqrt(diag(covariance_matrix)))^2
    asyn <- Vc / Vs
    c(Vc = Vc, Vs = Vs, asyn = asyn)
  })
  list(Vc = covariances[1,], Vs = covariances[2,], asyn = covariances[3,])
}

#' @title preprocess a multivariate timeseries
#' @description trim all data starting from the first negative value, or replace with 0
preproc_trim_negative <- function(timeseries, trim_or_replace = c('trim', 'replace')) {
  timeseries <- data.matrix(timeseries)
  # the index where the first negative value emerges
  idx <- min(which(timeseries < 0, arr.ind = TRUE)[,1])
  trim_or_replace <- match.arg(trim_or_replace)
  if (trim_or_replace == 'trim')
    timeseries <- timeseries[1:(idx-1), ]
  else if (trim_or_replace == 'replace')
    timeseries[idx:dim(timeseries)[1], ] <- 0
  timeseries
}

#' @title visualise a matrix
image2 <- function(m)
  image(t(m[nrow(m):1,]), axes=FALSE)

#' @title Emulate ggplot2 default color palette
#' @param n number of colors
#' @references http://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

gg_color_palette_default <- function(n) {
  require(graphics)
  palette('default') # reset back to the default
  rep(palette(), ceiling(n / length(palette())))[1:n]
}

#' @title another form of uniform distribution between [mean - sd, mean + sd]
runif2 <- function(n, mean, sd) {
  runif(n) * 2 * sd + (mean - sd)
}

#' @title transfer an incidence matrix to an adjacency matrix
#' @param inc, an incidency matrix
#' @return adj, an adiacency matrix
inc_to_adj <- function(inc){
  p <- dim(inc)[1]  # number of Plants
  a <- dim(inc)[2]  # number of Animals
  s <- p + a  # number of all Species
  adj <- matrix(0, s, s)  # initialize the adjacency matrix as a zero-matrix
  adj[1:p, (p + 1):s] <- inc  # the upper right sub-matrix is the incidence matrix
  adj <- adj + t(adj)  # the lower left sub-matrix is transpose of the incidence matrix
  return(adj)
}
adj_to_inc <- function(adj, n1, n2) {
  stopifnot(dim(adj)[1] == n1 + n2)
  inc <- adj[1:n1, (n1 + 1):(n1 + n2)]
}

