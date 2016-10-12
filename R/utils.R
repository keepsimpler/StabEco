#' @title ggplot a multivariate timeseries
#' @param timeseries a numeric matrix representing a multivariate timeseries, the first colnum is the index data
#' @param title, is_legend
ggplot_timeseries <- function(timeseries, title = NULL, xlab = NULL, ylab = NULL, is_legend = FALSE, size = 0.5) {
  timeseries_length <- dim(timeseries)[1]
  variable_num <- dim(timeseries)[2]
  timeseries <- data.frame(timeseries)
  colnames(timeseries)[1] <- 'time'  # the name of the index column
  require(reshape2)
  # reshape to a long format
  timeseries.long <- melt(data = timeseries, id.vars = 'time', variable.name = "variable", value.name = "value")
  p <- ggplot() + #, lty = variable
    geom_line(data = timeseries.long, aes(x = time,y = value, colour = variable), size = size) +
    theme_bw() +
    labs(title = title, x = xlab, y = ylab) +
    theme(legend.position = ifelse(test = is_legend,"right", "none"))
  p
}

