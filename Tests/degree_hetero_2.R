# correlation among measurements of degree heterogeneity
n1 = 50
km = 6
alpha_min = 0
alpha_max = 4
by = 0.04
alphas = seq(from = alpha_min, to = alpha_max, by = by)
gammas = 1/ alphas + 1 #seq(from = 2, to = 10, by = 0.1) # exponent 

library(bipartite)
graphs_hetero <- ddply(data.frame(idx = 1:5), .variables = c('idx'), function(i) {
  graphs <- get_graphs_sf(n1, km, alpha_min = alpha_min, alpha_max = alpha_max, by = by, ntry = 500)
  graphs_hetero = get_graphs_hetero(graphs)
  nodfs = sapply(graphs, function(graph) nested(graph, method = c('NODF', 'NODF2')))
  nodfs = data.frame(t(nodfs))
  nodfs$variance_avg = graphs_hetero$variance_avg
  nodfs$entropy_avg = graphs_hetero$entropy_avg
  nodfs$alpha = alphas
  nodfs$gamma = gammas
  nodfs
})
graphs_hetero$beta = graphs_hetero$alpha
graphs_hetero_2 <- aggregate(cbind(NODF, variance_avg) ~beta + gamma, graphs_hetero, mean) # NODF2 ,  entropy_avg
pairs(graphs_hetero_2)

require(GGally)
require(ggplot2)

my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=loess, fill="blue", color="blue", ...) 
  p
}

g = ggpairs(graphs_hetero_2, columns = c('beta', 'NODF', 'variance_avg'), columnLabels = c('beta', 'NODF', 'H(Gm)'), lower = list(continuous = my_fn), axisLabels="internal") + theme_bw() #  diag = list(continuous = "blankDiag")
g
