library(R6)
library(igraph)

#' @title graph represented by two blocks matrix
#' @examples 
#' TwoblocksGraph$new(type = 'two_blocks_regular', n1 = 50, n2 = 50, k = 5)
TwoblocksGraph <- R6Class('TwoblocksGraph', 
  inherit = Graph,
  public = list(
    n1 = NULL,
    n2 = NULL,
    k = NULL,
    initialize = function(type = 'two_blocks_regular', n1, n2, k) {
      graph = matrix(0, nrow = n1 + n2, ncol = n1 + n2)
      graph[1:n1, 1:n1] = as.matrix(as_adjacency_matrix(sample_k_regular(n1, k)))
      graph[(n1+1):(n1+n2), (n1+1):(n1+n2)] = as.matrix(as_adjacency_matrix(sample_k_regular(n2, k)))
      self$type <- type
      self$n1 = n1
      self$n2 = n2
      self$k = k
      self$n <- self$n1 + self$n2
      private$G <- graph_from_adjacency_matrix(graph)
    },
    get_graph = function() {
      graph <- as.matrix(as_adjacency_matrix(private$G))
      return(graph)
    },
    plot = function() {
      plot(private$G) # , layout = layout_as_bipartite
    }
  ))