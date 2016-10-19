library(R6)
library(igraph)

#' @title Graph of many different types
#' @field type graph types:
#' \describe{
#' \item{bipartite_regular}{regular bipartite graph}
#' \item{two_blocks_regular}{regular graph with two blocks}
#' }
#' @examples 
#' Graph$new(type = 'bipartite_regular', n1 = 50, k = 5, directed = T, is_adj = T)
#' Graph$new(type = 'two_blocks_regular', n1 = 50, n2 = 50, k = 5)
Graph <- R6Class('Graph',
  public = list(
    type = NULL,
    initialize = function(type = c('bipartite_regular', # n1 k directed is_adj
                                   'two_blocks_regular' # n1 n2 k
                                   ), 
                          n1 = NULL, n2 = NULL, k = NULL, 
                          directed=TRUE, is_adj = TRUE, ...) {
      type <- match.arg(type)
      switch (type,
        bipartite_regular  = {
          G = sample_k_regular(n1, k, directed = directed)
          graph <- as.matrix(as_adjacency_matrix(G))
          # shuffle the rows and cols, to simulate a bipartite regular graph,
          # rather than a unipartite regular graph which have 0 value diagonal elements
          # this permutation does not change the eigenvalues
          # graph <- graph[sample.int(s[1]), sample.int(s[1])]
          if (is_adj == TRUE)
            graph <- inc_to_adj(graph)
        },
        two_blocks_regular = {
          graph = matrix(0, nrow = n1 + n2, ncol = n1 + n2)
          graph[1:n1, 1:n1] = as.matrix(as_adjacency_matrix(sample_k_regular(n1, k)))
          graph[(n1+1):(n1+n2), (n1+1):(n1+n2)] = as.matrix(as_adjacency_matrix(sample_k_regular(n2, k)))
        }
      )
      private$graph <- graph
      self$type <- type
    },
    get_graph = function() {
      private$graph
    }
  ),
  private = list(
    graph = NULL
  ))

