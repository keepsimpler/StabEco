library(R6)
library(igraph)
library(bipartite) # use sortweb function

#' @title bipartite graph
#' @examples 
#' BiGraph$new(type = 'bipartite_regular', n1 = 50, k = 5, directed = T, is_adj = T)
BiGraph <- R6Class('BiGraph', 
  inherit = Graph,
  public = list(
    n1 = NULL,
    n2 = NULL,
    k = NULL,
    # directed = NULL,
    is_adj = NULL, # inner matrix is adjacency or incidency
    # set graph object from an incidency matrix
    set_graph_inc = function(graph_inc) {
      self$n1 = dim(graph_inc)[1]
      self$n2 = dim(graph_inc)[2]
      self$is_adj = FALSE
      private$graph = graph_inc
      self$k = 2 * sum(private$graph) / (self$n1 + self$n2)
      self$type = 'bipartite_degree_seq'
    },
    initialize = function(type = c('bipartite_regular'), n1, n2 = NULL, k,
                          directed = TRUE, is_adj = TRUE, ...) {
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
                self$n1 <- n1
                self$n2 <- n1
                self$type <- type
                self$k <- k
                self$n <- self$n1 + self$n2
                private$graph <- graph
#                self$directed <- directed
                self$is_adj <- is_adj
              }
      )
    },
    get_graph = function(is_adj = TRUE) {
      if (self$is_adj == is_adj) {
        return(private$graph)
      }
      else if (self$is_adj == TRUE && is_adj == FALSE) {
        return(adj_to_inc(private$graph, self$n1, self$n2))
      }
      else if (self$is_adj == FALSE && is_adj == TRUE) {
        return(inc_to_adj(private$graph))
      }
    },
    as_adj = function() {
      if (self$is_adj == FALSE) {
        private$graph <- inc_to_adj(private$graph)
      }
    },
    as_inc = function() {
      if (self$is_adj == FALSE) {
        private$graph <- adj_to_inc(private$graph, self$n1, self$n2)
      }
    }
  ))