library(R6)
library(igraph)

#' @title The root class of all graph objects
#' @field type graph type
#' @field n number of nodes
#' @field G inner represent of the graph, using \code{igraph} package
Graph <- R6Class('Graph',
  public = list(
    type = NULL,
    n = NULL,
    initialize = function() {
      cat('Initialize Graph object.\n')
    },
    get_graph = function() {
    },
    is_connected = function(mode = c('weak', 'strong')) {
      is_connected(private$G, mode)
    }
  ),
  private = list(
    G = NULL
  ))

