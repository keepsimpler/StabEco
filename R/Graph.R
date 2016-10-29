library(R6)
library(igraph)

#' @title The root class of all graph objects
#' @field type graph type
#' @field n number of nodes
#' @field graph inner matrix of the graph
Graph <- R6Class('Graph',
  public = list(
    type = NULL,
    n = NULL,
    initialize = function() {
      cat('Initialize Graph object.\n')
    },
    get_graph = function() {
      private$graph
    }
  ),
  private = list(
    graph = NULL
  ))

