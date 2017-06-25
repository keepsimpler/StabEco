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
      #cat('Graph')
      graph <- as.matrix(as_adjacency_matrix(private$G))
      return(graph)
    },
    get_hetero = function() {
      graph <- as.matrix(as_adjacency_matrix(private$G))
      degrees <- rowSums(graph)
      degree_avg <- mean(degrees) # average degree
      #variance <- sum(degrees^2) / sum(degrees)^2
      #entropy <- sum(degrees * log(degrees), na.rm = TRUE)
      variance_avg <- mean((degrees / degree_avg)^2)
      entropy_avg <- sum((degrees / degree_avg) * log((degrees / degree_avg)), na.rm = TRUE)
      eigenvalues <- sort(eigen(graph)$values, decreasing = TRUE)
      list(variance_avg = variance_avg, entropy_avg = entropy_avg, lambda1 = eigenvalues[1], lambda2 = eigenvalues[2]) # variance = variance, entropy = entropy, gap = eigenvalues[1] - eigenvalues[2]
    },
    is_connected = function(mode = c('weak', 'strong')) {
      is_connected(private$G, mode)
    }
  ),
  private = list(
    G = NULL
  ))

