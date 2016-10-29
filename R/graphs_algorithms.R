#' @title Generate bipartite graphs with different degree heterogeneity by rewiring links (or Nestedness optimization algrithm).
#' @description 
#' 1) start from a regular bipartite graph.
#' 2) repeat rewiring a link to a new species who has more neighbors than the old species
#' until failed for enough times.
#' 3) return to step 2)
#' @examples 
#' biGraph <- BiGraph$new(type = 'bipartite_regular', n1 = 10, k = 4, is_adj = FALSE)
#' graphs <- bigraphs_rewiring(biGraph$get_graph(is_adj = FALSE))
bigraphs_rewiring <- function(graph) {
  graphs.rewiring = list()  # initialize the generated graphs by rewiring links

  n1 <- dim(graph)[1]
  n2 <- dim(graph)[2]
  rownames(graph) <- 1:n1
  colnames(graph) <- (n1 + 1):(n1 + n2)
  count = 0
  graphs.rewiring[[1]] <- graph # list(count = count, graph = graph)
  repeat {
    count = count + 1
    shouldcontinue = FALSE
    ## rewiring one link to a random node which has more neighbors
    ## if tring [ntry]*5 times, and still fail to rewire, then [shouldcontinue] is false.
    for (i in 1:5) {
      tmp = bigraphs_rewiring_onestep(graph, ntry = 1000)
      if (tmp$flag == TRUE) {  # the rewiring is success
        shouldcontinue = TRUE
        break
      }
      else {
        graph = tmp$graph
      }
    }
    if (!shouldcontinue) break
    graph = tmp$graph  # the new graph which has more degree heterogeneity
    graphs.rewiring[[length(graphs.rewiring) + 1]] <-  graph
    print(count)
  }
  
  # recede the positions of species
  graphs <- lapply(graphs.rewiring, function(graph) {
    #seq.row <- rownames(graph)[order(as.numeric(rownames(graph)))]
    #seq.col <- colnames(graph)[order(as.numeric(colnames(graph)))]
    seq.row <- as.character(1:n1)
    seq.col <- as.character((n1+1):(n1+n2))
    require(bipartite)
    sortweb(graph, sort.order = 'seq', sequence = list(seq.higher = seq.col, seq.lower = seq.row))
  })
  
  # compare before rewired links, and after rewired links
  # there are exactly TWO different links before and after one rewiring step
  # This is to check if the rewiring algo. is correct
  #graphs_rewired_links <- sapply(2:length(graphs), function(i) {
  #  which(graphs[[i]] - graphs[[i-1]] != 0)  # , arr.ind = T
  #})
  #dim(graphs_rewired_links) # check the rewired links
  
  graphs
}

###############################################################################
#' @title rewiring a link to some node with more links. (richer get richer)
#'
#' @param B incidence matrix of bipartite network, rows and cols represent two groups of nodes/species
#' @param connected, if the new graph should be connected?
#' @param ntry, how many to try?
#' @return the incidence matrix whose nestedness has been optimized by rewiring a link.
#' @details .
#' @import bipartite
bigraphs_rewiring_onestep <- function(B, connected = TRUE, ntry = 100) {
  require(bipartite) # sortweb function
  require(igraph)
  B = bipartite::sortweb(B)  # sort rows and cols descending, ensure the chosen species later has more interactions
  NumP <- dim(B)[1]
  NumA <- dim(B)[2]
  flag = FALSE # is the rewiring succeed, or the max tried times approach but the rewiring still not succeed
  ## random choose another species (plant or animal with equal probability),
  ## and rewire the link to the new selectd species
  for (i in 1:ntry) {
    flag1 = FALSE  #  if this rewiring has succeed?
    flag2 = FALSE  #  if the new graph is connected?
    B2 = B  # copy the original graph, in order to try more than one times
    ## pick one interaction between two random species
    repeat {
      row1 <- sample(1:NumP, 1)
      col1 <- sample(1:NumA, 1)
      if (B2[row1, col1] != 0) break
    }
    if (runif(1) < 0.5) {  # choose another plant  #NumP/(NumP + NumA)
      row2 =  sample(1:row1, 1)  # choose random plant with more interactions which is ensured by [sortweb]
      # Three exceptions: 1. the new chosen species [row2] is same with the old species [row1]
      # 2. the new chosen species [row2] already has interaction with [col1]
      # 3. the old species [row1] has only one interaction.
      if (row2 < row1 && B2[row2, col1] == 0 && sum(B2[row1,]) > 1) {
        B2[row2, col1] = B2[row1, col1]
        B2[row1, col1] = 0
        flag1 = TRUE  # the link has been rewired to a new plant
      }
    }
    else {  # choose another animal
      col2 =  sample(1:col1, 1)
      if (col2 < col1 && B2[row1, col2] == 0 && sum(B2[,col1]) > 1 ) {
        B2[row1, col2] = B2[row1, col1]
        B2[row1, col1] = 0
        flag1 = TRUE  # the link has been rewired to a new animal
      }
    }
    ## if the new graph is connected, [flag2] is TRUE
    G = igraph::graph_from_incidence_matrix(B2, add.names = NA)
    if (igraph::is.connected(G)) flag2 = TRUE
    
    ## if the rewiring is success, and (the new graph is connected or that is not required)
    if (flag1 && (flag2 || !connected)) {
      flag = TRUE
      break;
    }
  }
  if (flag == FALSE) {  # if failed, return the original matrix
    res = list(graph = B, flag = flag, tried = i)
    warning(paste(ntry, 'times has been tried, but the rewiring still donot succeed.'))
  }
  else {
    res = list(graph = B2, flag = flag, tried = i)
  }
  #sortweb(B)  # sort rows and cols descending for the next link rewiring.
  res
}

#' @title get measurements of a bipartite graph
#' @param graph the incident matrix of a bipartite graph
#' @return a list of measurements:
#' two measurements of degree heterogeneity:
#'  1. degree variance
#'  2. degree entropy
#'  two measurements of eigenvalues:
#'  1. the largest eigenvalue
#'  2. the second eigenvalue
#'  and the degrees sequence
get_degree_heterogeneity <- function(graph) {
  degrees <- c(rowSums(graph), colSums(graph))
  degree_avg <- mean(degrees) # average degree
  variance <- sum(degrees^2) / sum(degrees)^2
  entropy <- sum(degrees * log(degrees), na.rm = TRUE)
  variance_avg <- sum((degrees / degree_avg)^2)
  entropy_avg <- sum((degrees / degree_avg) * log((degrees / degree_avg)), na.rm = TRUE)
  graphm <- inc_to_adj(graph)
  eigenvalues <- sort(eigen(graphm)$values, decreasing = TRUE)
  c(variance = variance, entropy = entropy, variance_avg = variance_avg, entropy_avg = entropy_avg, lambda1 = eigenvalues[1], lambda2 = eigenvalues[2], gap = eigenvalues[1] - eigenvalues[2]) # degrees = degrees, 
}

#' @title generate a sequence of bipartite graphs with increasing hetero
#' @param n1 km
#' @return multi objects: a list of graphs, a dataframe of hetero measurements, and a dataframe of degrees
#' @examples 
#' list[graphs, graphs_hetero, graphs_degrees] <- get_graphs(n1 = 50, km = 5)
get_graphs <- function(n1, km) {
  biGraph <- BiGraph$new(type = 'bipartite_regular', n1 = n1, k = km, is_adj = FALSE)
  graphs <- bigraphs_rewiring(biGraph$get_graph(is_adj = FALSE))
  # compute degree heterogeneitys of graphs
  graphs_hetero <- ldply(graphs, function(graph) {
    get_degree_heterogeneity(graph)
  })
  graphs_hetero$graphs_index <- 1:nrow(graphs_hetero)
  graphs_degrees <- ldply(graphs, function(graph) {
    degrees <- c(rowSums(graph), colSums(graph))
    degrees
  })
  graphs_degrees$graphs_index <- 1:nrow(graphs_degrees)
  list(graphs = graphs, graphs_hetero = graphs_hetero, graphs_degrees = graphs_degrees)
}

get_graphs_sf <- function(n1, km) {
  alphas = seq(from = 1, to = 0.01, length.out = 100)
  gammas = 1/ alphas + 1 #seq(from = 2, to = 10, by = 0.1) # exponent 
  graphs <- lapply(gammas, function(gamma) {
    G = sample_fitness_pl(n1, n1 * km, exponent.out = gamma, exponent.in = gamma)
    graph <- as.matrix(as_adjacency_matrix(G))
  })
  sapply(graphs, sum)
  # compute degree heterogeneitys of graphs
  graphs_hetero <- ldply(graphs, function(graph) {
    get_degree_heterogeneity(graph)
  })
  graphs_hetero$graphs_index <- 1:nrow(graphs_hetero)
  graphs_degrees <- ldply(graphs, function(graph) {
    degrees <- c(rowSums(graph), colSums(graph))
    degrees
  })
  graphs_degrees$graphs_index <- 1:nrow(graphs_degrees)
  list(graphs = graphs, graphs_hetero = graphs_hetero, graphs_degrees = graphs_degrees)
}
