library("amap")
library("mco")
library("cluster")
source("nsga2.r")

helper.normalize <- function(data) {
  return( (data - min(data)) / (max(data) - min(data)) )
}

evaluator.multiobjective.clustering <- function( results, dmatrix ) {
  
  num_clusters <- ncol(results$population) - 4
  solution_count <- length(results$clustering)
  
  silhouette_indices = double(solution_count)
  for (c_index in 1:solution_count ) {
    clustering <- results$clustering[[c_index]]
    
    sil <- silhouette( clustering, dmatrix = as.matrix(dmatrix) )
    silhouette_indices[c_index] <- summary(sil)$avg.width
  }
  
  pareto_points <- as.matrix( results$population[ , (num_clusters+1):(num_clusters+2) ] )
  hypervolume <- dominatedHypervolume(pareto_points)
  
  pareto_points[, 1] <- helper.normalize(pareto_points[, 1])
  pareto_points[, 2] <- helper.normalize(pareto_points[, 2])
  n_hypervolume <- dominatedHypervolume(pareto_points, c(1, 1))
  
  metrics <- list(
    silhoutte = silhouette_indices,
    hypervolume = hypervolume,
    normalized_hypervolume = n_hypervolume
  )
  return( metrics )
}

evaluator.biological.significance <- function() {
}

evaluator.metaheuristics <- function(metaheuristic, meta_params, run_evaluator = evaluator.multiobjective.clustering, runs = 10) {
  results <- sapply( 1:runs, function(n) run_evaluator(do.call(metaheuristic, meta_params)) )
  print( results )
  return( colMeans( results ) )
  #return( rowMeans( results ) )
}

evaluator.friedman <- function() {
}
