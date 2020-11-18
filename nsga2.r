library(nsga2R)
library(matrixcalc)
library(pdist)

packages = c("nsga2R")

helper.order.matrix <- function(matrix) {
  matrix[ order( as.numeric(rownames(matrix)) ), order( as.numeric(colnames(matrix)) ) ]
}

# The following code is part of the example of ?nsga2R::crowdingDist4frnt
helper.pareto.ranking <- function(popSize, ranking) {
  rnkIndex <- integer(popSize)
  i <- 1
  while (i <= length(ranking)) {
    rnkIndex[ranking[[i]]] <- i
    i <- i + 1
  }
  return(rnkIndex)
}

generate.initial.population <- function(genes, population_size, num_clusters) {
  as.data.frame( t(sapply( 1:population_size, function(x) sample( genes, num_clusters, replace=F ) )) )
  #print(population_initial)
  #rm(.Random.seed, envir=globalenv())
  #return( as.matrix(population_initial) )
}

objective.functions <- function(cluster_solutions, dmatrix) {
  objective_indices <- vector()
  for(i in 1:nrow(cluster_solutions)) {
    objective_indices[i] <- fitness.medoid( cluster_solutions[i, ], dmatrix )
  }
  
  return( objective_indices )
}

fitness.medoid <- function(cluster_solution, gene_dmatrix) {
  elements <- nrow( gene_dmatrix )
  medoids <- length( cluster_solution )
  distance_to_medoids <- matrix(0, medoids, elements)
  
  for (medoid in cluster_solution) {
    for (gene in rownames(gene_dmatrix)) {
      distance_to_medoids <- (gene_dmatrix[medoid, gene]) ^ 2
    }
  }
  # check if this distances are correct, shouldn't it sum only the cluster elements?
  # i.e the minimum in each one
  
  XB_numerator = sum( distance_to_medoids )
  
  medoid_pairs = utils::combn( unlist(cluster_solution), 2 )
  
  separation = vector()
  for( col in 1:ncol(medoid_pairs) ) {
    separation[col] <- ( gene_dmatrix[ t( medoid_pairs[, col] ) ] ) ^ 2
  }
  minimum_separation <- min(separation)
  
  XB_denominator = minimum_separation
  
  return( XB_numerator / XB_denominator )
}

operator.nsga2.sorting.and.crowding <- function(population_size, population, dmatrix_expression, dmatrix_biological) {
  
  objective_exp <- objective.functions( population, dmatrix_expression )
  objective_bio <- objective.functions( population, dmatrix_biological )
  
  objectives <- cbind( objective_exp, objective_bio )
  obj_range <- ( apply(objectives, 2, max) - apply(objectives, 2, min) )
  
  non_dominated_sorted <- nsga2R::fastNonDominatedSorting( objectives )
  solutions_rank <- helper.pareto.ranking( population_size, non_dominated_sorted )
  population <- cbind( population, objectives, solutions_rank )
  
  crowding_distance <- nsga2R::crowdingDist4frnt( population, solutions_rank, obj_range )
  print(population)
  print(crowding_distance)
  
  population <- cbind( population, crowding_distance = c(apply(crowding_distance,1,sum)) )
  population<-as.data.frame(population)
  population<-population[order(population$solutions_rank, -(population$crowding_distance)), ]
  population<-as.matrix(population)
  
  return( population )
}

operator.local.search <- function(should_apply, population_size, population, dmatrix_expression, dmatrix_biological, local_search, ls_params) {
  if( should_apply == TRUE && local_search != FALSE ) {
    data <- operator.nsga2.sorting.and.crowding(population_size, population, dmatrix_expression, dmatrix_biological)
  }
  return( data )
}

operator.crossover <- function(crossover_ratio) {
  
}

operator.mutation <- function(mutation_ratio) {
  
}

operator.selection <- function(population, pool_size, tour_size) {
  selection <- nsga2R::tournamentSelection(population, pool_size, tour_size)
}

nsga2.custom <- function(dmatrix_expression, dmatrix_biological, num_clusters=6, generations=50, population_size=10, crossover_ratio=0.80, mutation_ratio=0.01, tour_size=2, neighborhood = 0.10, local_search=FALSE, ls_pos=FALSE) {
  
  # Sanity checks
  if( nrow( dmatrix_expression ) != ncol( dmatrix_expression ) || nrow( dmatrix_biological ) != ncol( dmatrix_biological ) ) {
    stop( "Both distance matrices (dmatrix_expression and dmatrix_biological) should be square matrices" )
  }
  if( identical( dim(dmatrix_expression), dim(dmatrix_biological) ) == FALSE  ) {
    stop("Both matrices should have the same dimensions")
  }
  if( num_clusters <= 0 ) {
    stop("num_clusters should be a positive number greater than 0")
  }
  # Order matrices just in case
  dmatrix_expression <- helper.order.matrix(dmatrix_expression)
  dmatrix_biological <- helper.order.matrix(dmatrix_biological)
  
  # Sum distance matrix
  #gene_dmatrix <- dmatrix_expression + dmatrix_biological
  # Normalized euclidean merged matrix
  #gene_dmatrix <- as.data.frame( as.matrix( pdist( dmatrix_expression, dmatrix_biological ) ) )
  #gene_dmatrix <- ( gene_dmatrix - min(gene_dmatrix) ) / ( max(gene_dmatrix) - min(gene_dmatrix) )
  
  #colnames(gene_dmatrix) <- colnames(dmatrix_expression)
  #rownames(gene_dmatrix) <- rownames(dmatrix_expression)
  
  gene_list <- colnames(dmatrix_expression)
  
  population_parents <- generate.initial.population(gene_list, population_size, num_clusters)
  population_parents <- operator.nsga2.sorting.and.crowding(population_size, population_parents, dmatrix_expression, dmatrix_biological)
  
  #operator.local.search(ls_pos == 1, population_size, population_parents, dmatrix_expression, dmatrix_biological, local_search, ls_params)
  g = 1
  
  while( g <= generations ) {
    selection <- operator.selection(population_parents, population_size, tour_size)
    population_parents <- operator.nsga2.sorting.and.crowding()
    
    population_children <- operator.crossover()
    population_children <- operator.mutation()
    
    population_mix = population_children + population_parents
    population_parents = population_mix[,population_size]
    
    #operator.local.search(ls_pos == 2, population_size, population_parents, dmatrix_expression, dmatrix_biological, local_search, ls_params)
    
    g = g + 1
  }
  #operator.local.search(ls_pos == 3, population_size, population_parents, dmatrix_expression, dmatrix_biological, local_search, ls_params)
  results = c()
  
  return( results )
}