packages = c("nsga2R")

operator.nsga2.sorting.and.crowding <- function() {
  ranking <- nsga2R::fastNonDominatedSorting(objectiveValues)
  crowdingDistance <- nsga2R::crowdingDist4frnt(population, ranking, objRange)
}

operator.local.search <- function(apply, local_search, ls_params, nsga_params, data) {
  if( apply == TRUE && local_search != FALSE ) {
    operator.nsga2.sorting.and.crowding()
  }
  return( data )
}

operator.crossover <- function(crossover_ratio) {
  
}

operator.mutation <- function(mutation_ratio) {
  
}

operator.selection <- function(population_size, pool_size, tour_size) {
  selection <- nsga2R::tournamentSelection(population_size, pool_size, tour_size)
}

nsga2.custom <- function(dmatrix_expression, dmatrix_biological, num_clusters, generations, population_size, crossover_ratio, mutation_ratio, tour_size, neighborhood, local_search=FALSE, ls_pos=FALSE) {
  
  # Sanity checks
  if( is.square.matrix( dmatrix_expression ) == FALSE || is.square.matrix( dmatrix_biological ) == FALSE ) {
    stop( "Both distance matrices (dmatrix_expression and dmatrix_biological) should be square matrices" )
  }
  if( identical( dim(dmatrix_expression), dim(dmatrix_biological) ) == FALSE  ) {
    stop("Both matrices should have the same dimensions")
  }
  if( num_clusters <= 0 ) {
    stop("num_clusters should be a positive number greater than 0")
  }
  population_parents <- generate.initial.population(population_size)
  population_parents <- operator.nsga2.sorting.and.crowding()
  
  operator.local.search(ls_pos == 1, local_search, ls_params, nsga_params, data)
  g = 1
  
  while( g <= generations ) {
    selection <- operator.selection(population_size, pool_size, tour_size)
    population_parents <- operator.nsga2.sorting.and.crowding()
    
    population_children <- operator.crossover()
    population_children <- operator.mutation()
    
    population_mix = population_children + population_parents
    population_parents = population_mix[,population_size]
    
    operator.local.search(ls_pos == 2, local_search, ls_params, nsga_params, data)
    
    g = g + 1
  }
  operator.local.search(ls_pos == 3, local_search, ls_params, nsga_params, data)
  results = c()
  
  return( results )
}