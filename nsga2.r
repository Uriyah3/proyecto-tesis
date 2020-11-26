library(nsga2R)
library(matrixcalc)
library(pdist)
library(CLAV)

source("local_search.r")

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

helper.randomize.duplicates <- function(population, gene_list, num_clusters) {
  duplicated_population <- population[ duplicated(population[, 1:num_clusters]), ]
  unique_population <- population[ !duplicated(population[, 1:num_clusters]), ]
  
  for(i in 1:nrow(duplicated_population)) {
    duplicated_population[i, 1:num_clusters] = sample( gene_list, num_clusters, replace=F )
  }
  population <- rbind(unique_population, duplicated_population)
  
  return( population )
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
  rownames(distance_to_medoids) <- cluster_solution
  colnames(distance_to_medoids) <- rownames(gene_dmatrix)

  for (medoid in cluster_solution) {
    for (gene in rownames(gene_dmatrix)) {
      distance_to_medoids[medoid, gene] <- (gene_dmatrix[medoid, gene]) ^ 2
    }
  }
  
  medoid_index <- integer(medoids)  
  # Using CLAV::iv.xb
  for(medoid in 1:medoids) {
    medoid_index[medoid] = grep( paste("^", cluster_solution[medoid], "$", sep=""), colnames(gene_dmatrix) )
  }
  
  # Might need to pre-calculate every clustering for every medoid (and for both distance matrices)
  # UPDATE: Wait... I can't do that, that would be almost brute forcing the problem, since the clustering
  # depends on the other medoids and I would need to run every combination. On that note
  # do I even need to square the distance to a medoid? I don't think I do, since
  # I just need the closest medoid to each gene. That is enough to find the clustering.
  clustering <- apply(distance_to_medoids, 2, function(x) rownames(distance_to_medoids)[which.min(x)])
  
  XB <- tryCatch(
    CLAV::iv.xb(clustering, unlist(cluster_solution), gene_dmatrix),
    
    error=function(cond) {
      message(
        paste(
          "The following clustering solution had a problem: [", 
          paste(clustering, collapse=", "), "]", sep=""
        )
      )
      message("Here's the original error message:")
      message(cond)
      
      return(Inf)
    }
  ) 
  
  return( XB )
  
  ## Manual version below
  # # check if this distances are correct, shouldn't it sum only the cluster elements?
  # # i.e the minimum in each one
  # 
  # XB_numerator = sum( apply(distance_to_medoids, 2, min) )
  # 
  # medoid_pairs = utils::combn( unlist(cluster_solution), 2 )
  # 
  # separation = vector()
  # for( col in 1:ncol(medoid_pairs) ) {
  #   separation[col] <- ( gene_dmatrix[ t( medoid_pairs[, col] ) ] ) ^ 2
  # }
  # minimum_separation <- min(separation)
  # 
  # XB_denominator = elements * minimum_separation
  # 
  # return( XB_numerator / XB_denominator )
}

operator.nsga2.sorting.and.crowding <- function(population, dmatrix_expression, dmatrix_biological) {
  
  population_size <- nrow(population)
  
  objective_exp <- objective.functions( population, dmatrix_expression )
  objective_bio <- objective.functions( population, dmatrix_biological )
  
  objectives <- cbind( objective_exp, objective_bio )
  obj_range <- ( apply(objectives, 2, max) - apply(objectives, 2, min) )
  
  ranking <- nsga2R::fastNonDominatedSorting( objectives )
  solutions_rank <- helper.pareto.ranking( population_size, ranking )
  population <- cbind( population, objectives, solutions_rank )
  
  crowding_distance <- nsga2R::crowdingDist4frnt( population, ranking, obj_range )
  
  population <- cbind( population, crowding_distance = c(apply(crowding_distance,1,sum)) )
  population<-population[order(population$solutions_rank, -(population$crowding_distance)), ]
  
  return( population )
}

operator.local.search <- function(should_apply, population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, local_search) {
  if( should_apply == TRUE && !is.null(local_search) ) {
    population <- population[,1:num_clusters]
    
    if( local_search == "pls" ) 
    {
      new_population <- local.search.pareto.local.search(population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, operator.nsga2.sorting.and.crowding, fitness.medoid)
    } 
    else if( local_search == "nmols" ) 
    {
      new_population <- local.search.narrow.mols(population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, operator.nsga2.sorting.and.crowding, fitness.medoid)
    }
    else if( local_search == "lnols" )
    {
      new_population <- local.search.large.mols(population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, operator.nsga2.sorting.and.crowding, fitness.medoid)
    }
    else if( local_search == "mosa" )
    {
      new_population <- local.search.mosa(population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, operator.nsga2.sorting.and.crowding, fitness.medoid)
    }
    else if( local_search == "ensemble" )
    {
      new_population <- local.search.ensemble(population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, operator.nsga2.sorting.and.crowding, fitness.medoid)
    }
    else if( local_search == "pr" )
    {
      new_population <- local.search.path.relinking(population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, operator.nsga2.sorting.and.crowding, fitness.medoid)
    }
    
    new_population <- new_population[, 1:num_clusters]
    rownames(new_population) <- 1:nrow(new_population)
    
    population <- rbind(population, new_population)
    population <- helper.randomize.duplicates(population, gene_list, num_clusters)
    
    population <- operator.nsga2.sorting.and.crowding(population, dmatrix_expression, dmatrix_biological)
  }
  return( population )
}

operator.crossover.random <- function(population_size, num_clusters, mating_pool, crossover_ratio) {
  population_children <- as.data.frame( matrix(0, population_size * 2, num_clusters) )
  
  for( i in 1:population_size ) {
    pairs <- sample( 1:population_size, 2, replace=FALSE )
    
    for( medoid in 1:num_clusters ) {
      should_cross <- ( stats::runif(1, 0, 1) <= crossover_ratio )
      
      if( should_cross == TRUE ) {
        population_children[ i, medoid ] <- mating_pool[ pairs[2], medoid ]
        population_children[ i + population_size, medoid ] <- mating_pool[ pairs[1], medoid ]
      } else {
        population_children[ i, medoid ] <- mating_pool[ pairs[1], medoid ]
        population_children[ i + population_size, medoid ] <- mating_pool[ pairs[2], medoid ]
      }
    }
  }
  
  # Check if any solution has the same medoid multiple times
  # NOTE: It is extremely important to make sure that no solution has the same
  # medoid/gene multiple times, since it will screw with the xie_beni method
  # being used and it WILL crash the program.
  for( i in 1:nrow(population_children) ) {
    solution <- population_children[i, 1:num_clusters]
    if( sum( duplicated(t(solution)) ) > 0 ) {
      population_children <- population_children[-i, ]
    }
  }
  
  return( population_children )
}

operator.mutation.random <- function(gene_list, num_clusters, population_children, mutation_ratio) {
  for( i in nrow(population_children) ) {
    should_mutate <- stats::runif(1, 0, 1) <= mutation_ratio
    if( should_mutate ) {
      mutate_index <- sample(1:num_clusters, 1)
      
      repeat {
        mutated_medoid <- sample(gene_list, 1)
        
        if( !(mutated_medoid %in% population_children[i, ]) ) {
          break
        }
      }
      population_children[i, mutate_index] <- mutated_medoid
    }
  }
  
  return( population_children )
}

operator.selection <- function(population, pool_size, tour_size) {
  return( nsga2R::tournamentSelection(population, pool_size, tour_size) )
}

operator.join.populations <- function(parent_pop, child_pop, num_clusters, gene_list) {
  parent_pop <- parent_pop[, 1:num_clusters]
  child_pop <- child_pop[, 1:num_clusters]
  
  mixed_pop <- rbind(parent_pop, child_pop)
  helper.randomize.duplicates(mixed_pop, gene_list, num_clusters)
  return( mixed_pop )
}

# graph / parent child population
# use two parameters: number of childs per node, number of levels
# Crossover would happen between same population and between other populations
# Parent comunicates with child, parents comunicate between the same level
#
# Currently using normal population.
nsga2.custom <- function(dmatrix_expression, dmatrix_biological, num_clusters=6, generations=50, population_size=20, crossover_ratio=0.60, mutation_ratio=0.10, tour_size=2, neighborhood = 0.50, local_search=NULL, ls_pos=FALSE) {
  
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
  
  gene_list <- colnames(dmatrix_expression)
  
  if( !is.null(local_search) ) {
    # Sum distance matrix
    dmatrix_combined <- dmatrix_expression + dmatrix_biological
    # Find genes that are close to one another
    neighborhood_matrix <- sapply(gene_list, function(gene) {
      neighborhood_genes <- dmatrix_combined[gene, , drop=FALSE]
      apply(neighborhood_genes, 1, function(x) colnames(neighborhood_genes)[which(x > 0.00 & x < neighborhood)] )
    })
  }
  
  population_parents <- generate.initial.population(gene_list, population_size, num_clusters)
  population_parents <- operator.nsga2.sorting.and.crowding(population_parents, dmatrix_expression, dmatrix_biological)
  
  population_parents <- operator.local.search(ls_pos == 1, population_size, population_parents, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, local_search)
  g = 1
  
  while( g <= generations ) {
    mating_pool <- operator.selection(population_parents, population_size, tour_size)
    
    population_children <- operator.crossover.random(population_size, num_clusters, mating_pool, crossover_ratio)
    population_children <- operator.mutation.random(gene_list, num_clusters, population_children, mutation_ratio)
    
    population_mix <- operator.join.populations(population_parents, population_children, num_clusters, gene_list)
    population_mix <- operator.nsga2.sorting.and.crowding(population_mix, dmatrix_expression, dmatrix_biological)
    population_parents = population_mix[ 1:population_size, ]
    
    population_parents <- operator.local.search(ls_pos == 2, population_size, population_parents, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, local_search)
    
    g = g + 1
  }
  population_parents <- operator.local.search(ls_pos == 3, population_size, population_parents, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, local_search)
  results = population_parents
  
  return( results )
}