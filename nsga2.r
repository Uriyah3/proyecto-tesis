library(nsga2R)
library(hash)
library(ggplot2)

source("local_search.r")
source("globals.r")

fitness_counter = 0
fitness_hash <- hash()

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
  
  # Esta función tenía un bug donde siempre agregaba un medoide a una población sin duplicados.
  if(nrow(duplicated_population) > 0) {
    for(i in 1:nrow(duplicated_population)) {
      duplicated_population[i, 1:num_clusters] = sample( gene_list, num_clusters, replace=F )
    }
    
    population <- rbind(unique_population, duplicated_population)
  }
  
  return( population )
}

generate.initial.population <- function(genes, population_size, num_clusters) {
  as.data.frame( t(sapply( 1:population_size, function(x) as.character(sample( genes, num_clusters, replace=F )) )), stringsAsFactors = FALSE )
}

objective.functions <- function(cluster_solutions, dmatrix, type=NULL) {
  objective_indices <- vector()
  for(i in 1:nrow(cluster_solutions)) {
    objective_indices[i] <- fitness.medoid.wg( cluster_solutions[i, ], dmatrix, type )
  }
  
  return( objective_indices )
}

medoid.fix.representation <- function(gene_list, medoid_solution, num_clusters) {
  medoid_solution <- as.list(medoid_solution)
  while(sum(duplicated(medoid_solution)) > 0) {
    medoid_solution[duplicated(medoid_solution)] <- sample(gene_list, sum(duplicated(medoid_solution)), replace=F)
  }
  return(as.data.frame(medoid_solution, stringsAsFactors = FALSE))
}

# fitness.medoid <- function(cluster_solution, gene_dmatrix) {
#   elements <- nrow( gene_dmatrix )
#   medoids <- length( cluster_solution )
#   distance_to_medoids <- gene_dmatrix[rownames(gene_dmatrix) %in% cluster_solution, ]
#   
#   medoid_index <- integer(medoids)  
#   # Using CLAV::iv.xb
#   for(medoid in 1:medoids) {
#     medoid_index[medoid] = grep( paste("^", cluster_solution[medoid], "$", sep=""), colnames(gene_dmatrix) )
#   }
# 
#   clustering <- apply(distance_to_medoids, 2, function(x) rownames(distance_to_medoids)[which.min(x)])
# 
#   XB <- tryCatch(
#     CLAV::iv.xb(clustering, unlist(cluster_solution), gene_dmatrix, variant="weighted"),
# 
#     error=function(cond) {
#       message(
#         paste(
#           "The following medoid solution had a problem: [",
#           paste(cluster_solution, collapse=", "), "]", sep=""
#         )
#       )
#       message("Here's the original error message:")
#       message(cond)
# 
#       return(Inf)
#     }
#   )
# 
#   return( XB )
#   
#   # Manual version below
#   # check if this distances are correct, shouldn't it sum only the cluster elements?
#   # i.e the minimum in each one
#   # distance_to_medoids <- distance_to_medoids ^ 2
#   # XB_numerator = sum( apply(distance_to_medoids, 2, min) )
#   # medoid_pairs = utils::combn( unlist(cluster_solution), 2 )
#   # separation = vector()
#   # for( col in 1:ncol(medoid_pairs) ) {
#   #   separation[col] <- ( gene_dmatrix[ t( medoid_pairs[, col] ) ] ) ^ 2
#   # }
#   # minimum_separation <- min(separation)
#   # XB_denominator = elements * minimum_separation
#   # XB <- XB_numerator / XB_denominator
# 
#   return( XB )
# }

fitness.medoid.wg <- function(cluster_solution, gene_dmatrix, type=NULL) {
  if(is.character(type)) {
    # Revisar si se ha calculado esta solución con esta matriz anteriormente
    key <- paste(paste(sort(cluster_solution), collapse=","), type)
    if(has.key(key, fitness_hash)) {
      return( fitness_hash[[key]] )
    }
  }
  
  fitness_counter <<- fitness_counter + 1
  # Calcular distancia mínima de cada gene a su cluster.
  elements <- nrow( gene_dmatrix )
  medoids <- length( cluster_solution )
  distance_to_medoids <- gene_dmatrix[rownames(gene_dmatrix) %in% cluster_solution, ]
  
  # Promediar distancia de un gen a su cluster con la minima distancia a otro cluster
  Rm <- sapply(distance_to_medoids, function(x) {if(length(x[!x %in% min(x)]) > 0) { min(x) / min(x[!x %in% min(x)])} else { 1.0 } } )
  
  # Sacar el Jk de cada cluster
  clustering <- apply(distance_to_medoids, 2, function(x) rownames(distance_to_medoids)[which.min(x)])
  Jk <- as.data.frame(cbind(Rm, clustering), stringsAsFactors = FALSE)
  Jk <- transform(Jk, Rm = as.numeric(Rm))
  elements_k <- aggregate(Rm ~ clustering, Jk, length)
  Jk <- aggregate(Rm ~ clustering, Jk, mean)
  Jk <- sapply(Jk$Rm, function(x) max(0, x))
  
  # Sacar el C promediado (necesita la cantidad de elementos por cluster)
  WG <- sum(Jk * elements_k$Rm) / elements
  
  if(is.character(type)) {
    # Guardar en el hash
    fitness_hash[[key]] <<- WG
  }
  
  return( WG )
}

operator.nsga2.sorting.and.crowding <- function(population, dmatrix_expression, dmatrix_biological) {
  
  population_size <- nrow(population)
  
  objective_exp <- objective.functions( population, dmatrix_expression, 'expression' )
  objective_bio <- objective.functions( population, dmatrix_biological, 'biological' )
  
  objectives <- cbind( objective_exp, objective_bio )
  obj_range <- ( apply(objectives, 2, max) - apply(objectives, 2, min) )
  
  ranking <- nsga2R::fastNonDominatedSorting( objectives )
  solutions_rank <- helper.pareto.ranking( population_size, ranking )
  population <- cbind( population, objectives, solutions_rank )
  
  crowding_distance <- nsga2R::crowdingDist4frnt( population, ranking, obj_range )
  
  population <- cbind( population, crowding_distance = c(apply(crowding_distance,1,sum)) )
  population <- population[order(population$solutions_rank, -(population$crowding_distance)), ]
  
  return( population )
}

operator.local.search <- function(should_apply, population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, local_search, ls_params = NULL) {
  if( should_apply == TRUE && !is.null(local_search) ) {
    population <- population[,1:num_clusters]
    
    if( local_search == local_search_algorithms$pls ) 
    {
      if( is.null(ls_params) ) {
        new_population <- local.search.pareto.local.search(population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, operator.nsga2.sorting.and.crowding, fitness.medoid.wg)
      } else {
        new_population <- local.search.pareto.local.search(population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, operator.nsga2.sorting.and.crowding, fitness.medoid.wg, ls_params$acceptance_criteria_fn, ls_params$rank_cutoff, ls_params$max_generations)
      }
    } 
    else if( local_search == local_search_algorithms$nmols ) 
    {
      new_population <- local.search.narrow.mols(population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, operator.nsga2.sorting.and.crowding, fitness.medoid.wg)
    }
    else if( local_search == local_search_algorithms$lmols )
    {
      new_population <- local.search.large.mols(population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, operator.nsga2.sorting.and.crowding, fitness.medoid.wg)
    }
    else if( local_search == local_search_algorithms$mosa )
    {
      new_population <- local.search.mosa(population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, operator.nsga2.sorting.and.crowding, fitness.medoid.wg)
    }
    else if( local_search == local_search_algorithms$ensemble )
    {
      new_population <- local.search.ensemble(population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, operator.nsga2.sorting.and.crowding, fitness.medoid.wg)
    }
    else {
      stop("Wrong local search algorithm specified")
    }
    
    new_population <- new_population[, 1:num_clusters]
    rownames(new_population) <- 1:nrow(new_population)
    
    population <- rbind(population, new_population)
    population <- helper.randomize.duplicates(population, gene_list, num_clusters)
    
    population <- operator.nsga2.sorting.and.crowding(population, dmatrix_expression, dmatrix_biological)
  }
  return( population )
}

operator.crossover.random <- function(gene_list, population_size, num_clusters, mating_pool, crossover_ratio) {
  population_children <- as.data.frame( matrix(0, population_size * 2, num_clusters), stringsAsFactors = FALSE )
  
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
  
  # Randomize any duplicated medoids in solutions
  for( i in 1:nrow(population_children) ) {
    population_children[i, 1:num_clusters] <- medoid.fix.representation(gene_list, population_children[i, 1:num_clusters], num_clusters)
  } # Generar función para saber que tan similares son dos individuos d(s1, s2) = k - #elementos iguales
  
  return( population_children )
}
# Contabilizar mutación
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
  mixed_pop <- helper.randomize.duplicates(mixed_pop, gene_list, num_clusters)
  return( mixed_pop )
}

generate.results <- function(population_size, num_clusters, population, dmatrix_expression, dmatrix_biological) {
  
  population <- population[ 1:population_size, ]
  population <- population[ population$solutions_rank == 1, ]
  
  frontier_clustering <- list(1:nrow(population))
  
  elements <- nrow( dmatrix_expression )
  for( i in 1:nrow(population) ) {
    cluster_solution <- population[i, 1:num_clusters]
    medoids <- length( cluster_solution )
    distance_to_medoids <- dmatrix_expression[rownames(dmatrix_expression) %in% cluster_solution, ]
    
    medoid_index <- integer(medoids)  
    for(medoid in 1:medoids) {
      medoid_index[medoid] = grep( paste("^", cluster_solution[medoid], "$", sep=""), colnames(dmatrix_expression) )
    }
    
    clustering <- apply(distance_to_medoids, 2, function(x) which.min(x))
    frontier_clustering[[i]] <- as.vector(clustering)
  }
  
  return( 
    list(
      population = population,
      clustering = frontier_clustering
    )
  )
}

# Possible implementation.
# graph / parent child population
# use two parameters: number of childs per node, number of levels
# Crossover would happen between same population and between other populations
# Parent comunicates with child, parents comunicate between the same level
#
# Currently using normal population.
nsga2.custom <- function(dmatrix_expression, dmatrix_biological, num_clusters=5, generations=50, population_size=20, crossover_ratio=0.60, mutation_ratio=0.10, tour_size=2, neighborhood = 0.45, local_search=NULL, ls_pos=FALSE, ls_params=NULL, debug=FALSE, print_iteration=FALSE) {
  
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
  if( !is.null(local_search) && !(local_search %in% local_search_algorithms) ) {
    stop( paste("The local search algorithm must be one of the following:", paste(local_search_algorithms, collapse=", ")) )
  }
  
  # Order matrices just in case
  dmatrix_expression <- helper.order.matrix(dmatrix_expression)
  dmatrix_biological <- helper.order.matrix(dmatrix_biological)
  
  gene_list <- colnames(dmatrix_expression)
  
  # Contador de fitness y hash de fitness
  fitness_counter <<- 0
  clear(fitness_hash)
  if(exists("fitness_hash")) rm("fitness_hash", envir = globalenv())
  fitness_hash <<- hash()
  
  if( !is.null(local_search) ) {
    # Sum distance matrix
    #dmatrix_combined <- dmatrix_expression + dmatrix_biological
    dmatrix_combined <- sqrt(dmatrix_expression**2 + dmatrix_biological**2)
    # Find genes that are close to one another
    neighborhood_matrix <- sapply(gene_list, function(gene) {
      neighborhood_genes <- dmatrix_combined[gene, , drop=FALSE]
      apply(neighborhood_genes, 1, function(x) colnames(neighborhood_genes)[which(x > 0.00 & x < neighborhood)] )
    })
  }
  
  population_parents <- generate.initial.population(gene_list, population_size, num_clusters)
  population_parents <- operator.nsga2.sorting.and.crowding(population_parents, dmatrix_expression, dmatrix_biological)
  
  population_parents <- operator.local.search(ls_pos == 1, population_size, population_parents, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, local_search, ls_params)
  g = 1
  
  while( g <= generations ) {
    if (debug) {
      print(paste("running generation:", g))
    }
    
    mating_pool <- operator.selection(population_parents, population_size, tour_size)
    
    population_children <- operator.crossover.random(gene_list, population_size, num_clusters, mating_pool, crossover_ratio)
    population_children <- operator.mutation.random(gene_list, num_clusters, population_children, mutation_ratio)
    
    population_mix <- operator.join.populations(population_parents, population_children, num_clusters, gene_list)
    population_mix <- operator.nsga2.sorting.and.crowding(population_mix, dmatrix_expression, dmatrix_biological)
    population_parents <- population_mix[ 1:population_size, ]
    
    population_parents <- operator.local.search(ls_pos == 2, population_size, population_parents, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, local_search, ls_params)
    
    g = g + 1
    
    if (debug || print_iteration) {
      results <- generate.results(population_size, num_clusters, population_parents, dmatrix_expression, dmatrix_biological)
      metrics <- evaluator.multiobjective.clustering( results, dmatrix_expression )
      
      if (debug) {
        print(population_parents)
        print(metrics)
      }
      
      if (print_iteration) {
        d <- results$population[, colnames(results$population) %in% c('objective_exp', 'objective_bio')]
        ggplot() +
          geom_step(data=d, mapping=aes(x=objective_exp, y=objective_bio)) +
          geom_step(data=d, mapping=aes(x=objective_exp, y=objective_bio), direction="vh", linetype=3) +
          geom_point(data=d, mapping=aes(x=objective_exp, y=objective_bio), color="red") + xlim(0, 1) + ylim(0, 1)
        
        ggsave(str_interp("iteracion-${g}.png"), device="png", path="plots")
      }
    }
  }
  population_parents <- operator.local.search(ls_pos == 3, population_size, population_parents, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, local_search, ls_params)
  
  results <- generate.results(population_size, num_clusters, population_parents, dmatrix_expression, dmatrix_biological)
  
  if (debug) {
    print(paste("Se realizaron", fitness_counter, "calculos de la funcion de fitness"))
  }
  
  return( results )
}