library(nsga2R)
library(hash)
library(ggplot2)
library(stringr)

source("local_search.r")
source("globals.r")

# Create counters to test algorithm's performance
fitness_counter <<- 0
mutation_counter <<- 0
not_mutated_counter <<- 0
reset_crossover_counter <<- 0
crossover_counter <<- 0

# Hash used to avoid calculating a function's objective value more than once.
fitness_hash <<- hash()

#' Orders a matrix's rownames and colnames assuming they can be transformed using
#' as.numeric() from lower to greater.
#' 
#' @param matrix Matrix to be ordered.
#' @return ordered matrix (rownames and colnames ordered).
#' 
helper.order.matrix <- function(matrix) {
  matrix[ order( as.numeric(rownames(matrix)) ), order( as.numeric(colnames(matrix)) ) ]
}

#' Generate a pareto ranking using the results of nsga2R::fastNonDominatedSorting.
#' 
#' @param popSize Integer value. Size of the solution population.
#' @param ranking Results of nsga2R::fastNonDominatedSorting. 
#' @return Vector of ranking of each solution.
#' @note The following code is part of the example in ?nsga2R::crowdingDist4frnt
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
  while(sum(duplicated(medoid_solution)) > 0 || length(medoid_solution[medoid_solution == 0]) > 0) {
    reset_crossover_counter <<- reset_crossover_counter + 1
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

operator.local.search <- function(should_apply, population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, local_search, ls_params = NULL, debug = FALSE) {
  if( should_apply == TRUE && !is.null(local_search) ) {
    population <- population[,1:num_clusters]
    
    if( local_search == local_search_algorithms$pls ) 
    {
      if( is.null(ls_params) ) {
        new_population <- local.search.pareto.local.search(population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, operator.nsga2.sorting.and.crowding, fitness.medoid.wg, debug = debug)
      } else {
        new_population <- local.search.pareto.local.search(ls_params$exploration_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, operator.nsga2.sorting.and.crowding, fitness.medoid.wg, ls_params$acceptance_criteria_fn, ls_params$rank_cutoff, debug = debug)
      }
    } 
    else if( local_search == local_search_algorithms$nmols ) 
    {
      new_population <- local.search.narrow.mols(ls_params$exploration_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, operator.nsga2.sorting.and.crowding, fitness.medoid.wg, debug = debug)
    }
    else if( local_search == local_search_algorithms$lmols )
    {
      new_population <- local.search.large.mols(ls_params$exploration_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, operator.nsga2.sorting.and.crowding, fitness.medoid.wg, debug = debug)
    }
    else if( local_search == local_search_algorithms$mosa )
    {
      new_population <- local.search.mosa(ls_params$exploration_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, operator.nsga2.sorting.and.crowding, fitness.medoid.wg, ls_params$acceptance_criteria_fn, ls_params$rank_cutoff, ls_params$alfa, debug = debug)
    }
    else if( local_search == local_search_algorithms$ensemble )
    {
      new_population <- local.search.ensemble(ls_params$exploration_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, operator.nsga2.sorting.and.crowding, fitness.medoid.wg, ls_params$rank_cutoff, debug = debug)
    }
    else {
      stop("Wrong local search algorithm specified")
    }
    
    new_population <- new_population[, 1:num_clusters]
    rownames(new_population) <- 1:nrow(new_population)
    
    population <- unique(rbind(population, new_population))
    #population <- helper.randomize.duplicates(population, gene_list, num_clusters)
    
    population <- operator.nsga2.sorting.and.crowding(population, dmatrix_expression, dmatrix_biological)
  }
  return( population )
}

operator.crossover.random <- function(gene_list, population_size, num_clusters, mating_pool, crossover_ratio) {
  population_size = round(population_size/2)
  population_children <- as.data.frame( matrix(0, population_size * 2, num_clusters), stringsAsFactors = FALSE )
  
  for( i in 1:population_size ) {
    pairs <- sample( 1:population_size, 2, replace=FALSE )
    crossover_counter <<- crossover_counter + 1
    
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
  } 
  
  return( population_children )
}

#' 
#' 
#' 
#' 
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
      
      mutation_counter <<- mutation_counter + 1
    } else {
      not_mutated_counter <<- not_mutated_counter + 1
    }
  }
  
  return( population_children )
}

operator.diversify.population <- function(gene_list, num_clusters, population) {
  
  # Generar función para saber que tan similares son dos individuos d(s1, s2) = k - #elementos iguales
  similitud = integer( nrow(population) * ( nrow(population) - 1 ) / 2 )
  i <- 1
  
  for(first_idx in 1:nrow(population)) {
    for(second_idx in (first_idx + 1):nrow(population)) {
      if( second_idx > nrow(population) || second_idx == first_idx) next
      
      first_solution <- population[first_idx, 1:num_clusters, drop=FALSE]
      second_solution <- population[second_idx, 1:num_clusters, drop=FALSE]
      
      similitud[[i]] <- length(intersect(first_solution, second_solution)) / num_clusters
      i <- i + 1
    }
  }
  factor_similitud <- mean(similitud)
  # Hacer algo si la similitud es muy alta
  
  return( factor_similitud )
}

#' 
#' 
#' 
#' 
operator.selection <- function(population, pool_size, tour_size) {
  return( nsga2R::tournamentSelection(population, pool_size, tour_size) )
}

#' 
#' 
#' 
#' 
operator.join.populations <- function(parent_pop, child_pop, num_clusters, gene_list) {
  parent_pop <- parent_pop[, 1:num_clusters]
  child_pop <- child_pop[, 1:num_clusters]
  
  mixed_pop <- rbind(parent_pop, child_pop)
  mixed_pop <- helper.randomize.duplicates(mixed_pop, gene_list, num_clusters)
  return( mixed_pop )
}

#' 
#' 
#' 
#' 
generate.results <- function(population_size, num_clusters, population, dmatrix_expression, dmatrix_biological) {
  
  # Select pareto front
  population <- population[ 1:population_size, ]
  population <- population[ population$solutions_rank == 1, ]

  # Mark to which cluster does each gene belong to in each solution
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

#' NSGA-II
#' 
#' @param dmatrix_expression Matrix or Data.frame. n^n matrix with gene expression distance.
#' @param dmatrix_biological Matrix or Data.frame. n^n matrix with gene biological distance.
#' This distance can be calculated using different biological sources.
#' @param evaluations Integer. Budget for evaluating new solutions.
#' @param population_size Intger. Number of solutions kept in the NSGA-II population.
#' @param crossover_ratio Probability of using the crossover operator.
#' @param crossover_prob Probability of taking solutions in rank 1 over solutions in the
#' lowest rank. 0.5 means 50% less chance of taking solutions in rank 1. 2.0 means 2x the
#' chance of taking solutions in rank 1 over rank max. Probability goes from this number
#' to 1.0 linearly over each rank. CURRENTLY NOT USED
#' @param mutation_ratio Probability of mutating solutions.
#' @param tour_size Integer. Size of the selection tournament.
#' @param neighborhodd Float. Max distance to consider two genes as neighbors. [0, 1.414]
#' @return List. population and clustering. Approximation to the Pareto Front and cluster
#' to which each gene belongs to in each solution in the Pareto Front respectively.
#' 
nsga2.custom <- function(dmatrix_expression, dmatrix_biological, num_clusters=5, evaluations=1000, population_size=20, crossover_ratio=0.60, crossover_prob=1.0, mutation_ratio=0.10, tour_size=2, neighborhood = 0.45, local_search=NULL, ls_pos=FALSE, ls_budget=60.0, ls_params=NULL, debug=FALSE, message_iteration=FALSE) {
  
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
  
  # Debug counters
  mutation_counter <<- 0
  not_mutated_counter <<- 0
  reset_crossover_counter <<- 0
  crossover_counter <<- 0
  
  # Contador de fitness y hash de fitness
  fitness_counter <<- 0
  clear(fitness_hash)
  if(exists("fitness_hash")) rm("fitness_hash", envir = globalenv())
  fitness_hash <<- hash()
  
  if( !is.null(local_search) ) {
    # Sum distance matrix
    #dmatrix_combined <- dmatrix_expression + dmatrix_biological
    # Euclidean distance matrix
    dmatrix_combined <- sqrt(dmatrix_expression**2 + dmatrix_biological**2)
    # Find genes that are close to one another
    neighborhood_matrix <- sapply(gene_list, function(gene) {
      neighborhood_genes <- dmatrix_combined[gene, , drop=FALSE]
      apply(neighborhood_genes, 1, function(x) colnames(neighborhood_genes)[which(x > 0.00 & x < neighborhood)] )
    })
  }

  # Divide budget size by 2 to take into account that each solution is evaluated twice
  original_evaluations <- evaluations
  evaluations <- round( ((evaluations - population_size * 2) * 1.04) / 2)
  
  # Setup global evaluations budget
  if (is.null(local_search)) {
    nsga_budget <- evaluations
  } else {
    nsga_budget <- evaluations * ( 1 - (ls_budget / 100.0))
    ls_budget <- evaluations * (ls_budget / 100.0)
  }
  # Setup nsga budget
  max_pool_size <- max(round(population_size / 2), 10)
  min_pool_size <- 2
  pool_size <- min(max(round(sqrt(nsga_budget)), min_pool_size), max_pool_size)
  
  generations <- round(nsga_budget / pool_size)
  
  # Setup local search budget
  if( !is.null(local_search) ) {
    if(ls_pos == 2) {
      ls_params["exploration_size"] <- round(ls_budget / generations)
    } else {
      ls_params["exploration_size"] <- round(ls_budget)
    }
  }
  
  if (debug) {
    message(str_interp("Running with a total budget = ${original_evaluations}"))
    message(str_interp("Running NSGA-II with pool_size=${pool_size} for ${generations} generations with a budget of ${nsga_budget*2}"))
    if ( !is.null(local_search) ) {
      message(str_interp("Running local search with a budget of ${ls_budget*2} using ${ls_params['exploration_size']} each run"))
    }
  }
  
  # Start NSGA-II loop
  population_parents <- generate.initial.population(gene_list, population_size, num_clusters)
  population_parents <- operator.nsga2.sorting.and.crowding(population_parents, dmatrix_expression, dmatrix_biological)
  
  population_parents <- operator.local.search(ls_pos == 1, population_size, population_parents, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, local_search, ls_params, debug)
  g = 1
  
  similitudes <- list()
  while( g <= generations) { # && fitness_counter < evaluations) {
    if (debug) {
      message(paste("running generation:", g))
      message(paste("start of generation:", fitness_counter))
    }
    
    mating_pool <- operator.selection(population_parents, pool_size, tour_size)
    
    population_children <- operator.crossover.random(gene_list, pool_size, num_clusters, mating_pool, crossover_ratio)
    population_children <- operator.mutation.random(gene_list, num_clusters, population_children, mutation_ratio)
    
    population_mix <- operator.join.populations(population_parents, population_children, num_clusters, gene_list)
    population_mix <- operator.nsga2.sorting.and.crowding(population_mix, dmatrix_expression, dmatrix_biological)
    population_parents <- population_mix[ 1:population_size, ]
    similitud <- operator.diversify.population(gene_list, num_clusters, population_children)
    similitudes[[g]] <- similitud
    if (debug) message(similitud)
    
    if (debug) message(paste("Before local search:", fitness_counter))
    population_parents <- operator.local.search(ls_pos == 2, population_size, population_parents, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, local_search, ls_params, debug)
    
    if (debug) message(paste("after local search:", fitness_counter))
    g = g + 1
    
    # message graph and/or debugging values, i.e. results of each iteration
    if (debug && message_iteration) {
      results <- generate.results(population_size, num_clusters, population_parents, dmatrix_expression, dmatrix_biological)
      metrics <- evaluator.multiobjective.clustering.no.bio( results, dmatrix_expression )
      
      if (debug) {
        message(population_parents)
        message(metrics)
      }
      
      if (message_iteration) {
        d <- results$population[, colnames(results$population) %in% c('objective_exp', 'objective_bio')]
        ggplot() +
          geom_step(data=d, mapping=aes(x=objective_exp, y=objective_bio)) +
          geom_step(data=d, mapping=aes(x=objective_exp, y=objective_bio), direction="vh", linetype=3) +
          geom_point(data=d, mapping=aes(x=objective_exp, y=objective_bio), color="red") + xlim(0, 1) + ylim(0, 1)
        
        ggsave(str_interp("iteracion-${g}.png"), device="png", path="plots")
      }
    }
  }
  
  # Prepare final results
  population_parents <- operator.local.search(ls_pos == 3, population_size, population_parents, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, local_search, ls_params, debug)
  
  results <- generate.results(population_size, num_clusters, population_parents, dmatrix_expression, dmatrix_biological)
  
  # Show similarity change from generation to generation
  if (debug) {
    message(paste("Se realizaron", fitness_counter, "calculos de la funcion de fitness"))
    
    data <- cbind(g = 1:length(similitudes), sim=unlist(similitudes))
    ggplot(as.data.frame(data), aes(x = g)) +
      geom_line(aes(y = sim), color = 'cyan') +
      geom_area(aes(y = sim), fill = 'cyan', alpha = .1) +
      #geom_line(aes(y = g2), color="steelblue") +
      #geom_area(aes(y = g2), fill = 'darkred', alpha = .1) +
      xlab('Generacion') +
      ylab('Jaccard\npromedio') +
      ggtitle('Similitud promedio entre todos las soluciones de la poblacion') + ylim(0, 1) +
      theme(text = element_text(color = "#222222")
            ,panel.background = element_rect(fill = '#444B5A')
            ,panel.grid.minor = element_line(color = '#4d5566')
            ,panel.grid.major = element_line(color = '#586174')
            ,plot.title = element_text(size = 16)
            ,axis.title = element_text(size = 12, color = '#333333')
      )
    
    
    ggsave(paste('jaccard_sim-', round(stats::runif(1, 1, 10000)), '.png', sep=''), device="png", path="plots")
  }
  
  return( results )
}
