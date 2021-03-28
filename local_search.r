
#' Generates neighboring solutions (to a medoid representation) using a
#' pre-generated neighborhood matrix.
#' The generation method is as follows. From the medoid, 1 to num_clusters genes 
#' are randomly selected, each of these genes is replaced with a random neighboring
#' gene from the \code{neighborhood_matrix}. A gene is only replaced if the
#' replacement isn't currently in the solution.
#' 
#' @param exploration_size The amount of neighboring solutions to generate
#' @param num_clusters The number of clusters in a solution
#' @param solution A matrix with one solution.
#' @param neighborhood_matrix An adjacency matrix with all neighboring genes.
#' @param row_name_id An unique identifier that will be used to identify each row.
#' Its not changed inside this method, so when it is used it should be manually
#' updated after running this method.
#' @return A matrix with neighboring solutions.
#' 
helper.generate.neighborhood <- function(exploration_size, num_clusters, solution, neighborhood_matrix, row_name_id, debug = FALSE) {
  if (debug) message("Starting neighborhood exploration...")
  medoid_neighborhood <- as.data.frame( matrix(0, exploration_size, (num_clusters)), stringsAsFactors = FALSE )
  
  if (debug) message("Filtering genes with neigbors...")
  genes_with_neighbors <- solution[, 1:num_clusters]
  genes_with_neighbors <- genes_with_neighbors[(sapply(genes_with_neighbors, function(gene) {
    return( length(helper.get.neighborhood.gene(neighborhood_matrix, gene)) > 0 )
  })) ]
  
  if (debug) message( "The following genes have neighbors" )
  if (debug) message( genes_with_neighbors )
  for(medoid_generated in 1:nrow(medoid_neighborhood)) {
    medoid <- solution
    
    # Find genes that have at least 1 neighbor
    
    # Random neighborhood
    # Change a random amount of genes in the solution to a random neighboring gene.
    # genes_to_change <- sample(colnames(genes_with_neighbors), sample(1:length(genes_with_neighbors), 1))
    # for(gene_index in 1:length(genes_to_change)) {
    #   gene <- medoid[, genes_to_change[gene_index] ]
    #   neighbor_gene <- sample( neighborhood_matrix[[gene]], 1 )
    # 
    #   if( !(neighbor_gene %in% medoid[, 1:num_clusters]) ) {
    #     medoid[, gene_index] <- neighbor_gene
    #   }
    # }
    
    # 1 change to neighborhood
    # Changes a single gene to one of its neighbors
    if (length(genes_with_neighbors) > 0) {
      gene <- sample(genes_with_neighbors, 1)[1, ]
      neighbor_gene <- sample( helper.get.neighborhood.gene(neighborhood_matrix, gene), 1 )
      
      gene_index <- sample(1:num_clusters, 1)
      if( !(neighbor_gene %in% medoid[, 1:num_clusters]) ) {
        medoid[, gene_index] <- neighbor_gene
      }
    } else {
      # No neighbors, no change in the medoid
      if (debug) message(paste("Medoid ", medoid, "has no neighbors for current configuration"))
    }
    
    medoid_neighborhood[medoid_generated, ] <- medoid[, 1:num_clusters]
  }
  
  #message( paste("Neighbors found:", nrow(medoid_neighborhood), "without duplication")  )
  medoid_neighborhood <- medoid_neighborhood[ !duplicated(medoid_neighborhood[, 1:num_clusters]), ]
  #message( paste("Neighbors found:", nrow(medoid_neighborhood), "removing duplication")  )
  medoid_neighborhood <- cbind(medoid_neighborhood, add=rep( FALSE,nrow(medoid_neighborhood) ) )
  rownames(medoid_neighborhood) <- sapply(rownames(medoid_neighborhood), function(unused) {
    row_name_id <<- row_name_id + 1
    paste("V", row_name_id, sep="")
  })
  
  if (debug) message("neighborhood exploration DONE")
  return(medoid_neighborhood)
}

helper.get.neighborhood.gene <- function(neighborhood_matrix, gene) {
  if (is.list(neighborhood_matrix)) {
    if (gene %in% names(neighborhood_matrix)) {
      return(neighborhood_matrix[[gene]])
    } else {
      message("----BUG----BUG----BUG----BUG----BUG----BUG----")
      message(paste("Couldn't find the following gene in the neighborhood matrix:", gene))
      message("----BUG----BUG----BUG----BUG----BUG----BUG----")
      return( list() )
    } 
  } else {
    # Todos los genes son vecinos con todos
    return(neighborhood_matrix[neighborhood_matrix != gene])
  }
}

#' Helper for the Clustering ensemble local search method.
#' Generates all combinations between two medoids solutions without changing the
#' order of the genes in each medoid. In other words, creates 2^num_clusters -2
#' new solutions. The -2 is to skip creating the solutions equal to the medoids
#' being combined (\code{first_solution} and \code{second_solution}).
#' Since this calculation is exponential, it shouldn't be used with a num_clusters
#' greater than 10.
#' 
#' @param first_solution Vector or data.frame of genes. A medoid represented solution.
#' @param second_solution Vector or data.frame of genes. A medoid represented solution.
#' @param row_name_id An unique identifier that will be used to identify each row.
#' @return All medoids that can be generated by swapping the elements between
#' \code{first_solution} and \code{second_solution}.
#' @note It's not recommended to use ce local search with num_clusters greater than 8.
#' 
helper.generate.ensemble <- function(first_solution, second_solution, row_name_id) {
  num_clusters <- ncol(first_solution)
  medoid_ensemble <- as.data.frame( matrix(0, (2**num_clusters) - 2, num_clusters), stringsAsFactors = FALSE )
  
  for(row in 1:( nrow(medoid_ensemble) ) ) {
    for(col in 1:num_clusters) {
      if( bitwAnd(row, ( 2**(col-1) )) != 0 ) medoid_ensemble[row, col] <- first_solution[col]
      else medoid_ensemble[row, col] <- second_solution[col]        
    }
  }
  
  medoid_ensemble <- medoid_ensemble[ !duplicated(medoid_ensemble[, 1:num_clusters]), ]
  for( i in 1:nrow(medoid_ensemble) ) {
    solution <- medoid_ensemble[i, 1:num_clusters]
    if( sum( duplicated(t(solution)) ) > 0 ) {
      medoid_ensemble <- medoid_ensemble[-i, ]
    }
  }
  rownames(medoid_ensemble) <- sapply(rownames(medoid_ensemble), function(unused) {
    row_name_id <<- row_name_id + 1
    paste("V", row_name_id, sep="")
  })
  
  return(medoid_ensemble)
}

#' Helper for the Clustering ensemble local search method.
#' Generates some combinations between two medoids solutions. It may not keep the same
#' gene ordering. To be used when num_clusters is greater than 10.
#' 
#' @param first_solution Vector or data.frame of genes. A medoid represented solution.
#' @param second_solution Vector or data.frame of genes. A medoid represented solution.
#' @param exploration_size Integer. Number of solutions to generate.
#' @param row_name_id An unique identifier that will be used to identify each row.
#' @return Some solutions that can be generated by swapping the elements between
#' \code{first_solution} and \code{second_solution}.
#' 
helper.generate.partial.ensemble <- function(first_solution, second_solution, exploration_size, row_name_id) {
  num_clusters <- ncol(first_solution)
  medoid_ensemble <- as.data.frame( matrix(0, exploration_size, num_clusters), stringsAsFactors = FALSE )
  genes <- cbind(first_solution, second_solution)
  genes <- genes[!duplicated(as.list(genes))]
  
  for (medoid_solution in 1:exploration_size) {
    medoid_ensemble[medoid_solution, ] <- sort(sample(genes, num_clusters))
  }
  
  rownames(medoid_ensemble) <- sapply(rownames(medoid_ensemble), function(unused) {
    row_name_id <<- row_name_id + 1
    paste("V", row_name_id, sep="")
  })
  
  return(medoid_ensemble)
}

#' Calculate the probability based on the temperature and the difference of energy
#' between the new and old solution. Based on metaheuristics book by El-Ghazali Talbi
#' Uses boltzmann distribution as the probability.
#' 
#' @param objective_exp Float value. Expression based objective of solution being evaluated.
#' @param objective_bio Float value. Biologic based objective of solution being evaluated.
#' @param solution Solution compared to. Data frame or list that has to contain the
#' following columns: $objective_exp and $objective_bio
#' @return Value between 0 and 1. Probability of accepting this solution.
#'
helper.mosa.probability <- function(objective_exp, objective_bio, solution, temperature)
{
  # Energy should be fixed to consider other solutions. A new nondominated solution will only
  # be good if it also nondominates other solutions on the archive.
  #Energy <- sqrt((solution$objective_exp - objective_exp)**2 + (solution$objective_bio - objective_bio)**2)
  Energy <- (0.2) - (solution$objective_exp - objective_exp)**2 + (solution$objective_bio - objective_bio)
  probability <- exp(-(Energy/temperature)) # Probablidad boltzman
}

#' Modifies the archive adding the new found solutions and reorders the rank of the
#' solutions. It cuts off solutions with rank worse than \code{rank_cutoff}.
#' It maintains information referencing whether the solutions were explored or not.
#' 
#' @param archive Data frame of medoid solutions
#' @param add_to_archive Data frame of new medoid solutions
#' @param gene_list Vector of all available genes.
#' @param num_cluster Integer of the number of clusters.
#' @param ordering_fn Function used to order the archive
#' @param rank_cutoff Integer. Solutions with a rank > than this get deleted.
#' @return archive of solutions with the new solutions added as unexplored after
#' cutting off all solutions with a bad rank.
#' 
helper.cut.and.join.archive <- function(archive, add_to_archive, gene_list, num_clusters, ordering_fn, rank_cutoff) {
  explored <- archive[ , "explored", drop=FALSE]
  new_archive <- helper.randomize.duplicates( rbind( archive[, 1:num_clusters], add_to_archive), gene_list, num_clusters )
  archive <- ordering_fn(new_archive, dmatrix_expression, dmatrix_biological)
  archive <- archive[archive$solutions_rank <= rank_cutoff, ]
  
  new_solutions <- rownames( archive[!( rownames(archive) %in% rownames(explored) ), ] )
  new_solutions <- as.data.frame(rep( FALSE,length(new_solutions) ), stringsAsFactors = FALSE)
  rownames(new_solutions) <- rownames( archive[!( rownames(archive) %in% rownames(explored) ), ] )
  colnames(new_solutions) <- "explored"
  
  explored <- explored[ rownames(explored) %in% rownames(archive) , , drop=FALSE]
  explored <- rbind( new_solutions, explored )
  
  explored <- explored[ order(row.names(explored)), , drop=FALSE ]
  archive <- archive[ order(row.names(archive)), , drop=FALSE ]
  archive <- cbind( archive, explored )
  
  return( archive )
}

#' Non-dominance acceptance function. The new solution is better than the other solution
#' in at least one objective.
#' 
#' @param objective_exp Float value. Expression based objective of solution being evaluated.
#' @param objective_bio Float value. Biologic based objective of solution being evaluated.
#' @param solution Solution compared to. Data frame or list that has to contain the
#' following columns: $objective_exp and $objective_bio
#' @return boolean value indicating if the solution passes the non-dominance criteria.
#' 
helper.non.dominated <- function(objective_exp, objective_bio, solution) {
  return( objective_bio < solution$objective_bio || objective_exp < solution$objective_exp )
}

#' Dominance acceptance function. The new solution is better than the other solution in both 
#' objectives or, at worst, is equal in just one objective.
#' 
#' @param objective_exp Float value. Expression based objective of solution being evaluated.
#' @param objective_bio Float value. Biologic based objective of solution being evaluated.
#' @param solution Solution compared to. Data frame or list that has to contain the
#' following columns: $objective_exp and $objective_bio
#' @return boolean value indicating if the solution passes the dominance criteria.
#' 
helper.dominates <- function(objective_exp, objective_bio, solution) {
  return( ( objective_bio <= solution$objective_bio && objective_exp < solution$objective_exp ) ||
          ( objective_bio < solution$objective_bio && objective_exp <= solution$objective_exp ) )
}

helper.add.evaluation <- function(fitness_helper, evaluations) {
  if (fitness_counter > fitness_helper) {
    return (evaluations + 1)
  }
  return (evaluations)
}

#' Pareto Local Search
#' 
#' @param exploration_size Integer value. How many solutions should be explored each
#' iteration. Usually, the population_size of NSGA-II is used.
#' @param population Matrix / data.frame. Each row represents a solution.
#' @param num_clusters Integer value. How many medoids each solution has.
#' @param gene_list Vector. All the genes that can be used in the population. These
#' names should match the rows and cols of dmatrix_expression and dmatrix_biological. 
#' @param dmatrix_expression Gene expression distance matrix between each pair of genes.
#' @param dmatrix_biological Biological distance matrix between each pair of genes.
#' @param neighborhood_matrix Distance matrix. How close each pair of genes is based on
#' its expression and biological distance.
#' @param ordering_fn Function used to sort the population. 
#' Usually uses operator.nsga2.sorting.and.crowding .
#' @param fitness_fn Function used to calculate objectives. Usually uses fitness.medoid.wg
#' @param acceptance_criteria_fn Function used to accept a new solution into the
#' population archive. 3 parameters, new_objective_exp, new_objective_bio and old_solution
#' @param rank_cutoff Integer value. Solutions below this rank are removed from the archive.
#' @param max_generations Integer value. Maximum number of iterations.
#' 
#' @return matrix of new pareto front found by PLS.
#' 
local.search.pareto.local.search <- function(exploration_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, ordering_fn, fitness_fn, acceptance_criteria_fn = helper.non.dominated, rank_cutoff = 2, debug=FALSE) {
  
  # Create the archive of solutions
  archive <- ordering_fn(population[, 1:num_clusters], dmatrix_expression, dmatrix_biological)
  archive <- cbind( archive, explored=rep( FALSE,nrow(archive) ) )
  archive <- archive[archive$solutions_rank <= rank_cutoff, ]
  
  new_pool_size <- round(sqrt(exploration_size))
  evaluations <- 0
  row_name_id <- 0
  
  # Loop until alocated budget is used or until there are no more unexplored solutions
  while( sum(archive[, ncol(archive)]) < nrow(archive) && evaluations < exploration_size ) {
    if (debug) {
      message(paste("Evaluations: ", evaluations, "/", exploration_size, sep=""))
    }
    # Choose an unexplored solution
    unexplored <- archive[archive$explored == FALSE, , drop=FALSE]
    row_index <- rownames( unexplored[ sample(nrow(unexplored), 1), ] )[1]
    solution <- unexplored[row_index, , drop=FALSE]
    
    # Generate neighborhood
    medoid_neighborhood <- helper.generate.neighborhood(new_pool_size, num_clusters, solution, neighborhood_matrix, row_name_id, debug)
    row_name_id <- row_name_id + exploration_size + 1
    
    # Check if any neighboring solutions pass the acceptance criteria
    for(i in 1:nrow(medoid_neighborhood)) {
      if (evaluations < exploration_size) {
        fitness_helper <- fitness_counter
        
        objective_exp <- fitness_fn( medoid_neighborhood[i, 1:num_clusters], dmatrix_expression, 'expression' )
        objective_bio <- fitness_fn( medoid_neighborhood[i, 1:num_clusters], dmatrix_biological, 'biological' )
        
        evaluations <- helper.add.evaluation(fitness_helper, evaluations)
        
        if( acceptance_criteria_fn(objective_exp, objective_bio, solution) ) {
          medoid_neighborhood[i, "add"] <- TRUE
        }
      }
      
    }
    
    # Add solutions found in the neighborhood if there are any that passed the
    # acceptance criteria used
    if( sum(medoid_neighborhood[, ncol(medoid_neighborhood)]) > 1 ) {
      archive <- helper.cut.and.join.archive(
        archive, 
        medoid_neighborhood[ medoid_neighborhood$add == TRUE, 1:(num_clusters) ],
        gene_list, num_clusters, ordering_fn, rank_cutoff = rank_cutoff
      )
    }
    
    # Mark solution as explored if it still is in the archive
    if( row_index %in% rownames(archive) ) {
      archive[row_index, ncol(archive)] = TRUE
    }
  }
  
  # Clean the results and return them
  archive <- helper.randomize.duplicates( archive, gene_list, num_clusters )
  archive <- ordering_fn(archive[, 1:num_clusters], dmatrix_expression, dmatrix_biological)
  
  # Repeat PLS until all alocated budget is used. Avoid infinite recursion if PLS cant find new solutions.
  if (evaluations < exploration_size && evaluations > 0) {
    archive <- local.search.pareto.local.search(exploration_size - evaluations, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, ordering_fn, fitness_fn, acceptance_criteria_fn, rank_cutoff)
  }
  
  if (debug) message( paste("Pareto local search found", evaluations, "new solutions") )
  return( archive )
}

#' Large Multiobjective local search
#' 
#' @param exploration_size Integer value. How many solutions should be explored each
#' iteration. Usually, the population_size of NSGA-II is used.
#' @param population Matrix / data.frame. Each row represents a solution.
#' @param num_clusters Integer value. How many medoids each solution has.
#' @param gene_list Vector. All the genes that can be used in the population. These
#' names should match the rows and cols of dmatrix_expression and dmatrix_biological. 
#' @param dmatrix_expression Gene expression distance matrix between each pair of genes.
#' @param dmatrix_biological Biological distance matrix between each pair of genes.
#' @param neighborhood_matrix Distance matrix. How close each pair of genes is based on
#' its expression and biological distance.
#' @param ordering_fn Function used to sort the population. 
#' Usually uses operator.nsga2.sorting.and.crowding .
#' @param fitness_fn Function used to calculate objectives. Usually uses fitness.medoid.wg
#' 
#' @return Matrix of neighboring solutions that pass the non-dominance criteria.
#' 
local.search.large.mols <- function(exploration_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, ordering_fn, fitness_fn, debug=FALSE) {
  
  population <- ordering_fn(population[, 1:num_clusters], dmatrix_expression, dmatrix_biological)
  archive <- population[1, 1:num_clusters, drop=FALSE]
  
  new_pool_size <- round(exploration_size / nrow(population))
  evaluations <- 0
  row_name_id <- 0
  
  # Explore the neighbors of every solution in the population
  for(row_index in 1:nrow(population)) {
    if (debug) {
      message(paste("Evaluations: ", evaluations, "/", exploration_size, sep=""))
    }
    if (evaluations >= exploration_size) {
      break
    }
    # Actualizar pool_size
    new_pool_size <- round((exploration_size - evaluations) / (nrow(population) - row_index + 1))
    
    solution <- population[row_index, , drop=FALSE]
    medoid_neighborhood <- helper.generate.neighborhood(new_pool_size, num_clusters, solution, neighborhood_matrix, row_name_id, debug)
    row_name_id <- row_name_id + new_pool_size + 1
    
    for(i in 1:nrow(medoid_neighborhood)) {
      if (evaluations < exploration_size) {
        fitness_helper <- fitness_counter
      
        objective_exp <- fitness_fn( medoid_neighborhood[i, 1:num_clusters], dmatrix_expression, 'expression' )
        objective_bio <- fitness_fn( medoid_neighborhood[i, 1:num_clusters], dmatrix_biological, 'biological' )
      
        evaluations <- helper.add.evaluation(fitness_helper, evaluations)
      
        if( helper.non.dominated(objective_exp, objective_bio, solution) ) {
          archive <- rbind(archive, medoid_neighborhood[i, 1:num_clusters])
        }
      }
    }
  }
  
  return( archive )
}

#' Narrow Multiobjective local search
#' 
#' @param exploration_size Integer value. How many solutions should be explored each
#' iteration. Usually, the population_size of NSGA-II is used.
#' @param population Matrix / data.frame. Each row represents a solution.
#' @param num_clusters Integer value. How many medoids each solution has.
#' @param gene_list Vector. All the genes that can be used in the population. These
#' names should match the rows and cols of dmatrix_expression and dmatrix_biological. 
#' @param dmatrix_expression Gene expression distance matrix between each pair of genes.
#' @param dmatrix_biological Biological distance matrix between each pair of genes.
#' @param neighborhood_matrix Distance matrix. How close each pair of genes is based on
#' its expression and biological distance.
#' @param ordering_fn Function used to sort the population. 
#' Usually uses operator.nsga2.sorting.and.crowding .
#' @param fitness_fn Function used to calculate objectives. Usually uses fitness.medoid.wg
#' 
#' @return Matrix of neighboring solutions that pass the dominance criteria.
#' 
local.search.narrow.mols <- function(exploration_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, ordering_fn, fitness_fn, debug=FALSE) {
  
  population <- ordering_fn(population[, 1:num_clusters], dmatrix_expression, dmatrix_biological)
  archive <- population[1, 1:num_clusters, drop=FALSE]
  
  new_pool_size <- round(exploration_size / nrow(population))
  evaluations <- 0
  row_name_id <- 0
  
  # Explore the neighbors of every solution in the population
  for(row_index in 1:nrow(population)) {
    if (debug) {
      message(paste("Evaluations: ", evaluations, "/", exploration_size, sep=""))
    }
    if (evaluations >= exploration_size) {
      break
    }
    # Actualizar pool_size
    new_pool_size <- round((exploration_size - evaluations) / (nrow(population) - row_index + 1))
    
    solution <- population[row_index, , drop=FALSE]
    medoid_neighborhood <- helper.generate.neighborhood(new_pool_size, num_clusters, solution, neighborhood_matrix, row_name_id, debug)
    row_name_id <- row_name_id + new_pool_size + 1
    
    for(i in 1:nrow(medoid_neighborhood)) {
      if (evaluations < exploration_size) {
        fitness_helper <- fitness_counter
        
        objective_exp <- fitness_fn( medoid_neighborhood[i, 1:num_clusters], dmatrix_expression, 'expression' )
        objective_bio <- fitness_fn( medoid_neighborhood[i, 1:num_clusters], dmatrix_biological, 'biological' )
        
        evaluations <- helper.add.evaluation(fitness_helper, evaluations)
        
        if( helper.dominates(objective_exp, objective_bio, solution) ) {
          archive <- rbind(archive, medoid_neighborhood[i, 1:num_clusters])
        }
      }
    }
  }
  
  return( archive )
}

#' Multiobjective simulated annealing
#' 
#' @param exploration_size Integer value. How many solutions should be explored each
#' iteration. Usually, the population_size of NSGA-II is used.
#' @param population Matrix / data.frame. Each row represents a solution.
#' @param num_clusters Integer value. How many medoids each solution has.
#' @param gene_list Vector. All the genes that can be used in the population. These
#' names should match the rows and cols of dmatrix_expression and dmatrix_biological. 
#' @param dmatrix_expression Gene expression distance matrix between each pair of genes.
#' @param dmatrix_biological Biological distance matrix between each pair of genes.
#' @param neighborhood_matrix Distance matrix. How close each pair of genes is based on
#' its expression and biological distance.
#' @param ordering_fn Function used to sort the population. 
#' Usually uses operator.nsga2.sorting.and.crowding .
#' @param fitness_fn Function used to calculate objectives. Usually uses fitness.medoid.wg
#' @param rank_cutoff Integer value. Solutions below this rank are removed from the archive.
#' @param alfa Float value. Factor used for geomtric cooling.
#' @param max_steps Integer value. How many steps does the cooling goes on for. Iterations.
#' 
#' @return Matrix of neighboring solutions that pass the non-dominance criteria.
#' 
local.search.mosa <- function(exploration_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, ordering_fn, fitness_fn, acceptance_criteria_fn = helper.non.dominated, rank_cutoff = 5, alfa = 0.95, debug=FALSE) {
  
  archive <- ordering_fn(population[, 1:num_clusters], dmatrix_expression, dmatrix_biological)
  archive <- cbind( archive, explored=rep( FALSE,nrow(archive) ) )
  archive <- archive[archive$solutions_rank <= rank_cutoff, ]
  
  new_pool_size <- round(sqrt(sqrt(exploration_size)))
  step <- 0
  estimated_max_steps <- round(exploration_size * 1.2 / new_pool_size)
  evaluations <- 0
  row_name_id <- 0
  initial_temperature <- (0.1/alfa**(estimated_max_steps))
  cooling <- list() # DEBUGG
  energy <- list() # DEBUGG
  temperature_l <- list() # DEBUGG
  
  while(evaluations < exploration_size && step < exploration_size) {
    if (debug) {
      message(paste("Step: ", step, ". Evaluations: ", evaluations, "/", exploration_size, sep=""))
    }
    # Geometric cooling
    temperature <- (alfa**step) * initial_temperature
    
    row_index <- rownames( archive[ sample(nrow(archive), 1), ] )[1]
    solution <- archive[row_index, , drop=FALSE]
    
    # Generate neighborhood
    medoid_neighborhood <- helper.generate.neighborhood(new_pool_size, num_clusters, solution, neighborhood_matrix, row_name_id, debug)
    row_name_id <- row_name_id + new_pool_size + 1
    
    for(i in 1:nrow(medoid_neighborhood)) {
      if (evaluations < exploration_size) {
        fitness_helper <- fitness_counter
        
        objective_exp <- fitness_fn( medoid_neighborhood[i, 1:num_clusters], dmatrix_expression, 'expression' )
        objective_bio <- fitness_fn( medoid_neighborhood[i, 1:num_clusters], dmatrix_biological, 'biological' )
        
        evaluations <- helper.add.evaluation(fitness_helper, evaluations)
        
        cooling <- append(cooling, helper.mosa.probability(objective_exp, objective_bio, solution, temperature)) # DEBUGG
        energy <- append(energy, (0.2) - (solution$objective_exp - objective_exp) + (solution$objective_bio - objective_bio)) # DEBUGG 
        temperature_l <- append(temperature_l, temperature) # DEBUGG
        # If it passes the acceptance criteria or if it has enough probability
        if( acceptance_criteria_fn(objective_exp, objective_bio, solution) ||
            runif(1) < helper.mosa.probability(objective_exp, objective_bio, solution, temperature) ) {
          medoid_neighborhood[i, num_clusters+1] <- TRUE
        }
      }
    }
    
    if( sum(medoid_neighborhood[, ncol(medoid_neighborhood)]) > 1 ) {
      archive <- helper.cut.and.join.archive(
        archive, 
        medoid_neighborhood[ medoid_neighborhood$add == TRUE, 1:num_clusters ],
        gene_list, num_clusters, ordering_fn, rank_cutoff = rank_cutoff 
      )
    }
    
    step <- step + 1
  }
  archive <- helper.randomize.duplicates( archive, gene_list, num_clusters )
  archive <- ordering_fn(archive[, 1:num_clusters], dmatrix_expression, dmatrix_biological)
  
  if(debug) {
    data <- cbind(g = 1:length(energy), energy=unlist(energy), prob=unlist(cooling), temperature=(unlist(temperature_l) / max(unlist(temperature_l))))
    ggplot(as.data.frame(data), aes(x = g)) +
      geom_line(aes(y = energy), color = 'cyan') +
      geom_area(aes(y = energy), fill = 'cyan', alpha = .1) +
      geom_line(aes(y = prob), color="steelblue") +
      geom_area(aes(y = prob), fill = 'darkred', alpha = .1) +
      geom_line(aes(y = temperature), color="red") +
      geom_area(aes(y = temperature), fill = 'white', alpha = .1) +
      xlab('Generacion') +
      ylab('Jaccard\npromedio') +
      ggtitle('Similitud promedio entre todos las soluciones de la poblacion') + ylim(-0.5, 1.5) +
      theme(text = element_text(color = "#222222")
            ,panel.background = element_rect(fill = '#444B5A')
            ,panel.grid.minor = element_line(color = '#4d5566')
            ,panel.grid.major = element_line(color = '#586174')
            ,plot.title = element_text(size = 16)
            ,axis.title = element_text(size = 12, color = '#333333')
      )
    
    ggsave(paste('mosa-prob-energy-', round(stats::runif(1, 1, 10000)), '.png', sep=''), device="png", path="plots")
    message( paste("Multiobjective simulated annealing ran for", step, "iterations") )
  }
  return( archive )
}

#' Multiobjective clustering ensemble local search. Combines all pairs of solutions in
#' the pareto frontier. If the local search budget is lacking, it samples less elements
#' from each pair of solutions.
#' 
#' @param exploration_size Integer value. How many solutions should be explored each
#' iteration. Usually, the population_size of NSGA-II is used.
#' @param population Matrix / data.frame. Each row represents a solution.
#' @param num_clusters Integer value. How many medoids each solution has.
#' @param gene_list Vector. All the genes that can be used in the population. These
#' names should match the rows and cols of dmatrix_expression and dmatrix_biological. 
#' @param dmatrix_expression Gene expression distance matrix between each pair of genes.
#' @param dmatrix_biological Biological distance matrix between each pair of genes.
#' @param neighborhood_matrix Distance matrix. How close each pair of genes is based on
#' its expression and biological distance.
#' @param ordering_fn Function used to sort the population. 
#' Usually uses operator.nsga2.sorting.and.crowding .
#' @param fitness_fn Function used to calculate objectives. Usually uses fitness.medoid.wg
#' @param rank_cutoff Integer value. Solutions below this rank are removed from the archive.
#' 
#' @return Matrix of solutions created by combining the top solutions in the population.
#' 
local.search.ensemble <- function(exploration_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, ordering_fn, fitness_fn, rank_cutoff = 1, debug=FALSE) {
  
  archive <- ordering_fn(population[, 1:num_clusters], dmatrix_expression, dmatrix_biological)
  archive <- archive[archive$solutions_rank <= rank_cutoff, ]
  archive <- cbind( archive, explored=rep( FALSE,nrow(archive) ) )
  
  new_archive <- archive[1, 1:num_clusters, drop=FALSE]
  
  # Cant apply ensemble if pareto front is only one solution.
  if (nrow(new_archive) == 1) {
    return(new_archive)
  }
  
  evaluations <- 0
  step <- 0
  row_name_id <- 0
  
  if (exploration_size > (nrow(archive) * nrow(archive) / 2 * (2**num_clusters))) {
    # Combine each pair of solutions in the top rank_cutoff ranks.
    for(first_idx in 1:nrow(archive)) {
      for(second_idx in (first_idx + 1):nrow(archive)) {
        #message(paste("[", first_idx, ", ", second_idx, "]", sep=""))
        if( second_idx > nrow(archive) || second_idx == first_idx) next
        
        medoid_combinations <- helper.generate.ensemble(archive[first_idx, 1:num_clusters], archive[second_idx, 1:num_clusters], row_name_id)
        row_name_id <- row_name_id + (2**num_clusters) + 1
        
        for (i in 1:nrow(medoid_combinations)) {
          objective_exp <- fitness_fn( medoid_combinations[i, 1:num_clusters], dmatrix_expression, 'expression' )
          objective_bio <- fitness_fn( medoid_combinations[i, 1:num_clusters], dmatrix_biological, 'biological' )
          
          if( helper.non.dominated(objective_exp, objective_bio, archive[first_idx, ]) ||
              helper.non.dominated(objective_exp, objective_bio, archive[second_idx, ]) ) {
            
            # This process may find duplicates. At the end of this algorithm this is addressed.
            new_archive <- rbind(new_archive, medoid_combinations[i, 1:num_clusters, drop=FALSE])
          }
        }
      }
    }
  } else {
    new_pool_size <- sqrt(exploration_size)
    # Randomly choose pairs and combine them until exploration_size solutions are created.
    # Also, avoid entering an infinite loop when no new solutions are found for some iterations
    while (evaluations < exploration_size && step < exploration_size * 2) {
      
      indexes <- sample(1:nrow(archive), 2, replace = FALSE)
      
      medoid_combinations <- helper.generate.partial.ensemble(archive[indexes[1], 1:num_clusters], archive[indexes[2], 1:num_clusters], new_pool_size, row_name_id)
      row_name_id <- row_name_id + new_pool_size + 1
      
      for (i in 1:nrow(medoid_combinations)) {
        if (evaluations < exploration_size) {
          fitness_helper <- fitness_counter
          
          objective_exp <- fitness_fn( medoid_combinations[i, 1:num_clusters], dmatrix_expression, 'expression' )
          objective_bio <- fitness_fn( medoid_combinations[i, 1:num_clusters], dmatrix_biological, 'biological' )
          
          evaluations <- helper.add.evaluation(fitness_helper, evaluations)
          
          if( helper.non.dominated(objective_exp, objective_bio, archive[indexes[1], ]) ||
              helper.non.dominated(objective_exp, objective_bio, archive[indexes[2], ]) ) {
            
            new_archive <- rbind(new_archive, medoid_combinations[i, 1:num_clusters, drop=FALSE])
          }
        }
      }
      
      if (debug) {
        message(paste("Step: ", step, ". Evaluations: ", evaluations, "/", exploration_size, sep=""))
      }
      step <- step + 1
    }
  }
  
  
  # new_archive <- helper.randomize.duplicates(new_archive, gene_list, num_clusters)
  new_archive <- unique(new_archive)
  
  if (debug) {
    message( paste("Clustering ensemble ran for", step, "iterations") )
  }
  return( new_archive )
}
