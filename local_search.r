
#' Generates neighboring solutions (to a medoid representation) using a
#' pre-generated neighborhood matrix.
#' The generation method is as follows. From the medoid, 1 to num_clusters genes 
#' are randomly selected, each of these genes is replaced with a random neighboring
#' gene from the \code{neighborhood_matrix}. A gene is only replaced if the
#' replacement isn't currently in the solution.
#' 
#' @param population_size The amount of neighboring solutions to generate
#' @param num_clusters The number of clusters in a solution
#' @param solutions A matrix of one or more solutions. Could be the pareto front
#' or just a random group of solutions.
#' @param neighborhood_matrix An adjacency matrix with all neighboring genes.
#' @param row_name_id An unique identifier that will be used to identify each row.
#' Its not changed inside this method, so when it is used it should be manually
#' updated after running this method.
#' @return A matrix with neighboring solutions.
#' 
helper.generate.neighborhood <- function(population_size, num_clusters, solutions, neighborhood_matrix, row_name_id) {
  medoid_neighborhood <- as.data.frame( matrix(0, population_size, num_clusters) )
  for(medoid_generated in 1:nrow(medoid_neighborhood)) {
    medoid <- solutions[ sample(1:nrow(solutions), 1), ]
    
    # Random neighborhood
    # Change a random amount of genes in the solution to a random neighboring gene.
    # genes_to_change <- sample(colnames(medoid[, 1:num_clusters]), sample(1:num_clusters, 1))
    # for(gene_index in 1:length(genes_to_change)) {
    #   gene <- medoid[, genes_to_change[gene_index] ]
    #   if( length( neighborhood_matrix[[gene]] ) > 0 ) {
    #     neighbor_gene <- sample( neighborhood_matrix[[gene]], 1 )
    #     
    #     if( !(neighbor_gene %in% medoid[, 1:6]) ) {
    #       medoid[, gene_index] <- neighbor_gene
    #     }
    #   }
    # }
    # medoid_neighborhood[medoid_generated, ] <- medoid[, 1:num_clusters]
    
    # 1 change to neighborhood
    genes_with_neighbors <- medoid[, 1:num_clusters]
    # Find genes that have at least 1 neighbor
    gene <- sample(genes_with_neighbors, 1)
    neighbor_gene <- sample( neighborhood_matrix[[gene]], 1 )

    if( !(neighbor_gene %in% medoid[, 1:6]) ) {
      medoid[, gene_index] <- neighbor_gene
    }
    medoid_neighborhood[medoid_generated, ] <- medoid[, 1:num_clusters]
  }
  
  #print( paste("Neighbors found:", nrow(medoid_neighborhood), "without duplication")  )
  medoid_neighborhood <- medoid_neighborhood[ !duplicated(medoid_neighborhood[, 1:num_clusters]), ]
  #print( paste("Neighbors found:", nrow(medoid_neighborhood), "removing duplication")  )
  medoid_neighborhood <- cbind(medoid_neighborhood, add=rep( FALSE,nrow(medoid_neighborhood) ) )
  rownames(medoid_neighborhood) <- sapply(rownames(medoid_neighborhood), function(unused) {
    row_name_id <<- row_name_id + 1
    paste("V", row_name_id, sep="")
  })
  
  return(medoid_neighborhood)
}

helper.generate.ensemble <- function(first_solution, second_solution, row_name_id) {
  num_clusters <- ncol(first_solution)
  medoid_ensemble <- as.data.frame( matrix(0, (2**num_clusters) - 2, num_clusters) )
  
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
  #medoid_ensemble <- cbind(medoid_ensemble, add=rep( FALSE,nrow(medoid_ensemble) ) )
  rownames(medoid_ensemble) <- sapply(rownames(medoid_ensemble), function(unused) {
    row_name_id <<- row_name_id + 1
    paste("V", row_name_id, sep="")
  })
  
  return(medoid_ensemble)
}

helper.cut.and.join.archive <- function(archive, add_to_archive, gene_list, num_clusters, ordering_fn, rank_cutoff) {
  explored <- archive[ , "explored", drop=FALSE]
  
  new_archive <- helper.randomize.duplicates( rbind( archive[, 1:num_clusters], add_to_archive), gene_list, num_clusters )
  archive <- ordering_fn(new_archive, dmatrix_expression, dmatrix_biological)
  archive <- archive[archive$solutions_rank <= rank_cutoff, ]
  
  new_solutions <- rownames( archive[!( rownames(archive) %in% rownames(explored) ), ] )
  new_solutions <- as.data.frame(rep( FALSE,length(new_solutions) ))
  rownames(new_solutions) <- rownames( archive[!( rownames(archive) %in% rownames(explored) ), ] )
  colnames(new_solutions) <- "explored"
  
  explored <- explored[ rownames(explored) %in% rownames(archive) , , drop=FALSE]
  explored <- rbind( new_solutions, explored )
  
  explored <- explored[ order(row.names(explored)), , drop=FALSE ]
  archive <- archive[ order(row.names(archive)), , drop=FALSE ]
  archive <- cbind( archive, explored )
  
  return( archive )
}

helper.non.dominated <- function(objective_exp, objective_bio, solution) {
  return( objective_bio < solution$objective_bio || objective_exp < solution$objective_exp )
}

helper.dominates <- function(objective_exp, objective_bio, solution) {
  return( ( objective_bio <= solution$objective_bio && objective_exp < solution$objective_exp ) ||
          ( objective_bio < solution$objective_bio && objective_exp <= solution$objective_exp ) )
}

#' This became part of normal pls when rank_cutoff is equal to 1.
local.search.frontier.pareto.local.search <- function(population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, ordering_fn, fitness_fn) {
  
  archive <- ordering_fn(population[, 1:num_clusters], dmatrix_expression, dmatrix_biological)
  archive <- archive[archive$solutions_rank == 1, ]
  archive <- cbind( archive, explored=rep( FALSE,nrow(archive) ) )
    
  max_generations <- 50
  g <- 1
  row_name_id <- 0
  
  while( sum(archive[, ncol(archive)]) < nrow(archive) && g <= max_generations ) {
    archive[, ncol(archive)] <- TRUE
    
    # Generate neighborhood
    medoid_neighborhood <- helper.generate.neighborhood(population_size * nrow(archive) * 2, num_clusters, archive, neighborhood_matrix, row_name_id)
    row_name_id <- row_name_id + population_size * nrow(archive) * 2 + 1
    
    for(i in 1:nrow(medoid_neighborhood)) {
      objective_exp <- fitness_fn( medoid_neighborhood[i, 1:num_clusters], dmatrix_expression )
      objective_bio <- fitness_fn( medoid_neighborhood[i, 1:num_clusters], dmatrix_biological )
      
      if( sum( archive[ archive$objective_bio > objective_bio, "objective_bio" ] ) > 0 || 
          sum( archive[ archive$objective_exp > objective_exp, "objective_exp"  ] ) > 0 ) {
        medoid_neighborhood[i, num_clusters+1] <- TRUE
      }
    }
    
    # Add solutions found in the neighborhood if any non-dominated solution was found
    if( sum(medoid_neighborhood[, ncol(medoid_neighborhood)]) > 1 ) {
      archive <- helper.cut.and.join.archive(
        archive, 
        medoid_neighborhood[ medoid_neighborhood$add == TRUE, 1:num_clusters ],
        gene_list, num_clusters, ordering_fn, rank_cutoff = 1
      )
    }
    
    g <- g + 1
    
    #print(g)
    #print(archive)
  }
  
  print( paste("Frontier Pareto local search ran for", g, "iterations") )
  return( archive )
}

#' 
#' 
#' 
#' 
#' 
#' 
local.search.pareto.local.search <- function(population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, ordering_fn, fitness_fn, acceptance_criteria_fn = helper.non.dominated, rank_cutoff = 10, max_generations = 150) {
  
  archive <- ordering_fn(population[, 1:num_clusters], dmatrix_expression, dmatrix_biological)
  archive <- cbind( archive, explored=rep( FALSE,nrow(archive) ) )
  
  g <- 1
  row_name_id <- 0
  
  archive <- archive[archive$solutions_rank <= rank_cutoff, ]
  
  while( sum(archive[, ncol(archive)]) < nrow(archive) && g <= max_generations ) {
    unexplored <- archive[archive$explored == FALSE, , drop=FALSE]
    row_index <- rownames( unexplored[ sample(nrow(unexplored), 1), ] )[1]
    solution <- unexplored[row_index, , drop=FALSE]
    
    # Generate neighborhood
    medoid_neighborhood <- helper.generate.neighborhood(population_size, num_clusters, solution, neighborhood_matrix, row_name_id)
    row_name_id <- row_name_id + population_size + 1
    
    for(i in 1:nrow(medoid_neighborhood)) {
      objective_exp <- fitness_fn( medoid_neighborhood[i, 1:num_clusters], dmatrix_expression )
      objective_bio <- fitness_fn( medoid_neighborhood[i, 1:num_clusters], dmatrix_biological )
      
      if( acceptance_criteria_fn(objective_exp, objective_bio, solution) ) {
        medoid_neighborhood[i, num_clusters+1] <- TRUE
      }
    }
    
    # Add solutions found in the neighborhood if there are any that passed the
    # acceptance criteria used
    if( sum(medoid_neighborhood[, ncol(medoid_neighborhood)]) > 1 ) {
      archive <- helper.cut.and.join.archive(
        archive, 
        medoid_neighborhood[ medoid_neighborhood$add == TRUE, 1:num_clusters ],
        gene_list, num_clusters, ordering_fn, rank_cutoff = rank_cutoff
      )
    }
    
    
    if( row_index %in% rownames(archive) ) {
      archive[row_index, ncol(archive)] = TRUE
    }
    g <- g + 1
    
    print(g)
    print(archive)
  }
  
  print( paste("Pareto local search ran for", g, "iterations") )
  return( archive )
}

local.search.large.mols <- function(population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, ordering_fn, fitness_fn) {
  
  population <- ordering_fn(population[, 1:num_clusters], dmatrix_expression, dmatrix_biological)
  archive <- population[1, 1:num_clusters, drop=FALSE]
  
  row_name_id <- 0
  
  for(row_index in 1:nrow(population)) {
    
    solution <- population[row_index, , drop=FALSE]
    medoid_neighborhood <- helper.generate.neighborhood(population_size, num_clusters, solution, neighborhood_matrix, row_name_id)
    row_name_id <- row_name_id + population_size + 1
    
    for(i in 1:nrow(medoid_neighborhood)) {
      objective_exp <- fitness_fn( medoid_neighborhood[i, 1:num_clusters], dmatrix_expression )
      objective_bio <- fitness_fn( medoid_neighborhood[i, 1:num_clusters], dmatrix_biological )
      
      if( helper.non.dominated(objective_exp, objective_bio, solution) ) {
        archive <- rbind(archive, medoid_neighborhood[i, 1:num_clusters])
      }
    }
  }
  
  return( archive )
}

local.search.narrow.mols <- function(population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, ordering_fn, fitness_fn) {
  
  population <- ordering_fn(population[, 1:num_clusters], dmatrix_expression, dmatrix_biological)
  archive <- population[1, 1:num_clusters, drop=FALSE]
  
  row_name_id <- 0
  
  for(row_index in 1:nrow(population)) {
    
    solution <- population[row_index, , drop=FALSE]
    medoid_neighborhood <- helper.generate.neighborhood(population_size, num_clusters, solution, neighborhood_matrix, row_name_id)
    row_name_id <- row_name_id + population_size + 1
    
    for(i in 1:nrow(medoid_neighborhood)) {
      objective_exp <- fitness_fn( medoid_neighborhood[i, 1:num_clusters], dmatrix_expression )
      objective_bio <- fitness_fn( medoid_neighborhood[i, 1:num_clusters], dmatrix_biological )
      
      if( helper.dominates(objective_exp, objective_bio, solution) ) {
        archive <- rbind(archive, medoid_neighborhood[i, 1:num_clusters])
      }
    }
  }
  
  return( archive )
}

# Multiobjective simulated annealing
local.seach.mosa <- function(population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, ordering_fn, fitness_fn, acceptance_criteria_fn = helper.non.dominated, rank_cutoff = 10, min_temperature = 0.001, max_temperature = 1000) {
  
  archive <- ordering_fn(population[, 1:num_clusters], dmatrix_expression, dmatrix_biological)
  archive <- cbind( archive, explored=rep( FALSE,nrow(archive) ) )
  
  step <- 0
  max_steps <- 500
  row_name_id <- 0
  archive <- archive[archive$solutions_rank <= max_rank, ]
  
  while( temperature > 0 ) {
    temperature <- (step + 1) / max_steps
    
    row_index <- rownames( archive[ sample(nrow(archive), 1), ] )[1]
    solution <- archive[row_index, , drop=FALSE]
    
    # Generate neighborhood
    medoid_neighborhood <- helper.generate.neighborhood(population_size, num_clusters, solution, neighborhood_matrix, row_name_id)
    row_name_id <- row_name_id + population_size + 1
    
    for(i in 1:nrow(medoid_neighborhood)) {
      objective_exp <- fitness_fn( medoid_neighborhood[i, 1:num_clusters], dmatrix_expression )
      objective_bio <- fitness_fn( medoid_neighborhood[i, 1:num_clusters], dmatrix_biological )
      
      # If it passes the acceptance criteria or if it 
      if( acceptance_criteria_fn(objective_exp, objective_bio, solution) ||
          TRUE) {
        medoid_neighborhood[i, num_clusters+1] <- TRUE
      }
    }
    
    if( sum(medoid_neighborhood[, ncol(medoid_neighborhood)]) > 1 ) {
      archive <- helper.cut.and.join.archive(
        archive, 
        medoid_neighborhood[ medoid_neighborhood$add == TRUE, 1:num_clusters ],
        gene_list, num_clusters, ordering_fn, rank_cutoff = rank_cutoff 
      )
    }
    
    
    if( row_index %in% rownames(archive) ) {
      archive[row_index, ncol(archive)] = TRUE
    }
    step <- step + 1
    
    #print(paste(step, temperature, sep=" - "))
    #print(archive)
  }
  
  print( paste("Multiobjective simulated annealing ran for", step, "iterations") )
  return( archive )
}

# Multiobjective clustering ensemble local search
local.search.ensemble <- function(population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, ordering_fn, fitness_fn, rank_cutoff = 1) {
  
  archive <- ordering_fn(population[, 1:num_clusters], dmatrix_expression, dmatrix_biological)
  archive <- archive[archive$solutions_rank <= rank_cutoff, ]
  archive <- cbind( archive, explored=rep( FALSE,nrow(archive) ) )
  
  new_archive <- archive[1, 1:num_clusters, drop=FALSE]
  
  step <- 0
  row_name_id <- 0
  
  for(first_idx in 1:nrow(archive)) {
    for(second_idx in (first_idx + 1):nrow(archive)) {
      #print(paste("[", first_idx, ", ", second_idx, "]", sep=""))
      if( second_idx > nrow(archive) || second_idx == first_idx) next
      
      medoid_combinations <- helper.generate.ensemble(archive[first_idx, 1:num_clusters], archive[second_idx, 1:num_clusters], row_name_id)
      row_name_id <- row_name_id + (2**num_clusters) + 1
      
      for (i in 1:nrow(medoid_combinations)) {
        objective_exp <- fitness_fn( medoid_combinations[i, 1:num_clusters], dmatrix_expression )
        objective_bio <- fitness_fn( medoid_combinations[i, 1:num_clusters], dmatrix_biological )
        
        if( helper.non.dominated(objective_exp, objective_bio, archive[first_idx, ]) ||
            helper.non.dominated(objective_exp, objective_bio, archive[second_idx, ]) ) {
          
          # This process may find duplicates. At the end of this algorithm this is addressed.
          new_archive <- rbind(new_archive, medoid_combinations[i, 1:num_clusters, drop=FALSE])
        }
      }
      
      step <- step + 1
      #print(step)
    }
  }
  
  new_archive <- helper.randomize.duplicates(new_archive, gene_list, num_clusters)
  
  print( paste("Clustering ensemble ran for", step, "iterations") )
  return( new_archive )
}