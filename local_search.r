
local.search.pareto.local.search <- function(population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, ordering_fn, fitness_fn, neighborhood) {
  
  archive <- ordering_fn(population[, 1:num_clusters], dmatrix_expression, dmatrix_biological)
  archive <- archive[archive$solutions_rank == 1, ]
  archive <- cbind( archive, explored=rep( FALSE,nrow(archive) ) )
    
  max_generations <- 100
  g <- 1
  row_name_id <- 0
  
  while( sum(archive[, ncol(archive)]) < nrow(archive) && g <= max_generations ) {
    unexplored <- archive[archive$explored == FALSE, ]
    row_index <- rownames( unexplored[ sample(nrow(unexplored), 1), ] )[1]
    
    # Generate neighborhood
    
    medoid_neighborhood <- as.data.frame( matrix(0, population_size, num_clusters) )
    for(medoid_generated in 1:nrow(medoid_neighborhood)) {
      medoid <- archive[ sample(nrow(archive), 1), ]
      
      genes_to_change <- sample(colnames(medoid[, 1:num_clusters]), num_clusters-1)
      for(gene_index in 1:length(genes_to_change)) {
        gene <- medoid[, gene_index]
        if( length( unlist(neighborhood_matrix[gene]) ) > 0 ) {
          neighbor_gene <- sample( unlist(neighborhood_matrix[gene]), 1 )
          
          if( !(neighbor_gene %in% medoid[, 1:6]) ) {
            medoid[, gene_index] <- neighbor_gene
          }
        }
      }
      medoid_neighborhood[medoid_generated, ] <- medoid[, 1:num_clusters]
    }
    
    medoid_neighborhood <- medoid_neighborhood[ !duplicated(medoid_neighborhood[, 1:num_clusters]), ]
    medoid_neighborhood <- cbind(medoid_neighborhood, add=rep( FALSE,nrow(medoid_neighborhood) ) )
    rownames(medoid_neighborhood) <- sapply(rownames(medoid_neighborhood), function(unused) {
      row_name_id <<- row_name_id + 1
      paste("V", row_name_id, sep="")
    })
    
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
      add_archive <- medoid_neighborhood[ medoid_neighborhood$add == TRUE, 1:num_clusters ]
      explored <- archive[ , "explored", drop=FALSE]
      
      new_archive <- helper.randomize.duplicates( rbind( archive[, 1:num_clusters], add_archive), gene_list, num_clusters )
      archive <- ordering_fn(new_archive, dmatrix_expression, dmatrix_biological)
      archive <- archive[archive$solutions_rank == 1, ]
      
      new_solutions <- rownames( archive[!( rownames(archive) %in% rownames(explored) ), ] )
      new_solutions <- as.data.frame(rep( FALSE,length(new_solutions) ))
      rownames(new_solutions) <- rownames( archive[!( rownames(archive) %in% rownames(explored) ), ] )
      colnames(new_solutions) <- "explored"
      
      explored <- explored[ rownames(explored) %in% rownames(archive) , , drop=FALSE]
      explored <- rbind(  new_solutions, explored )
      
      explored <- explored[ order(row.names(explored)), , drop=FALSE ]
      archive <- archive[ order(row.names(archive)), , drop=FALSE ]
      archive <- cbind( archive, explored )
    }
      
    
    if( row_index %in% rownames(archive) ) {
      archive[row_index, ncol(archive)] = TRUE
    }
    g <- g + 1
  }
  
  print( paste("Pareto local search ran for", g, "iterations") )
  rownames(archive) <- 1:nrow(archive)
  return( archive )
}

local.search.path.relinking <- function(population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, ordering_fn, fitness_fn, neighborhood) {
  
}

local.search.large.mols <- function(population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, ordering_fn, fitness_fn, neighborhood) {
  
}

local.search.narrow.mols <- function(population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, ordering_fn, fitness_fn, neighborhood) {
  
}

# Multiobjective simulated annealing
local.seach.mosa <- function(population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, ordering_fn, fitness_fn, neighborhood) {
  
}

# Multiobjective clustering ensemble local search
local.search.ensemble <- function(population_size, population, num_clusters, gene_list, dmatrix_expression, dmatrix_biological, neighborhood_matrix, ordering_fn, fitness_fn, neighborhood) {
  
}