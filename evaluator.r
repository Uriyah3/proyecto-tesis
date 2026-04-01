library(amap)
library(mco)
library(cluster)
library(stringr)
library(stats)
source("nsga2.r")
source("david_client.r")

#' Normalize data to [0, 1] range using min-max normalization
#' 
#' @param data Numeric vector to normalize
#' @return Normalized vector with values between 0 and 1
#' 
helper.normalize <- function(data) {
  return( (data - min(data)) / (max(data) - min(data)) )
}

#' Find indices of median values in a vector
#' 
#' Similar to which.max/which.min but for median. Handles both odd and even
#' length vectors.
#' 
#' @param x Numeric vector
#' @return Index or indices of median value(s)
#' @note Code from https://stackoverflow.com/questions/10256503/
#' 
which.median = function(x) {
  if (length(x) %% 2 != 0) {
    which(x == median(x))
  } else if (length(x) %% 2 == 0) {
    a = sort(x)[c(length(x)/2, length(x)/2+1)]
    c(which(x == a[1]), which(x == a[2]))
  }
}

#' Choose the solution that has silhouette higher than or equal to the median
#' and minimum biological optimization index.
which.best <- function(which.x, silhouette_results, population = NULL) {
  if(is.character(which.x) && which.x == 'best'){ 
    median_silhouette <- silhouette_results[which.median(silhouette_results)[1]]
    min_bio = 2
    best_solution_index = 0
    for(i in 1:length(silhouette_results)) {
      if(silhouette_results[i] >= median_silhouette) {
        if(population$objective_bio[i] < min_bio ) {
          best_solution_index <- i
          min_bio <- population$objective_bio[i]
        }
      }
    }
    return(best_solution_index)
  } else {
    return(which.x(silhouette_results)[1])
  }
}

#' Evaluate multi-objective clustering results with biological significance
#' 
#' Calculates silhouette, hypervolume, and biological enrichment metrics.
#' Selects best solution based on median silhouette and minimum biological index.
#' 
#' @param results List with $population and $clustering from NSGA-II
#' @param dmatrix Distance matrix for silhouette calculation
#' @param debug Boolean for debug output
#' @param which.x Function to select solution (which.max, which.min, or 'best')
#' @param dataset_name String dataset identifier for caching
#' @param bio String biological source identifier for caching
#' @param iter Integer iteration number for caching
#' @return List with silhouette, hypervolume, biological metrics
#' 
evaluator.multiobjective.clustering <- function( results, dmatrix, debug = FALSE, which.x = which.max, dataset_name=NULL, bio=NULL, iter=NULL ) {
  
  silhouette_results <- evaluator.silhouette( results$clustering, dmatrix, debug=debug, dataset_name=dataset_name, bio=bio, iter=iter)
  hypervolume_results <- evaluator.hypervolume( results$population, debug=debug, dataset_name=dataset_name, bio=bio, iter=iter )
  
  # cond_sil <- which.x(silhouette_results$silhouette)[1]
  cond_sil <- which.best(which.x, silhouette_results$silhouette, results$population)
  if (debug) {
    message(paste("Se analiza biologicamente solucion con silueta =", silhouette_results$silhouette[[cond_sil]]))
  }
  biology_results <- evaluator.biological.significance( results$clustering[[cond_sil]], colnames(dmatrix), dataset_name=dataset_name, bio=bio, debug=debug, iter=iter )
  biology_summary <- unlist(biology_results)
  biology_summary <- c(by(biology_summary, names(biology_summary), mean, na.rm = TRUE))
  
  metrics <- list(
    silhouette = silhouette_results,
    hypervolume = hypervolume_results,
    biological = biology_results,
    biological_summary = biology_summary
  )
  return( metrics )
}

#' Evaluate multi-objective clustering without biological validation
#' 
#' Faster version that only calculates silhouette and hypervolume metrics.
#' Used during algorithm execution when biological validation is too slow.
#' 
#' @param results List with $population and $clustering from NSGA-II
#' @param dmatrix Distance matrix for silhouette calculation
#' @param debug Boolean for debug output
#' @param dataset_name String dataset identifier for caching
#' @param bio String biological source identifier for caching
#' @param iter Integer iteration number for caching
#' @param which.x Unused parameter (kept for compatibility)
#' @return List with silhouette and hypervolume metrics
#' 
evaluator.multiobjective.clustering.no.bio <- function( results, dmatrix, debug = FALSE, dataset_name=NULL, bio=NULL, iter=NULL, which.x=NULL ) {
  
  silhouette_results <- evaluator.silhouette( results$clustering, dmatrix, debug=debug, dataset_name=dataset_name, bio=bio, iter=iter )
  hypervolume_results <- evaluator.hypervolume( results$population, debug=debug, dataset_name=dataset_name, bio=bio, iter=iter )
  
  metrics <- list(
    silhouette = silhouette_results,
    hypervolume = hypervolume_results
  )
  
  return( metrics )
}

#' Calculate hypervolume indicator for multi-objective optimization
#' 
#' Computes raw, normalized, and centered hypervolume of Pareto front.
#' Results are cached to avoid redundant calculations.
#' 
#' @param population Matrix with objective columns
#' @param debug Boolean for debug output
#' @param dataset_name String dataset identifier for caching
#' @param bio String biological source identifier for caching
#' @param iter Integer iteration number for caching
#' @return List with hypervolume, normalized_hypervolume, centered_hypervolume
#' 
evaluator.hypervolume <- function( population, debug = FALSE, dataset_name=NULL, bio=NULL, iter=NULL ) {
  if (!is.null(metrics <- load.evaluation.from.cache(dataset_name, bio, iter, 'hypervolume'))) {
    return(metrics)
  }
  num_clusters <- ncol(population) - 4
  
  pareto_points <- as.matrix( population[ , (num_clusters+1):(num_clusters+2) ] )
  hypervolume <- dominatedHypervolume(pareto_points)
  centered_hypervolume <- dominatedHypervolume(pareto_points, c(1, 1))
  
  pareto_points[, 1] <- helper.normalize(pareto_points[, 1])
  pareto_points[, 2] <- helper.normalize(pareto_points[, 2])
  n_hypervolume <- dominatedHypervolume(pareto_points, c(1, 1))
  
  metrics <- list(
    hypervolume = hypervolume, # Raw hypervolume
    normalized_hypervolume = n_hypervolume, # Hypervolume using normalized data to [0,1].
    centered_hypervolume = centered_hypervolume # Hypervolume with (1,1) as reference point
  )
  store.evaluation.to.cache(metrics, dataset_name, bio, iter, 'hypervolume')
  return( metrics )
}

#' Calculate silhouette coefficient for clustering solutions
#' 
#' Evaluates cluster quality for all solutions in the Pareto front.
#' Higher silhouette values indicate better-defined clusters.
#' 
#' @param clustering List of cluster assignment vectors
#' @param dmatrix Distance matrix between genes
#' @param debug Boolean for debug output
#' @param dataset_name String dataset identifier for caching
#' @param bio String biological source identifier for caching
#' @param iter Integer iteration number for caching
#' @return List with silhouette statistics (max, mean, min, sd)
#' 
evaluator.silhouette <- function( clustering, dmatrix, debug = FALSE, dataset_name=NULL, bio=NULL, iter=NULL ) {
  if (!is.null(metrics <- load.evaluation.from.cache(dataset_name, bio, iter, 'silhouette'))) {
    return(metrics)
  }
  
  solution_count <- length(clustering)
  
  silhouette_indices = double(solution_count)
  for (c_index in 1:solution_count ) {
    if (is.data.table(dmatrix)) {
      dmatrix_numeric <- as.matrix(dmatrix[, !names(dmatrix) %in% "rn", with=FALSE])
    } else {
      dmatrix_numeric <- as.matrix(dmatrix)
    }
    sil <- tryCatch(
      silhouette( clustering[[c_index]], dmatrix = dmatrix_numeric ),
      error = function(e) NULL
    )
    if (!is.null(sil) && inherits(sil, "silhouette")) {
      sil_summary <- summary(sil)
      silhouette_indices[c_index] <- sil_summary$avg.width
    } else {
      silhouette_indices[c_index] <- NA_real_
    }
  }
  
  metrics <- list(
    silhouette = silhouette_indices,
    max_silhouette = max(silhouette_indices),
    mean_silhouette = mean(silhouette_indices),
    min_silhouette = min(silhouette_indices),
    sd_silhouette = sd(silhouette_indices)
  )
  store.evaluation.to.cache(metrics, dataset_name, bio, iter, 'silhouette')
  return( metrics )
}

evaluator.biological.significance <- function( clustering, full_gene_list, dataset_name=NULL, bio=NULL, debug = FALSE, iter=NULL, id=NULL ) {
  if (!is.null(metrics <- load.evaluation.from.cache(dataset_name, bio, iter, 'biological'))) {
    return(metrics)
  }
  message("Calculating manually")
  # Agregar los nombres de los genes a la cabecera
  clustering <- as.data.frame(t(clustering))
  colnames(clustering) <- full_gene_list
  
  num_clusters <- max(clustering)
  
  results <- list()
  for (cluster in 1:num_clusters) {
    gene_list <- colnames(clustering[ , clustering == cluster, drop=FALSE ])
    
    if (debug) {
      message(paste("Procesando lista#", id, cluster, " con ", length(gene_list), " genes...", sep=""))
    }
    if(length(gene_list) < 3) {
      if (debug) {
        message(paste("Saltando el listado 1-2 genes:",paste0(gene_list, collapse=",")))
      }
      next
    }
    
    # For evaluation purposes, if any cluster has more elements than the 3000 DAVID 
    # limit. Apply kmeans over the genes and average the results
    if (length(gene_list) >= 3000) {
      if (!is.null(dataset_name) && !is.null(bio)) {
        if(bio == 'base') {
          bio.source = 'go'
        } else {
          bio.source = bio
        }
        
        if (bio.source %in% names(biological_databases)) {
          dmatrix <- biological.matrix(NULL, biological_databases[[bio.source]], dataset=dataset_name)
        } else  {
          dmatrix <- expression.matrix(NULL, dataset=dataset_name)
        }
        dmatrix <- dmatrix[rownames(dmatrix) %in% gene_list, colnames(dmatrix) %in% gene_list, drop=FALSE]
        
        intra_clustering <- kmeans(dmatrix, 5, iter.max=50, nstart=10)
        if(max(intra_clustering$size) >= length(gene_list) * 0.95 && max(intra_clustering$size) > 3000)  {
          message(str_interp("Run k-means with more clusters because the big group is not getting separated: ${length(gene_list)} -> ${max(intra_clustering$size)}"))
          
          intra_clustering <- kmeans(dmatrix, 10, iter.max=20, nstart=5)
          if(max(intra_clustering$size) >= length(gene_list) * 0.95 && max(intra_clustering$size) > 3000)  {
            message("Run k-means using expression matrix because bio matrix is not separating at all")
            dmatrix <- expression.matrix(NULL, dataset=dataset_name)
            dmatrix <- dmatrix[rownames(dmatrix) %in% gene_list, colnames(dmatrix) %in% gene_list, drop=FALSE]
            intra_clustering <- kmeans(dmatrix, 5, iter.max=50, nstart=10)
          }
        }
        
        temp_results <- evaluator.biological.significance(intra_clustering$cluster, gene_list, dataset_name, bio, debug, id=paste(id,cluster,'.',sep=""))
        enrichment <- unlist(lapply(temp_results, function(result) {
          result$cluster_count <- NULL
          return(unlist(result))
        }))
        results[[cluster]] <- list()
        results[[cluster]]$cluster_count <- sum( unlist(lapply(temp_results, '[[', 'cluster_count')) )
        results[[cluster]]$enrichment = enrichment
        
      } else {
        if (debug) {
          message(paste("No se puede utilizar DAVID con esta lista de genes, dado que son", length(gene_list), "genes"))
        }
        results[[cluster]] <- evaluator.biological.anotate.list( gene_list, debug )
      }
    } else {
      results[[cluster]] <- NA
      attempt <- 1
      while( identical(results[[cluster]], NA) && attempt <= 5 ) {
        attempt <- attempt + 1
        try(
          results[[cluster]] <- evaluator.biological.anotate.list( gene_list, debug )
        )
      }
    }
    
    if (!is.list(results[[cluster]]) ||
        (!is.list(results[[cluster]]$enrichment) && all(is.na(results[[cluster]]$enrichment)))) {
      results[[cluster]] <- NULL
      next
    }
    
    if (is.null(id)) {
      results[[cluster]]$max_enrichment <- max(results[[cluster]]$enrichment)
      results[[cluster]]$mean_enrichment <- mean(results[[cluster]]$enrichment)
      results[[cluster]]$min_enrichment <- min(results[[cluster]]$enrichment)
      results[[cluster]]$sd_enrichment <- sd(results[[cluster]]$enrichment)
    }
  }
  
  if (is.null(id)) {
    # Only iterate over numeric-indexed cluster results (not summary entries)
    cluster_indices <- which(sapply(results, is.list))
    enrichment <- unlist(lapply(results[cluster_indices], function(result) {
      return(unlist(result$enrichment))
    }))
    results$cluster_count <- sum( unlist(lapply(results[cluster_indices], '[[', 'cluster_count')) )
    results$max_enrichment <- max(enrichment)
    results$mean_enrichment <- mean(enrichment)
    results$min_enrichment <- min(enrichment)
    results$sd_enrichment <- sd(enrichment)
    
    store.evaluation.to.cache(results, dataset_name, bio, iter, 'biological')
  }
  
  return( results )
}

source("credentials.r")

#' Annotate gene list using DAVID web service
#'
#' Performs functional annotation and enrichment analysis via DAVID API.
#' Returns cluster count and enrichment scores.
#'
#' @param gene_list Vector of ENTREZ gene IDs
#' @param debug Boolean for debug output
#' @return List with cluster_count and enrichment values
#'
evaluator.biological.anotate.list <- function( gene_list, debug = FALSE ) {
  david_email <- get.david.email()
  if(debug) {
    message(paste("Using DAVID email:", david_email))
  }

  if (length(gene_list) >= 3000) {
    return(
      list(
        cluster_count = min(length(gene_list), 100),
        enrichment = list(0)
      )
    )
  }

  result <- david.annotate.genes(david_email, gene_list, debug = debug)
  return(result)
}

#' Run metaheuristic algorithm multiple times and aggregate results
#' 
#' Executes the specified metaheuristic (usually nsga2.custom) for multiple
#' runs and computes mean metrics across all runs.
#' 
#' @param metaheuristic Function to execute (e.g., nsga2.custom)
#' @param meta_params List of parameters to pass to metaheuristic
#' @param run_evaluator Function to evaluate results (default: no bio)
#' @param runs Integer number of independent runs
#' @param debug Boolean for debug output
#' @param dataset_name String dataset identifier
#' @param bio String biological source identifier
#' @return List with aggregated metrics and individual run results
#' 
evaluator.metaheuristics <- function(metaheuristic, meta_params, run_evaluator = evaluator.multiobjective.clustering.no.bio, runs = 13, debug = FALSE, dataset_name = NULL, bio=NULL) {
  
  results <- lapply( 1:runs, function(n) {
    if (debug) {
      message("--------------------------------------------------")
      message(paste("Running metaheuristic on iteration:", n))
      message("\n\n\n")
    }
    
    start.time <- Sys.time()
    
    metaheuristic_results <- do.call(metaheuristic, meta_params)
    iteration_results = run_evaluator(metaheuristic_results, meta_params$dmatrix_expression, debug=debug, dataset_name=dataset_name, bio=bio)
    
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    iteration_results$time <- time.taken
    iteration_results$fitness_counter <- fitness_counter
    
    # Delete data that can't be averaged at the end
    iteration_results$silhouette$silhouette <- NULL
    iteration_results
  })
  
  full_results <- results
  # Average only numeric metrics across runs (skip difftime, NULL, etc.)
  mean_results <- tryCatch({
    unlisted <- unlist(lapply(results, function(r) {
      r$time <- as.numeric(r$time)
      r
    }))
    c(by(unlisted, names(unlisted), function(x) mean(as.numeric(x), na.rm = TRUE)))
  }, error = function(e) {
    warning(paste("Could not compute mean_results:", e$message))
    NA
  })
  return( list(
    full_results = full_results,
    mean_results = mean_results
  ))
}

#' Save metaheuristic results to cache for later analysis
#' 
#' Stores population and clustering results from each run to .rda files.
#' 
#' @param dataset.name String dataset identifier
#' @param identifier String biological source or algorithm variant
#' @param metaheuristic Function that was executed
#' @param meta_params List of parameters used
#' @param runs Integer number of runs to save
#' @param debug Boolean for debug output
#' @return NULL (saves to files)
#' 
save.metaheuristic.results <- function(dataset.name, identifier, metaheuristic, meta_params, runs = 13, debug = FALSE) {
  lapply( 1:runs, function(n) {
    if (debug) {
      message("--------------------------------------------------")
      message(paste("Running metaheuristic on iteration:", n))
      message("\n\n\n")
    }
    
    start.time <- Sys.time()
    
    iteration_results <- list()
    iteration_results$nsga <- do.call(metaheuristic, meta_params)
    
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    iteration_results$time <- time.taken
    iteration_results$fitness_counter <- fitness_counter
    
    filename <- build.saved.results.filename(dataset.name, identifier, n)
    saveRDS(iteration_results, filename)
  })
}

#' Load and re-evaluate saved metaheuristic results
#' 
#' Reads previously saved results and applies evaluation metrics.
#' Useful for post-hoc analysis without re-running experiments.
#' 
#' @param dataset.name String dataset identifier
#' @param identifier String biological source or algorithm variant
#' @param run_evaluator Function to evaluate results
#' @param runs Integer number of runs to load
#' @param debug Boolean for debug output
#' @param skip.loading Boolean to skip loading files (use existing)
#' @return List with aggregated metrics
#' 
reconstruct.metaheuristic.saved.results <- function(dataset.name, identifier, run_evaluator = evaluator.multiobjective.clustering.no.bio, runs = 13, debug = FALSE, skip.loading=FALSE) {
  if (!skip.loading) {
    dmatrix_expression <- expression.matrix(NULL, dataset=dataset.name)
  }
  results <- lapply( 1:runs, function(n) {
    if (debug) {
      message("--------------------------------------------------")
      message(paste("Loading metaheuristic results on iteration:", n))
      message("\n")
    }
    
    if(identical(run_evaluator, evaluator.multiobjective.clustering)) {
      n <- find.best.solution.for.david(dataset.name, identifier)
      if(debug) {
        message(str_interp("Using execution ${n} to load biological results"))
      }
    }
    filename <- build.saved.results.filename(dataset.name, identifier, n)
    iteration_results <- readRDS(filename)
    results = run_evaluator(iteration_results$nsga, dmatrix_expression, which.x = 'best', debug=debug, dataset_name=dataset.name, bio=identifier, iter=n)
    iteration_results$nsga <- NULL
    iteration_results <- c(results, iteration_results)
    
    iteration_results$silhouette$silhouette <- NULL
    iteration_results
  })
  
  full_results <- results
  
  mean_results <- unlist(results)
  mean_results <- c(by(mean_results, names(mean_results), mean, na.rm = TRUE))
  return( list(
    full_results = full_results,
    mean_results = mean_results
  ))
}
  
build.saved.results.filename <- function(dataset.name, identifier, iteration) {
  return(
    str_interp("cache/metaheuristic-${dataset.name}-${identifier}-${iteration}.rda")
  )
}

build.evaluation.cache.filename <- function(dataset.name, identifier, iteration, evaluation) {
  return(
    str_interp("cache/metaheuristic-${dataset.name}-${identifier}-${iteration}-${evaluation}.rda")
  )
}

load.evaluation.from.cache <- function(dataset.name, identifier, iteration, evaluation) {
  file.name <- build.evaluation.cache.filename(dataset.name, identifier, iteration, evaluation)
  if( file.exists(file.name) ) {
    return( readRDS(file.name) )
  } else {
    return( NULL )
  }
}

store.evaluation.to.cache <- function(data, dataset.name, identifier, iteration, evaluation) {
  file.name <- build.evaluation.cache.filename(dataset.name, identifier, iteration, evaluation)
  saveRDS(data, file.name)
}

#' Find solution with best silhouette for DAVID analysis
#' 
#' Loads cached results and identifies the solution with highest silhouette
#' coefficient for biological validation.
#' 
#' @param dataset.name String dataset identifier
#' @param identifier String biological source identifier
#' @param evaluation String evaluation type (default: 'biological')
#' @param runs Integer number of runs to search
#' @return List with run number and solution index
#' 
find.best.solution.for.david <- function(dataset.name, identifier, evaluation = 'biological', runs=13) {
  dmatrix_expression <- NULL
  solutions <- lapply(1:runs, function(iteration) {
    filename <- build.saved.results.filename(dataset.name, identifier, iteration)
    iteration_results <- readRDS(filename)
    
    silhouette_results <- NULL
    try(
      silhouette_results <- evaluator.silhouette( iteration_results$nsga$clustering, NULL, dataset_name=dataset.name, bio=identifier, iter=iteration )
    )
    if (is.null(silhouette_results) && is.null(dmatrix_expression)) {
      # esta línea no funciona porque no se ha definido el objeto dataset
      dmatrix_expression <- expression.matrix(NULL, dataset=dataset$name)
      silhouette_results <- evaluator.silhouette( iteration_results$nsga$clustering, dmatrix_expression, dataset_name=dataset.name, bio=identifier, iter=iteration )
    }
    
    solution_index <- which.best('best', silhouette_results$silhouette, iteration_results$nsga$population)
    return(list(
      objective_bio = iteration_results$nsga$population$objective_bio[[solution_index]],
      silhouette = silhouette_results$silhouette[[solution_index]]
    ))
  })
  solutions <- do.call(rbind, lapply(solutions, data.frame))
  #solutions$silhouette <- -solutions$silhouette
  #ranking <- nsga2R::fastNonDominatedSorting( solutions )
  #solutions_rank <- helper.pareto.ranking( nrow(solutions), ranking )
  #solutions <- cbind(solutions, solutions_rank)
  #solutions$silhouette <- -solutions$silhouette
  
  #solutions <- solutions[order(solutions$solutions_rank), ]
  
  best_solution <- which.best('best', solutions$silhouette, solutions)
}

#results <- evaluator.metaheuristics(nsga2.custom, list(dmatrix_expression=dmatrix_expression, dmatrix_biological=dmatrix_biological, population_size = 2, generations = 1, num_clusters = 3, ls_pos=NULL, local_search = NULL, debug=TRUE))
#results <- nsga2.custom(dmatrix_expression, dmatrix_biological, population_size = 80, generations = 10, num_clusters = 5, ls_pos=NULL, local_search = NULL, debug=TRUE)

