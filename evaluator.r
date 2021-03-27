library(amap)
library(mco)
library(cluster)
library(RDAVIDWebService)
source("nsga2.r")

helper.normalize <- function(data) {
  return( (data - min(data)) / (max(data) - min(data)) )
}

evaluator.multiobjective.clustering <- function( results, dmatrix ) {
  
  silhouette_results <- evaluator.silhouette( results$clustering, dmatrix )
  hypervolume_results <- evaluator.hypervolume( results$population )
  
  best_sil <- which.max(silhouette_results$silhouette)
  biology_results <- evaluator.biological.significance( results$clustering[[best_sil]], colnames(dmatrix) )
  
  metrics <- list(
    silhouette = silhouette_results,
    hypervolume = hypervolume_results,
    biological = biology_results
  )
  return( metrics )
}

evaluator.multiobjective.clustering.no.bio <- function( results, dmatrix ) {
  silhouette_results <- evaluator.silhouette( results$clustering, dmatrix )
  hypervolume_results <- evaluator.hypervolume( results$population )
  
  metrics <- list(
    silhouette = silhouette_results,
    hypervolume = hypervolume_results
  )
  return( metrics )
}

evaluator.hypervolume <- function( population ) {
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
  return( metrics )
}

evaluator.silhouette <- function( clustering, dmatrix ) {
  solution_count <- length(clustering)
  
  silhouette_indices = double(solution_count)
  for (c_index in 1:solution_count ) {
    sil <- silhouette( clustering[[c_index]], dmatrix = as.matrix(dmatrix) )
    silhouette_indices[c_index] <- summary(sil)$avg.width
  }
  
  metrics <- list(
    silhouette = silhouette_indices,
    max_silhouette = max(silhouette_indices),
    mean_silhouette = mean(silhouette_indices),
    min_silhouette = min(silhouette_indices),
    sd_silhouette = sd(silhouette_indices)
  )
  return( metrics )
}

evaluator.biological.significance <- function( clustering, full_gene_list ) {
  # Agregar los nombres de los genes a la cabecera
  clustering <- as.data.frame(t(clustering))
  colnames(clustering) <- full_gene_list
  
  num_clusters <- max(clustering)
  
  results <- list()
  for (cluster in 1:num_clusters) {
    gene_list <- colnames(clustering[ , clustering == cluster ])
    
    results[[cluster]] <- evaluator.biological.anotate.list( gene_list )
  }
  
  return( results )
}

evaluator.biological.anotate.list <- function( gene_list ) {
  #https://david.ncifcrf.gov/webservice/services/DAVIDWebService/authenticate?args0=nicolas.mariangel@usach.cl
  david <- DAVIDWebService$new(email='nicolas.mariangel@usach.cl', url='https://david.ncifcrf.gov/webservice/services/DAVIDWebService')
  # Se cae a veces con el siguiente error: 
  # [INFO] Unable to sendViaPost to url[https://david.ncifcrf.gov/webservice/services/DAVIDWebService]
  # java.net.SocketTimeoutException: Read timed out
  # Considerar que también tira este error cuando la lista de genes es > 3000
  setTimeOut(david, 120000)
  
  if (length(gene_list) >= 3000) {
    return(
      list(
        cluster_count = 100,
        max_enrichment = 0,
        mean_enrichment = 0,
        min_enrichment = 0,
        sd_enrichment = 0
      )
    )
  }
  
  
  david$addList(gene_list, "ENTREZ_GENE_ID", listName = "Prueba de anotación chart", listType = "Gene")
  davidCluster <- david$getClusterReport()
  
  clusters <- davidCluster@cluster
  enrichment <- sapply(clusters, `[[`, 'EnrichmentScore')
  
  if (length(enrichment) == 0) {
    metrics <- list(
      cluster_count = 100,
      max_enrichment = 0,
      mean_enrichment = 0,
      min_enrichment = 0,
      sd_enrichment = 0
    )
  } else {
    metrics <- list(
      cluster_count = length(clusters),
      max_enrichment = max(enrichment),
      mean_enrichment = mean(enrichment),
      min_enrichment = min(enrichment),
      sd_enrichment = sd(enrichment)
    )
  }
  
  return( metrics )
  
  # davidFunctionalChart <- david$getFunctionalAnnotationChart()
  # davidFunctionalChart <- DAVIDFunctionalAnnotationChart(davidFunctionalChart)
  
  
  # Valores importantes
  #tests$PValue
  #tests$Count
  # y tests$X. que posiblemente es el % que aparece en la página al usar functional annotation chart
  
  # Falta agregar más categorías a este análisis. Probé agregando todas pero crashea
  #david$setAnnotationCategories(david$getAllAnnotationCategoryNames())
  #categories(davidFunChart2)
  
  #return(
  #  annotation_chart <- davidFunctionalChart,
  #  annotation_cluster <- davidCluster
  #)
}

evaluator.metaheuristics <- function(metaheuristic, meta_params, run_evaluator = evaluator.multiobjective.clustering.no.bio, runs = 13, debug = FALSE) {
  
  results <- lapply( 1:runs, function(n) {
    if (debug) {
      message("--------------------------------------------------")
      message(paste("Running metaheuristic on iteration:", n))
      message("\n\n\n")
    }
    
    start.time <- Sys.time()
    
    metaheuristic_results <- do.call(metaheuristic, meta_params)
    iteration_results = run_evaluator(metaheuristic_results, meta_params$dmatrix_expression)
    
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    iteration_results$time <- time.taken
    iteration_results$fitness_counter <- fitness_counter
    
    # Delete data that can't be averaged at the end
    iteration_results$silhouette$silhouette <- NULL
    iteration_results
  })
  
  full_results <- results
  #message( results )
  mean_results <- unlist(results)
  mean_results <- c(by(mean_results, names(mean_results), mean, na.rm = TRUE))
  return( list(
    full_results = full_results,
    mean_results = mean_results
  ))
}

#results <- evaluator.metaheuristics(nsga2.custom, list(dmatrix_expression=dmatrix_expression, dmatrix_biological=dmatrix_biological, population_size = 2, generations = 1, num_clusters = 3, ls_pos=NULL, local_search = NULL, debug=TRUE))
#results <- nsga2.custom(dmatrix_expression, dmatrix_biological, population_size = 80, generations = 10, num_clusters = 5, ls_pos=NULL, local_search = NULL, debug=TRUE)

evaluator.friedman <- function() {
}
