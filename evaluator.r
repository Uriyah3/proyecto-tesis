library(amap)
library(mco)
library(cluster)
library(RDAVIDWebService)
source("nsga2.r")

helper.normalize <- function(data) {
  return( (data - min(data)) / (max(data) - min(data)) )
}

evaluator.multiobjective.clustering <- function( results, dmatrix ) {
  
  num_clusters <- ncol(results$population) - 4
  solution_count <- length(results$clustering)
  
  silhouette_indices = double(solution_count)
  for (c_index in 1:solution_count ) {
    clustering <- results$clustering[[c_index]]
    
    sil <- silhouette( clustering, dmatrix = as.matrix(dmatrix) )
    silhouette_indices[c_index] <- summary(sil)$avg.width
  }
  
  pareto_points <- as.matrix( results$population[ , (num_clusters+1):(num_clusters+2) ] )
  hypervolume <- dominatedHypervolume(pareto_points)
  centered_hypervolume <- dominatedHypervolume(pareto_points, c(1, 1))
  
  pareto_points[, 1] <- helper.normalize(pareto_points[, 1])
  pareto_points[, 2] <- helper.normalize(pareto_points[, 2])
  n_hypervolume <- dominatedHypervolume(pareto_points, c(1, 1))
  
  metrics <- list(
    silhouette = silhouette_indices,
    max_silhouette = max(silhouette_indices),
    mean_silhouette = mean(silhouette_indices),
    min_silhouette = min(silhouette_indices),
    sd_silhouette = sd(silhouette_indices),
    hypervolume = hypervolume, # Raw hypervolume
    normalized_hypervolume = n_hypervolume, # Hypervolume using normalized data to [0,1].
    centered_hypervolume = centered_hypervolume # Hypervolume with (1,1) as reference point
  )
  return( metrics )
}

evaluator.biological.significance <- function(gene_list) {
  #https://david.ncifcrf.gov/webservice/services/DAVIDWebService/authenticate?args0=nicolas.mariangel@usach.cl
  david <- DAVIDWebService$new(email='nicolas.mariangel@usach.cl', url='https://david.ncifcrf.gov/webservice/services/DAVIDWebService')
  
  #david$setAnnotationCategories(david$getAllAnnotationCategoryNames())
  david$addList(gene_list, "ENTREZ_GENE_ID", listName = "Prueba de anotación chart", listType = "Gene")
  tests <- david$getFunctionalAnnotationChart()
  
  # Valores importantes
  #tests$PValue
  #tests$Count
  # y tests$X. que posiblemente es el % que aparece en la página al usar functional annotation chart
  
  davidFunChart2<-DAVIDFunctionalAnnotationChart(tests)
  categories(davidFunChart2)
  
  
}

evaluator.metaheuristics <- function(metaheuristic, meta_params, run_evaluator = evaluator.multiobjective.clustering, runs = 10) {
  results <- lapply( 1:runs, function(n) {
    test <- do.call(metaheuristic, meta_params)
    iteration_results = run_evaluator(test, meta_params$dmatrix_expression)
    iteration_results$silhouette <- NULL
    iteration_results
  })
  #print( results )
  results <- unlist(results)
  results <- c(by(x, names(x), mean, na.rm = TRUE))
}

#resuls <- evaluator.metaheuristics(nsga2.custom, list(dmatrix_expression=dmatrix_expression, dmatrix_biological=dmatrix_biological, population_size = 2, generations = 1, num_clusters = 3, ls_pos=NULL, local_search = NULL, debug=TRUE))
#results <- nsga2.custom(dmatrix_expression, dmatrix_biological, population_size = 80, generations = 10, num_clusters = 5, ls_pos=NULL, local_search = NULL, debug=TRUE)

evaluator.friedman <- function() {
}
