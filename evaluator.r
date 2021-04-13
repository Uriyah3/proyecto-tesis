library(amap)
library(mco)
library(cluster)
library(stringr)
library(RDAVIDWebService)
source("nsga2.r")

helper.normalize <- function(data) {
  return( (data - min(data)) / (max(data) - min(data)) )
}

# https://stackoverflow.com/questions/10256503/function-for-median-similar-to-which-max-and-which-min-extracting-median-r
which.median = function(x) {
  if (length(x) %% 2 != 0) {
    which(x == median(x))
  } else if (length(x) %% 2 == 0) {
    a = sort(x)[c(length(x)/2, length(x)/2+1)]
    c(which(x == a[1]), which(x == a[2]))
  }
}

evaluator.multiobjective.clustering <- function( results, dmatrix, debug = FALSE, which.x = which.max, dataset_name=NULL ) {
  
  silhouette_results <- evaluator.silhouette( results$clustering, dmatrix, debug)
  hypervolume_results <- evaluator.hypervolume( results$population, debug )
  
  cond_sil <- which.x(silhouette_results$silhouette)[1]
  if (debug) {
    message(paste("Se analiza biologicamente solucion con silueta =", silhouette_results$silhouette[[cond_sil]]))
  }
  biology_results <- evaluator.biological.significance( results$clustering[[cond_sil]], colnames(dmatrix), debug )
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

evaluator.multiobjective.clustering.custom.bio <- function( results, dmatrix_expression, dataset_name, debug=FALSE, which.x = which.max ) {
  bio_sources = list(
    go = "gene ontology",
    string = "STRING",
    kegg = "KEGG-pathway",
    disgenet_dis = "disgenet-disease"
  )
  
  metrics <- evaluator.multiobjective.clustering.no.bio(results, dmatrix_expression, debug)
  
  gene_list <- colnames(dmatrix_expression)
  biology_silhouette <- list()
  for(source in bio_sources) {
    dmatrix_bio <- biological.matrix(gene_list, source, dataset=dataset_name)
    biology_silhouette[[source]] <- evaluator.silhouette(results$clustering, dmatrix_bio, debug)
  }
  dmatrix_bio <- NULL
  
  metrics$biological <- biology_silhouette
  return(metrics)
}

evaluator.multiobjective.clustering.no.bio <- function( results, dmatrix, debug = FALSE, dataset_name=NULL ) {
  silhouette_results <- evaluator.silhouette( results$clustering, dmatrix, debug )
  hypervolume_results <- evaluator.hypervolume( results$population, debug )
  
  metrics <- list(
    silhouette = silhouette_results,
    hypervolume = hypervolume_results
  )
  return( metrics )
}

evaluator.hypervolume <- function( population, debug = FALSE ) {
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

evaluator.silhouette <- function( clustering, dmatrix, debug = FALSE ) {
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

evaluator.biological.significance <- function( clustering, full_gene_list, debug = FALSE ) {
  # Agregar los nombres de los genes a la cabecera
  clustering <- as.data.frame(t(clustering))
  colnames(clustering) <- full_gene_list
  
  num_clusters <- max(clustering)
  
  results <- list()
  for (cluster in 1:num_clusters) {
    gene_list <- colnames(clustering[ , clustering == cluster, drop=FALSE ])
    
    if (debug) {
      message(paste("Procesando lista#", cluster, " con ", length(gene_list), " genes...", sep=""))
    }
    
    results[[cluster]] <- evaluator.biological.anotate.list( gene_list, debug )
  }
  
  return( results )
}

evaluator.biological.anotate.list <- function( gene_list, debug = FALSE ) {
  #https://david.ncifcrf.gov/webservice/services/DAVIDWebService/authenticate?args0=nicolas.mariangel@usach.cl
  david <- DAVIDWebService$new(email='nicolas.mariangel@usach.cl', url='https://david.ncifcrf.gov/webservice/services/DAVIDWebService')
  # Se cae a veces con el siguiente error: 
  # [INFO] Unable to sendViaPost to url[https://david.ncifcrf.gov/webservice/services/DAVIDWebService]
  # java.net.SocketTimeoutException: Read timed out
  # Considerar que también tira este error cuando la lista de genes es > 3000
  setTimeOut(david, 300000)
  
  if (length(gene_list) >= 3000 || length(gene_list) <= 2) {
    if (debug) {
      message(paste("No se puede utilizar DAVID con esta lista de genes, dado que son", length(gene_list), "genes"))
      if (length(gene_list) <= 2) {
        message(paste0(gene_list, collapse=","))
      }
    }
    return(
      list(
        cluster_count = min(length(gene_list), 100),
        max_enrichment = 0,
        mean_enrichment = 0,
        min_enrichment = 0,
        sd_enrichment = 0
      )
    )
  }
  
  
  david$addList(gene_list, "ENTREZ_GENE_ID", listName = paste("Prueba de anotacion chart", sample(1:10000, 1)), listType = "Gene")
  david$setAnnotationCategories(c("ENTREZ_GENE_ID", "BIOCARTA", "BBID", "BIOGRID_INTERACTION", "CGAP_EST_QUARTILE", "CGAP_SAGE_QUARTILE", "CHROMOSOME", "ENSEMBL_GENE_ID", "ENTREZ_GENE_SUMMARY", "GAD_DISEASE", "GAD_DISEASE_CLASS", "GENERIF_SUMMARY", "GNF_U133A_QUARTILE", "GOTERM_BP_ALL", "GOTERM_BP_DIRECT", "GOTERM_CC_ALL", "GOTERM_CC_DIRECT",  "GOTERM_MF_ALL", "GOTERM_MF_DIRECT", "HIV_INTERACTION", "HIV_INTERACTION_CATEGORY", "HIV_INTERACTION_PUBMED_ID", "KEGG_PATHWAY", "MINT", "OMIM_DISEASE", "PFAM", "PIR_SEQ_FEATURE", "PIR_SUMMARY", "PIR_SUPERFAMILY", "PRINTS", "PRODOM", "PROSITE", "PUBMED_ID", "REACTOME_PATHWAY", "SMART", "SP_COMMENT", "SP_COMMENT_TYPE", "SUPFAM", "TIGRFAMS", "UCSC_TFBS", "UNIGENE_EST_QUARTILE", "UP_KEYWORDS", "UP_SEQ_FEATURE", "UP_TISSUE"))
  davidCluster <- david$getClusterReport()
  
  clusters <- davidCluster@cluster
  enrichment <- sapply(clusters, `[[`, 'EnrichmentScore')
  
  if (length(enrichment) == 0) {
    if (debug) {
      message("DAVID falló en encontrar enrichment o se cayó la biblioteca")
    }
    metrics <- list(
      cluster_count = min(100, length(gene_list)),
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

evaluator.metaheuristics <- function(metaheuristic, meta_params, run_evaluator = evaluator.multiobjective.clustering.no.bio, runs = 13, debug = FALSE, dataset_name = NULL) {
  
  results <- lapply( 1:runs, function(n) {
    if (debug) {
      message("--------------------------------------------------")
      message(paste("Running metaheuristic on iteration:", n))
      message("\n\n\n")
    }
    
    start.time <- Sys.time()
    
    metaheuristic_results <- do.call(metaheuristic, meta_params)
    iteration_results = run_evaluator(metaheuristic_results, meta_params$dmatrix_expression, debug=debug, dataset_name=dataset_name)
    
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

reconstruct.metaheuristic.saved.results <- function(dataset.name, identifier, run_evaluator = evaluator.multiobjective.clustering.no.bio, runs = 13, debug = FALSE) {
  dmatrix_expression <- expression.matrix(NULL, dataset=dataset.name)
  results <- lapply( 1:runs, function(n) {
    if (debug) {
      message("--------------------------------------------------")
      message(paste("Loading metaheuristic results on iteration:", n))
      message("\n")
    }
    
    filename <- build.saved.results.filename(dataset.name, identifier, n)
    iteration_results <- readRDS(filename)
    results = run_evaluator(iteration_results$nsga, dmatrix_expression, debug=debug, dataset_name=dataset.name)
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

#results <- evaluator.metaheuristics(nsga2.custom, list(dmatrix_expression=dmatrix_expression, dmatrix_biological=dmatrix_biological, population_size = 2, generations = 1, num_clusters = 3, ls_pos=NULL, local_search = NULL, debug=TRUE))
#results <- nsga2.custom(dmatrix_expression, dmatrix_biological, population_size = 80, generations = 10, num_clusters = 5, ls_pos=NULL, local_search = NULL, debug=TRUE)



evaluator.friedman <- function() {
}
