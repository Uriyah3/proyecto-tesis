library(amap)
library(mco)
library(cluster)
library(stringr)
library(RDAVIDWebService)
library(stats)
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

evaluator.multiobjective.clustering.no.bio <- function( results, dmatrix, debug = FALSE, dataset_name=NULL, bio=NULL, iter=NULL, which.x=NULL ) {
  
  silhouette_results <- evaluator.silhouette( results$clustering, dmatrix, debug=debug, dataset_name=dataset_name, bio=bio, iter=iter )
  hypervolume_results <- evaluator.hypervolume( results$population, debug=debug, dataset_name=dataset_name, bio=bio, iter=iter )
  
  metrics <- list(
    silhouette = silhouette_results,
    hypervolume = hypervolume_results
  )
  
  return( metrics )
}

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

evaluator.silhouette <- function( clustering, dmatrix, debug = FALSE, dataset_name=NULL, bio=NULL, iter=NULL ) {
  if (!is.null(metrics <- load.evaluation.from.cache(dataset_name, bio, iter, 'silhouette'))) {
    return(metrics)
  }
  
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
    if(length(gene_list) <= 3) {
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
      while( is.na(results[[cluster]]) && attempt <= 5 ) {
        attempt <- attempt + 1
        try(
          results[[cluster]] <- evaluator.biological.anotate.list( gene_list, debug )
        )
      }
    }
    
    if (!is.list(results[[cluster]]$enrichment) && is.na(results[[cluster]]$enrichment)) {
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
    enrichment <- unlist(lapply(results, function(result) {
      return(unlist(result$enrichment))
    }))
    results$cluster_count <- sum( unlist(lapply(results, '[[', 'cluster_count')) )
    results$max_enrichment <- max(enrichment)
    results$mean_enrichment <- mean(enrichment)
    results$min_enrichment <- min(enrichment)
    results$sd_enrichment <- sd(enrichment)
    
    store.evaluation.to.cache(results, dataset_name, bio, iter, 'biological')
  }
  
  return( results )
}

#' Used to "bypass" the daily job limit in DAVID's web service
jobs <<- 190
email.list <<- NULL
email.rotator <- function() {
  if(is.null(email.list)) {
    email.list <<- readLines('mail_list.txt')
    email.start <<- Sys.time()
  }
  
  jobs <<- jobs - 1
  
  if(jobs <= 0) {
    email.list <<- email.list[-1]
    jobs <<- 190
    
    wait.time <- sample(100:900, 1)
    message(str_interp("Waiting ${wait.time} seconds to continue using next email"))
    Sys.sleep(wait.time)
  }
  
  # No quedan elementos al eliminar
  if (length(email.list) == 0) {
    message("Waiting until a day has passed to reuse the email list...")
    email.end <- Sys.time()
    time.to.day <- max(0, round(96000 - difftime(email.end,email.start,units="secs")))
    message(str_interp("Wating ${time.to.day} seconds"))
    Sys.sleep(time.to.day)
    email.list <<- NULL
    return(email.rotator())
  }
  
  return(email.list[1])
}

evaluator.biological.anotate.list <- function( gene_list, debug = FALSE ) {
  current.email = email.rotator()
  if(debug) {
    message(str_interp("Current mail usage: ${current.email} ${jobs}/190"))
  }
  #https://david.ncifcrf.gov/webservice/services/DAVIDWebService/authenticate?args0=nicolas.mariangel@usach.cl
  david <- DAVIDWebService$new(email=current.email, url='https://david.ncifcrf.gov/webservice/services/DAVIDWebService')
  # Se cae a veces con el siguiente error: 
  # [INFO] Unable to sendViaPost to url[https://david.ncifcrf.gov/webservice/services/DAVIDWebService]
  # java.net.SocketTimeoutException: Read timed out
  # Considerar que también tira este error cuando la lista de genes es > 3000
  setTimeOut(david, 1500000)
  
  if (length(gene_list) >= 3000) {
    return(
      list(
        cluster_count = min(length(gene_list), 100),
        enrichment = list(0)
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
      message(paste("DAVID falló en encontrar enrichment o se cayó la biblioteca para la lista:", paste0(gene_list, collapse=', ')))
    }
    metrics <- list(
      cluster_count = min(100, length(gene_list)),
      enrichment = NA
    )
  } else {
    metrics <- list(
      cluster_count = length(clusters),
      enrichment = enrichment
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

#results <- evaluator.metaheuristics(nsga2.custom, list(dmatrix_expression=dmatrix_expression, dmatrix_biological=dmatrix_biological, population_size = 2, generations = 1, num_clusters = 3, ls_pos=NULL, local_search = NULL, debug=TRUE))
#results <- nsga2.custom(dmatrix_expression, dmatrix_biological, population_size = 80, generations = 10, num_clusters = 5, ls_pos=NULL, local_search = NULL, debug=TRUE)

