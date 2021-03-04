source("nsga2.r")
source("matrices_evaluation.r")
source("globals.r")
source("evaluator.r")
source("file_utils.r")
library(stringr)
library(ggplot2)

test_file_name <- 'Prostate_13532_31.csv'
test_file <- paste("data/evaluation/", test_file_name, sep="")
data <- read.dataset(test_file)
data <- as.data.frame( t(data) )
colnames(data) <- sapply(colnames(data), function(colname) gsub("[_].*", "", gsub("G", "", colname)))

data <- clean.and.translate.entrez.id(data, "GPL3834")
data <- data[ , order( as.numeric(colnames(data)) ) ]

gene_list <- colnames(data)
dmatrix_expression = expression.matrix(t(data), dataset=test_file_name)
dmatrix_biological = biological.matrix(gene_list, biological_databases$string, dataset=test_file_name)

options(digits.secs = 6)
nsga_results <- vector("list", 11)
pls_results <- vector("list", 11)
mosa_results <- vector("list", 11)
for(i in 1:11) {
  start.time <- Sys.time()
  
  print(paste("Procesando nsga numero", i))
  nsga_results[[i]] <- nsga2.custom(dmatrix_expression, dmatrix_biological, population_size = 50, generations = 20, num_clusters = 4)
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  nsga_results[[i]]$metrics <- evaluator.multiobjective.clustering( nsga_results[[i]], dmatrix_expression )
  nsga_results[[i]]$fitness_counter <- fitness_counter
  nsga_results[[i]]$time <- time.taken
  
  print(nsga_results[[i]])
}

neighborhood <- 0.48
dmatrix_combined <- sqrt(dmatrix_expression**2 + dmatrix_biological**2)
# Find genes that are close to one another
neighborhood_matrix <- sapply(gene_list, function(gene) {
  neighborhood_genes <- dmatrix_combined[gene, , drop=FALSE]
  apply(neighborhood_genes, 1, function(x) colnames(neighborhood_genes)[which(x > 0.00 & x < neighborhood)] )
})
for(i in 1:11) {
  start.time <- Sys.time()
  
  print(paste("Procesando pls numero", i))
  initial_population <- generate.initial.population(gene_list, 80, 7)
  fitness_hash <<- hash()
  fitness_counter <<- 0
  
  pls_results[[i]]$raw <- local.search.pareto.local.search(population_size = 80, population = initial_population, num_clusters = 7, gene_list = gene_list, dmatrix_expression = dmatrix_expression, dmatrix_biological = dmatrix_biological, neighborhood_matrix = neighborhood_matrix, ordering_fn = operator.nsga2.sorting.and.crowding, fitness_fn = fitness.medoid.wg, acceptance_criteria_fn = helper.non.dominated, rank_cutoff = 2, max_generations = 150)
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  pls_results[[i]]$results <- generate.results(nrow(pls_results[[i]]$raw), 7, pls_results[[i]]$raw, dmatrix_expression, dmatrix_biological)
  pls_results[[i]]$metrics <- evaluator.multiobjective.clustering( pls_results[[i]]$results, dmatrix_expression )
  pls_results[[i]]$fitness_counter <- fitness_counter
  pls_results[[i]]$time <- time.taken
  
  clear(fitness_hash)
  print(pls_results[[i]]$metrics)
  print(pls_results[[i]]$time)
}

for(i in 1:11) {
  start.time <- Sys.time()
  
  print(paste("Procesando mosa numero", i))
  initial_population <- generate.initial.population(gene_list, 60, 7)
  fitness_hash <<- hash()
  fitness_counter <<- 0
  
  mosa_results[[i]]$raw <- local.search.mosa(population_size = 60, population = initial_population, num_clusters = 7, gene_list = gene_list, dmatrix_expression = dmatrix_expression, dmatrix_biological = dmatrix_biological, neighborhood_matrix = neighborhood_matrix, ordering_fn = operator.nsga2.sorting.and.crowding, fitness_fn = fitness.medoid.wg, max_steps=850)
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  mosa_results[[i]]$results <- generate.results(nrow(mosa_results[[i]]$raw), 7, mosa_results[[i]]$raw, dmatrix_expression, dmatrix_biological)
  mosa_results[[i]]$metrics <- evaluator.multiobjective.clustering( mosa_results[[i]]$results, dmatrix_expression )
  mosa_results[[i]]$fitness_counter <- fitness_counter
  mosa_results[[i]]$time <- time.taken
  
  clear(fitness_hash)
  print(mosa_results[[i]]$metrics)
  print(mosa_results[[i]]$time)
}

process.metric.data <- function(data) {
  max_silhouette <- vector("numeric", 11)
  mean_silhouette <- vector("numeric", 11)
  min_silhouette <- vector("numeric", 11)
  c_hypervolume <- vector("numeric", 11)
  n_hypervolume <- vector("numeric", 11)
  hypervolume <- vector("numeric", 11)
  data_fitness_counter <- 0
  mean_time <- 0
  for(i in 1:11) {
    metrics <- data[[i]]$metrics
    
    max_silhouette[[i]] <- max(metrics$silhouette)
    mean_silhouette[[i]] <- mean(metrics$silhouette)
    min_silhouette[[i]] <- min(metrics$silhouette)
    c_hypervolume[[i]] <- metrics$centered_hypervolume
    n_hypervolume[[i]] <- metrics$normalized_hypervolume
    hypervolume[[i]] <- metrics$hypervolume
    
    data[[i]]$avg_silhouette <- mean(metrics$silhouette)
    data[[i]]$sd_silhouette <- sd(metrics$silhouette)
    
    data_fitness_counter <- data_fitness_counter + data[[i]]$fitness_counter
    mean_time <- mean_time + data[[i]]$time
  }
  data_fitness_counter <- data_fitness_counter / 11
  mean_time <- mean_time / 11
  if(mean_time < 100) mean_time <- mean_time * 60
  
  silhouette_data <- as.data.frame(cbind(min_silhouette, mean_silhouette, max_silhouette))
  colnames(silhouette_data) <- c('Silueta mínima', 'Silueta promedio', 'Silueta máxima')
  silhouette_data <- melt(silhouette_data)
  
  data$silhouette_data = silhouette_data
  data$c_hypervolume = c_hypervolume
  data$n_hypervolume = n_hypervolume
  data$hypervolume = hypervolume
  data$mean_time = mean_time
  data$mean_fitness_counter= data_fitness_counter
  
  return(data)
}

show_results <- function(data, type="NSGA-II") {
  data <- process.metric.data(data)
  silhouette_data <- data$silhouette_data
  c_hypervolume <- data$c_hypervolume
  
  my.plot(
    silhouette_data,
    str_interp("Resultados índice silueta para 11 ejecuciones de ${type}"),
    "Tipo",
    "Índice Silueta",
    str_interp("${type}_resultados_silueta.png")
  )
  
  my.plot(
    melt(as.data.frame(c_hypervolume)),
    str_interp("Hypervolumen usando (1,1) como punto de referencia (${type})"),
    "",
    "Hypervolumen normalizado",
    str_interp("${type}_resultados_hypervolumen.png")
  )
  
  print(str_interp("Tiempo promedio de ${type}: ${format(.POSIXct(data$mean_time,tz='GMT'), '%H:%M:%S')}"))
  print(str_interp("Calculos de fitness en promedio: ${data$mean_fitness_counter}"))
  
  return(data)
}

my.plot <- function(data, title, xlab, ylab, filename) {
  ggplot(data, aes(x=variable, y=value)) + 
    geom_boxplot(
      
      # custom boxes
      color="blue",
      fill="blue",
      alpha=0.25,
      
      # custom outliers
      outlier.colour="red",
      outlier.fill="red",
      outlier.size=4
    ) +
    ggtitle(title) +
    xlab(xlab) +
    ylab(ylab)
  
  ggsave(str_interp(filename), device="png", path="plots")
}

compare_results <- function(data1, data2, data3, data4) {
  for(i in 1:11) {data1[[i]]$metrics <- evaluator.multiobjective.clustering( data1[[i]], dmatrix_expression )}
  for(i in 1:11) {data2[[i]]$metrics <- evaluator.multiobjective.clustering( data2[[i]], dmatrix_expression )}
  for(i in 1:11) {data3[[i]]$metrics <- evaluator.multiobjective.clustering( data3[[i]]$results, dmatrix_expression )}
  for(i in 1:11) {data4[[i]]$metrics <- evaluator.multiobjective.clustering( data4[[i]]$results, dmatrix_expression )}
  data1 <- process.metric.data(data1)
  data2 <- process.metric.data(data2)
  data3 <- process.metric.data(data3)
  data4 <- process.metric.data(data4)
  
  print(str_interp("Tiempo promedio de NSGA-II GO (p=80, g=40, k=7): ${format(.POSIXct(data1$mean_time,tz='GMT'), '%H:%M:%S')}"))
  print(str_interp("Calculos de fitness en promedio: ${data1$mean_fitness_counter}"))
  
  print(str_interp("Tiempo promedio de NSGA-II STRING (p=50, g=20, k=4): ${format(.POSIXct(data2$mean_time,tz='GMT'), '%H:%M:%S')}"))
  print(str_interp("Calculos de fitness en promedio: ${data2$mean_fitness_counter}"))
  
  print(str_interp("Tiempo promedio de PLS GO (p=80, k=7, iter=150): ${format(.POSIXct(data3$mean_time,tz='GMT'), '%H:%M:%S')}"))
  print(str_interp("Calculos de fitness en promedio: ${data3$mean_fitness_counter}"))
  
  print(str_interp("Tiempo promedio de MOSA GO (p=80, k=7, iter=850): ${format(.POSIXct(data4$mean_time,tz='GMT'), '%H:%M:%S')}"))
  print(str_interp("Calculos de fitness en promedio: ${data4$mean_fitness_counter}"))
  
  c_hypervolume <- cbind(data1$c_hypervolume, data2$c_hypervolume, data3$c_hypervolume, data4$c_hypervolume)
  n_hypervolume <- cbind(data1$n_hypervolume, data2$n_hypervolume, data3$n_hypervolume, data4$n_hypervolume)
  hypervolume <- cbind(data1$hypervolume, data2$hypervolume, data3$hypervolume, data4$hypervolume)
  colnames(c_hypervolume) <- c("NSGA-II", "NSGA-II STRING", "PLS", "MOSA")
  colnames(n_hypervolume) <- c("NSGA-II", "NSGA-II STRING", "PLS", "MOSA")
  colnames(hypervolume) <- c("NSGA-II", "NSGA-II STRING", "PLS", "MOSA")
  
  
  levels(data1$silhouette_data$variable)[match("Silueta máxima",levels(data1$silhouette_data$variable))] <- "NSGA-II Max"
  levels(data1$silhouette_data$variable)[match("Silueta promedio",levels(data1$silhouette_data$variable))] <- "NSGA-II Avg"
  levels(data1$silhouette_data$variable)[match("Silueta mínima",levels(data1$silhouette_data$variable))] <- "NSGA-II Min"
  levels(data2$silhouette_data$variable)[match("Silueta máxima",levels(data2$silhouette_data$variable))] <- "NSGA-II STRING Max"
  levels(data2$silhouette_data$variable)[match("Silueta promedio",levels(data2$silhouette_data$variable))] <- "NSGA-II STRING Avg"
  levels(data2$silhouette_data$variable)[match("Silueta mínima",levels(data2$silhouette_data$variable))] <- "NSGA-II STRING Min"
  levels(data3$silhouette_data$variable)[match("Silueta máxima",levels(data3$silhouette_data$variable))] <- "PLS Max"
  levels(data3$silhouette_data$variable)[match("Silueta promedio",levels(data3$silhouette_data$variable))] <- "PLS Avg"
  levels(data3$silhouette_data$variable)[match("Silueta mínima",levels(data3$silhouette_data$variable))] <- "Pls Min"
  levels(data4$silhouette_data$variable)[match("Silueta máxima",levels(data4$silhouette_data$variable))] <- "MOSA Max"
  levels(data4$silhouette_data$variable)[match("Silueta promedio",levels(data4$silhouette_data$variable))] <- "MOSA Avg"
  levels(data4$silhouette_data$variable)[match("Silueta mínima",levels(data4$silhouette_data$variable))] <- "MOSA Min"
  silhouette_data <- as.data.frame(rbind(data1$silhouette_data, data2$silhouette_data, data3$silhouette_data, data4$silhouette_data))
  
  my.plot(silhouette_data[grep("Max", silhouette_data$variable), ], "Comparacion entre la silueta máxima de 4 algoritmos", "Algoritmo", "Silueta", "4way-comparacion-silueta_max.png")
  my.plot(silhouette_data[grep("Avg", silhouette_data$variable), ], "Comparacion entre la silueta promedio de 4 algoritmos", "Algoritmo", "Silueta", "4way-comparacion-silueta_avg.png")
  my.plot(silhouette_data[grep("Min", silhouette_data$variable), ], "Comparacion entre la silueta mínima de 4 algoritmos", "Algoritmo", "Silueta", "4way-comparacion-silueta_min.png")
  my.plot(melt(as.data.frame(c_hypervolume)), "Comparacion del hypervolumen centrado en (1,1) de 4 algoritmos", "Algoritmo", "Hypervolumen", "4way-comparacion-c_hypervolume.png")
  my.plot(melt(as.data.frame(n_hypervolume)), "Comparacion del hypervolumen normalizado de 4 algoritmos", "Algoritmo", "Hypervolumen", "4way-comparacion-n_hypervolume.png")
  my.plot(melt(as.data.frame(hypervolume)), "Comparacion del hypervolumen de 4 algoritmos", "Algoritmo", "Hypervolumen", "4way-comparacion-hypervolume.png")
}

nsga_results <- show_results(nsga_results, "NSGA-II-STRING")
pls_results <- show_results(pls_results, "PLS")
mosa_results <- show_results(mosa_results, "MOSA")

compare_results(readRDS('cache/nsga_results_p.rda'), readRDS('cache/nsga_string_results.rda'), readRDS('cache/pls-results.rda'), readRDS('cache/mosa_results.rda'))


#plot(nsga_results[[1]]$population[, colnames(nsga_results[[1]]$population) %in% c('objective_exp', 'objective_bio')])
# d <- nsga_results[[1]]$population[, colnames(nsga_results[[1]]$population) %in% c('objective_exp', 'objective_bio')]
# ggplot() +
#     geom_step(data=d, mapping=aes(x=objective_exp, y=objective_bio)) +
#     geom_step(data=d, mapping=aes(x=objective_exp, y=objective_bio), direction="vh", linetype=3) +
#     geom_point(data=d, mapping=aes(x=objective_exp, y=objective_bio), color="red") + xlim(0, 1) + ylim(0, 1)