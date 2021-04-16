source("package_installer.r")
source("nsga2.r")
source("matrices_evaluation.r")
source("globals.r")
source("evaluator.r")
source("file_utils.r")
library(future)
library(profvis)
options(bitmapType='cairo')

datasets <- list(
  GSE89116 = list(
    name = "Breast_GSE89116",
    chip = "GPL6947",
    type = "evaluation",
    cancer = "Breast"
  ),
  GSE53757 = list(
    name = "Renal_GSE53757",
    chip = "GPL570",
    type = "evaluation",
    cancer = "Renal"
  ),
  GSE31189 = list(
    name = "Bladder_GSE31189",
    chip = "GPL570",
    type = "evaluation",
    cancer = "Bladder"
  ),
  GSE50161 = list(
    name = "Brain_GSE50161",
    chip = "GPL570",
    type = "evaluation",
    cancer = "Brain"
  ),
  GSE6919_U95Av2 = list(
    name = "Prostate_GSE6919_U95Av2",
    chip = "GPL8300",
    type = "evaluation",
    cancer = "Prostate"
  ),
  GSE6919_U95B = list(
    name = "Prostate_GSE6919_U95B",
    chip = "GPL92",
    type = "training",
    cancer = "Prostate"
  ),
  GSE6919_U95C = list(
    name = "Prostate_GSE6919_U95C",
    chip = "GPL93",
    type = "training",
    cancer = "Prostate"
  ),
  GSE28497 = list(
    name = "Leukemia_GSE28497",
    chip = "GPL96",
    type = "training",
    cancer = "Leukemia"
  ),
  # GSE45827 = list(
  #   name = "Breast_GSE45827",
  #   chip = "GPL570",
  #   type = "evaluation",
  #   cancer = "Breast"
  # ),
  # GSE19804 = list(
  #   name = "Lung_GSE19804",
  #   chip = "GPL570",
  #   type = "evaluation",
  #   cancer = "Lung"
  # ),
  GSE22405 = list(
    name = "Liver_GSE22405",
    chip = "GPL96",
    type = "training",
    cancer = "Liver"
  ),
  GSE60502 = list(
    name = "Liver_GSE60502",
    chip = "GPL96",
    type = "training",
    cancer = "Liver"
  ),
  GSE6008 = list(
    name = "Ovary_GSE6008",
    chip = "GPL96",
    type = "training",
    cancer = "Ovary"
  ),
  # GSE15824 = list(
  #   name = "Brain_GSE15824",
  #   chip = "GPL570",
  #   type = "training",
  #   cancer = "Brain"
  # ),
  GSE44861 = list(
    name = "Colorectal_GSE44861",
    chip = "GPL3921",
    type = "training",
    cancer = "Colorectal"
  ),
  GSE77953 = list(
    name = "Colorectal_GSE77953",
    chip = "GPL96",
    type = "training",
    cancer = "Colorectal"
  ),
  GSE10797 = list(
    name = "Breast_GSE10797",
    chip = "GPL571",
    type = "training",
    cancer = "Breast"
  ),
  # GSE7904 = list(
  #   name = "Breast_GSE7904",
  #   chip = "GPL570",
  #   type = "training",
  #   cancer = "Breast"
  # ),
  # GSE16515 = list(
  #   name = "Pancreatic_GSE16515",
  #   chip = "GPL570",
  #   type = "training",
  #   cancer = "Pancreatic"
  # ),
  GSE9476 = list(
    name = "Leukemia_GSE9476",
    chip = "GPL96",
    type = "training",
    cancer = "Leukemia"
  )
)


run.default <- function() {
  test_file_name <- 'Renal_GSE53757'
  test_file <- paste("data/evaluation/", test_file_name, '.csv.gz', sep="")
  data <- read.dataset(test_file)
  
  data <- clean.and.translate.entrez.id(data, "GPL570")
  data <- data[ , order( as.numeric(colnames(data)) ) ]
  
  gene_list <- colnames(data)
  dmatrix_expression = expression.matrix(t(data), dataset=test_file_name)
  dmatrix_biological = biological.matrix(gene_list, biological_databases$go, dataset=test_file_name)
  
  #set.seed(1048)
  #Rprof("profile1.out", line.profiling=TRUE, memory.profiling=TRUE)
  results <- nsga2.custom(dmatrix_expression, dmatrix_biological, population_size = 80, generations = 40, num_clusters = 5, ls_pos=NULL, local_search = NULL)
  #Rprof(NULL)
  
  #summaryRprof("profile1.out", lines = "show")
}

#' Used along run.r.sample.test to check if there is any bias in R's sample method.
#' 
#' @param data Array of data to be tested for random selections
#' @param times Integer. How many times is sample() run.
#' @param test_n Integer. Test n° identifier. Used to save png to different file.
#' @param test_type String. What kind of test. Used to save png to different file.
#' 
r.sample.test <- function(data, times, test_n, test_type) {
  appearances <- hash()
  for (i in data) { appearances[[paste(i)]] <- 0 }
  
  for (i in 1:times) {
    n <- sample(data, 1)
    appearances[[paste(n)]] <- appearances[[paste(n)]] + 1
  }
  
  data.string <- str_interp("${min(data)}:${max(data)}")
  png(filename=str_interp("plots/r_${test_type}_sample_test_${min(data)}_to_${max(data)}_executed_${times}_times_${test_n}.png"))
  plot(as.list(appearances), data, xlab=str_interp("N° apariciones en sample(${data.string}, 1)"), ylab=data.string, main=str_interp("Prueba de distribución de sample de R (${times} veces)."))
  dev.off()
}

#' Run multiple sample tests and save results in png files 
#' under the 'plots' directory
#' 
run.r.sample.test <- function() {
  
  set.seed(12521)
  
  for(big_test in 1:3) {
    r.sample.test(1:10000, 10000000, big_test, "big")
  }
  
  for(small_test in 1:10) {
    r.sample.test(1:1000, 1000000, small_test, "small")
  }
}

#' Load a dataset and translate its genes (colnames) to ENTREZ_GENE_ID
#' 
#' @param dataset A list with the attributes: name (file_name without extension),
#'   chip (Name of the platform used, example: "GPL_570"), type ("evaluation" or
#'   "training") and cancer
load.dataset <- function(dataset) {
  
  file_name <- paste(dataset$name, '.csv.gz', sep='')
  file_path <- paste('data/', dataset$type, '/', file_name, sep="")
  data <- read.dataset(file_path)
  
  data <- clean.and.translate.entrez.id(data, dataset$chip)
  data <- data[ , order( as.numeric(colnames(data)) ) ]
}

#' Calculate all biological distance matrices and store them in .rda files under
#' the 'cache' directory to save time in the future.
#' 
#' @param workers Integer. Number of cpu cores / processes to run in parallel.
#' It's not recommended to use more than 10 since each process can eat around 
#' 4-8% of memory (~8 GB RAM)
#' 
precalculate.biological.dmatrix <- function(workers = 7, expression_nbproc = 8) {
  plan(multisession, gc = TRUE, workers = workers)
  
  important_biological_databases <- biological_databases
  important_biological_databases$mesh <- NULL
  important_biological_databases$disgenet_pw <- NULL
  important_biological_databases$disease <- NULL
  for (dataset in datasets) {
    message(dataset)
    data <- load.dataset(dataset)
    gene_list <- colnames(data)
    for (biological_source in important_biological_databases) {
      future({biological.matrix(gene_list, biological_source, dataset=dataset$name)}, gc = TRUE)
    }
    expression.matrix(t(data), dataset=dataset$name, nbproc = expression_nbproc)
  }
}

#' For each dataset, check how many genes are there before cleaning, after
#' eliminating genes with no ENTREZ_GENE_ID translation and after translating
#' each gene.
#' 
check.gene.translator.effectiveness <- function() {
  
  translation.results <- list()
  for (dataset in datasets) {
    dataset.results <- list()
    
    # Code copied from the load.dataset method
    file_name <- paste(dataset$name, '.csv.gz', sep='')
    file_path <- paste('data/', dataset$type, '/', file_name, sep="")
    data <- read.dataset(file_path)
    
    dataset.results[["all_genes"]] <- ncol(data)
    
    chip = dataset$chip
    # Code copied from the clean.and.translate.entrez.id method
    gene_list = colnames(data)
    if( substring(gene_list[1], 1, 1) == "X" ) {
      gene_list = sapply(gene_list, function(gene) return(substring(gene, 2)), USE.NAMES = FALSE)
      
      colnames(data) <- gene_list
    }
    
    gene_translator <- process.gpl(chip)
    genes_to_keep <- gene_translator[gene_translator$ENTREZ_GENE_ID != "", , drop=FALSE]
    
    # Translate to ENTREZ_GENE_ID
    data <- data[, colnames(data) %in% rownames(genes_to_keep), drop=FALSE]
    
    dataset.results[["genes_with_entrez_id"]] <- ncol(data)
    
    colnames(data) <- sapply(colnames(data), function(gene_name) {
      unlist(strsplit(as.character(gene_translator[gene_name, ]), " /// "))[1]
    })
    # remove duplicates
    data <- data[, !duplicated(colnames(data))]
    
    dataset.results[["non_duplicated_genes"]] <- ncol(data)
    
    translation.results[[dataset$name]] <- dataset.results
  }
  
  saveRDS(translation.results, 'cache/translation_ineffectiveness.rda')
  translation.results
}

#' Compare how similar are each GPL's group of entrez_gene_id in respect
#' to each other.
#'
compare.gpl.platforms <- function() {
  chips <- unique( unlist( lapply(datasets, '[[', 'chip') ) )
  
  chips_data <- list()
  chips_comparison <- as.data.frame(matrix(NA, length(chips), length(chips)))
  colnames(chips_comparison) <- chips
  rownames(chips_comparison) <- chips
  
  for (chip in chips) {
    gene_to_entrez <- process.gpl(chip)
    gene_to_entrez <- gene_to_entrez[gene_to_entrez$ENTREZ_GENE_ID != "", , drop=FALSE]
    for (i in 1:nrow(gene_to_entrez)) {
      gene_to_entrez[i, ] <- unlist(strsplit(as.character(gene_to_entrez[i, ]), " /// "))[1]
    }
    
    chips_data[[chip]] <- gene_to_entrez
  }
  
  for (first_chip in chips) {
    for (second_chip in chips) {
      if (first_chip == second_chip) { next }
      
      first_chip_data <- chips_data[[first_chip]]
      second_chip_data <- chips_data[[second_chip]]
      
      chips_comparison[first_chip, second_chip] <- round(length(intersect(first_chip_data$ENTREZ_GENE_ID, second_chip_data$ENTREZ_GENE_ID)) / length(union(first_chip_data$ENTREZ_GENE_ID, second_chip_data$ENTREZ_GENE_ID)), 3)
    }
  }
  
  saveRDS(chips_comparison, 'cache/chips_entrez_gene_id_jaccard.rda')
  chips_comparison
}

bio.test <- function(debug = TRUE)
{
  dataset <- datasets$GSE6919_U95C
  dmatrix_expression = expression.matrix(t(data), dataset=dataset$name)
  dmatrix_biological = biological.matrix(gene_list, biological_databases[['kegg']], dataset=dataset$name)
  
  results_processed <- results_processed <- nsga2.custom(dmatrix_expression, dmatrix_biological, population_size = 40, num_clusters = 12, crossover_ratio = 0.7950, mutation_ratio = 0.0714, tour_size = 6, local_search = local_search_algorithms[['lmols']], neighborhood = 0.3214, ls_budget = 44.5125, ls_pos = 2, debug = TRUE, evaluations=8000)
  results_random <- nsga2.custom(dmatrix_expression, dmatrix_biological, evaluations=100, population_size = 40, num_clusters = 12, crossover_ratio = 0.7950, mutation_ratio = 0.0714, tour_size = 6, debug=TRUE)
  
  metrics_processed <- evaluator.multiobjective.clustering.custom.bio(results_processed, dmatrix_expression, dataset$name, which.x=which.median, debug=debug)
  metrics_random <- evaluator.multiobjective.clustering.custom.bio(results_random, dmatrix_expression, dataset$name, which.x=which.median, debug=debug)
  
  return(list(
    res_random = results_random,
    res_processed = results_processed,
    random = metrics_random,
    processed = metrics_processed
  ))
}

best_params <- list(
  go = list(
    evaluations=10000,
    population_size=83,
    num_clusters=10,
    crossover_ratio=0.5710,
    mutation_ratio=0.1205,
    tour_size=4,
    local_search=local_search_algorithms$lmols,
    neighborhood=0.6756,
    ls_budget=51.21,
    ls_pos=2
  ),
  string = list(
    evaluations=10000,
    population_size=74,
    num_clusters=11,
    crossover_ratio=0.6820,
    mutation_ratio=0.1498,
    tour_size=6,
    local_search=local_search_algorithms$pls,
    neighborhood= 0.5099,
    ls_budget= 65.45,
    ls_pos=2,
    ls_params = list(
      acceptance_criteria_fn=helper.non.dominated,
      rank_cutoff = 3
    )
  ),
  kegg = list(
    evaluations=10000,
    population_size=40,
    num_clusters=12,
    crossover_ratio=0.7950,
    mutation_ratio=0.0714,
    tour_size=6,
    local_search=local_search_algorithms$lmols,
    neighborhood=0.3214,
    ls_budget=44.51,
    ls_pos=2
  ),
  disgenet_dis = list(
    evaluations=10000,
    population_size=61,
    num_clusters=10,
    crossover_ratio=0.9372,
    mutation_ratio=0.1649,
    tour_size=5,
    local_search=local_search_algorithms$pls,
    neighborhood=0.3809,
    ls_budget=73.03,
    ls_pos=2,
    ls_params = list(
      acceptance_criteria_fn=helper.non.dominated,
      rank_cutoff = 1
    )
    
  ),
  base = list(
    evaluations=10000,
    population_size=40,
    num_clusters=12,
    crossover_ratio=0.70,
    mutation_ratio=0.07,
    tour_size=2,
    local_search=NULL
  )
)

calculate.results <- function(datasets_to_process, debug=FALSE) {
  for(process_dataset in datasets_to_process) {
    dataset <- datasets[[process_dataset]]
    for (type in names(best_params)) {
      if (type == 'base') {
        biological_source <- 'go'
      } else {
        biological_source <- type
      }
      
      dmatrix_expression = expression.matrix(NULL, dataset=dataset$name)
      dmatrix_biological = biological.matrix(gene_list, biological_databases[[biological_source]], dataset=dataset$name) 
      if (debug) {
        message(paste("Processing", dataset$name, type, '...'))
      }
      params <- best_params[[type]]
      params$dmatrix_expression <- dmatrix_expression
      params$dmatrix_biological <- dmatrix_biological
      params$debug <- FALSE
      
      save.metaheuristic.results(dataset$name, type, nsga2.custom, params, runs = 13, debug = debug)
    }
  }
}

# calculate.results(list('GSE89116', 'GSE53757'), debug=TRUE)
# calculate.results(list('GSE31189', 'GSE50161', 'GSE6919_U95Av2'), debug=TRUE)

profiler <- function(fn, times=1000)
{
  start.time <- Sys.time()
  for (i in 1:times) {
    fn()
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  print(time.taken)
  print(paste("Average execution time:", time.taken / times))
}

source("moc_gapbk.r")

moc.gapbk.evaluate <- function(dataset_key = 'GSE6919_U95Av2', bio = 'go', local_search=TRUE, pop_size=10) {
  start.time <- Sys.time()
  
  dataset <- datasets[[dataset_key]]
  
  dmatrix_expression = expression.matrix(NULL, dataset=dataset$name)
  dmatrix_biological = biological.matrix(NULL, biological_databases[[bio]], dataset=dataset$name)
  
  message(paste("Running moc.gapbk over dataset:", dataset_key, bio))
  results <- moc.gapbk(dmatrix_expression, dmatrix_biological, 10, local_search=local_search, generation=50, pop_size=pop_size)
  saveRDS(results, str_interp("cache/moc_gapbk_${dataset$name}_${bio}_results_pop${pop_size}_g50_ls.rds"))
  metrics <- evaluator.multiobjective.clustering(results, dmatrix_expression, debug=TRUE)
  saveRDS(metrics, str_interp("cache/moc_gapbk_${dataset$name}_${bio}_metrics_pop${pop_sizex}_g50_ls.rds"))
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
}

#moc.gapbk.evaluate('GSE89116')
#moc.gapbk.evaluate('GSE53757')
#moc.gapbk.evaluate('GSE31189')
#moc.gapbk.evaluate('GSE50161')
#moc.gapbk.evaluate('GSE6919_U95Av2')