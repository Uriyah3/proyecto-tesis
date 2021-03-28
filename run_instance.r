#!/usr/bin/env Rscript
library("optparse")
library(future.apply)
source('main.r')
options(future.globals.maxSize= 2091289600) # 2GB

option_list = list(
  make_option(c("--seed"), type = "integer",
              help = "random seed"),
  make_option(c("--nbproc"), type = "integer",
              help = "Number of process forks to use", default=3),
  make_option(c("--evaluations"), type = "integer",
              help = "fitness_evaluations"),
  make_option(c("--population"), type = "integer",
              help = "population"),
  make_option(c("--num_clusters"), type = "integer",
              help = "num_clusters"),
  make_option(c("--crossover"), type = "double",
              help = "crossover_ratio"),
  make_option(c("--crossover_prob"), type = "double",
              help = "crossover_prob_modifier"),
  make_option(c("--mutation"), type = "double",
              help = "mutation_ratio"),
  make_option(c("--tour_size"), type = "integer",
              help = "tour_size"),
  make_option(c("--neighborhood"), type = "double",
              help = "neighborhood", default=NULL),
  make_option(c("--local_search"), type = "character",
              help = "local_search"),
  make_option(c("--ls_budget"), type = "double",
              help = "local_search_budget"),
  make_option(c("--ls_pos"), type = "character",
              help = "ls_pos"),
  make_option(c("--acc_fn"), type = "character",
              help = "acceptance_criteria_fn", default="helper.dominates"),
  make_option(c("--pls_rank_cutoff"), type = "integer",
              help = "pls_rank_cutoff"),
  make_option(c("--ce_rank_cutoff"), type = "integer",
              help = "ce_rank_cutoff"),
  make_option(c("--mosa_rank_cutoff"), type = "integer",
              help = "mosa_rank_cutoff"),
  make_option(c("--alfa"), type = "double",
              help = "alfa"),
  make_option(c("-i", "--input"), type = "character",
              help = "Camino al archivo de prueba"),
  make_option(c("-b", "--biological_source"), type = "character",
              help = "Fuente biológica que se debe utilizar"),
  make_option(c("-g", "--gpl"), type = "character",
              help = "Plataforma GPL que utiliza el dataset"),
  make_option(c("-d", "--debug"), type = "logical",
              help = "Imprimir debugging output", default=FALSE)
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
set.seed(opt$seed)

message(opt)

dataset <- tools::file_path_sans_ext(tools::file_path_sans_ext(basename(opt$input)))
dataset <- str_split_fixed(dataset, "_", 2)[2]
dataset <- datasets[[dataset]]

message(dataset)

# No es necesario cargar los datos originales si tengo todas las matrices de 
# distancias precalculadas en archivos .rda
# data <- load.dataset(dataset)
data <- NULL

gene_list <- colnames(data)
dmatrix_expression = expression.matrix(t(data), dataset=dataset$name)
dmatrix_biological = biological.matrix(gene_list, biological_databases[[opt$biological_source]], dataset=dataset$name)
gene_list <- colnames(dmatrix_expression)

rank_cutoff <- NULL
if (!is.null(opt$pls_rank_cutoff)) {
  rank_cutoff <- opt$pls_rank_cutoff
} else if (!is.null(opt$ce_rank_cutoff)) {
  rank_cutoff <- opt$ce_rank_cutoff
} else if (!is.null(opt$mosa_rank_cutoff)) {
  rank_cutoff <- opt$mosa_rank_cutoff
}

neighborhood_matrix <- NULL
if (!is.null(opt$neighborhood)) {
  if(opt$debug) message("Precalculating neighborhood_matrix")
  plan(multicore)
  options(mc.cores=opt$nbproc)
  
  # Save a bit of time precalculating neighborhood_matrix
  dmatrix_combined <- sqrt(dmatrix_expression**2 + dmatrix_biological**2)
  # Find genes that are close to one another
  neighborhood_matrix <- future_sapply(gene_list, function(gene) {
    neighborhood_genes <- dmatrix_combined[gene, , drop=FALSE]
    apply(neighborhood_genes, 1, function(x) colnames(neighborhood_genes)[which(x > 0.000000 & x < opt$neighborhood)] )
  })
  # Such a neighborhood matrix uses too much memory, just consider it as all genes are neighbors
  if (length(unlist(neighborhood_matrix)) > length(gene_list) * (length(gene_list) - 1) *0.9) {
    neighborhood_matrix <- gene_list
  }
  dmatrix_combined <- NULL
  invisible(gc())
}

if (opt$debug) message("Starting metaheuristics evaluation...")

params <- list(dmatrix_expression=dmatrix_expression, dmatrix_biological=dmatrix_biological, num_clusters=opt$num_clusters, evaluations=opt$evaluations, population_size=opt$population, crossover_ratio=opt$crossover, crossover_prob=opt$crossover_prob, mutation_ratio=opt$mutation, tour_size=opt$tour_size, neighborhood = opt$neighborhood, local_search=local_search_algorithms[[opt$local_search]], ls_pos=opt$ls_pos, ls_budget=opt$ls_budget, debug=opt$debug, ls_params=list(acceptance_criteria_fn=get(opt$acc_fn), rank_cutoff=rank_cutoff, alfa=opt$alfa), neighborhood_matrix=neighborhood_matrix)

results <- evaluator.metaheuristics(nsga2.custom, params, debug = opt$debug)

if (opt$debug) {
  message(results)
  message("Returning mean hypervolume (1,1) reference point")
}

cat(results$mean_results$centered_hypervolume * -1)
