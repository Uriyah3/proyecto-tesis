#!/usr/bin/env Rscript
library("optparse")
source('main.r')

option_list = list(
  make_option(c("--seed"), type = "integer",
              help = "random seed"),
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
              help = "neighborhood"),
  make_option(c("--local_search"), type = "character",
              help = "local_search"),
  make_option(c("--ls_budget"), type = "double",
              help = "local_search_budget"),
  make_option(c("--ls_pos"), type = "character",
              help = "ls_pos"),
  make_option(c("--acc_fn"), type = "character",
              help = "acceptance_criteria_fn"),
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

dataset <- tools::file_path_sans_ext(tools::file_path_sans_ext(basename(opt$input)))
dataset <- str_split_fixed(dataset, "_", 2)[2]
dataset <- datasets[[dataset]]

data <- load.dataset(dataset)

gene_list <- colnames(data)
dmatrix_expression = expression.matrix(t(data), dataset=dataset$name)
dmatrix_biological = biological.matrix(gene_list, biological_databases[[opt$biological_source]], dataset=dataset$name)

rank_cutoff <- NULL
if (!is.null(opt$pls_rank_cutoff)) {
  rank_cutoff <- opt$pls_rank_cutoff
} else if (!is.nul(opt$ce_rank_cutoff)) {
  rank_cutoff <- opt$ce_rank_cutoff
} else if (!is.nul(opt$mosa_rank_cutoff)) {
  rank_cutoff <- opt$mosa_rank_cutoff
}

params <- list(dmatrix_expression=dmatrix_expression, dmatrix_biological=dmatrix_biological, num_clusters=opt$num_clusters, evaluations=opt$evaluations, population_size=opt$population, crossover_ratio=opt$crossover, crossover_prob=opt$crossover_prob, mutation_ratio=opt$mutation, tour_size=opt$tour_size, neighborhood = opt$neighborhood, local_search=local_search_algorithms[[opt$local_search]], ls_pos=opt$ls_pos, ls_budget=opt$ls_budget, debug=opt$debug, ls_params=list(acceptance_criteria_fn=get(opt$acc_fn), rank_cutoff=rank_cutoff, alfa=opt$alfa))

results <- evaluator.metaheuristics(nsga2.custom, params)

cat(results$mean_results$centered_hypervolume)