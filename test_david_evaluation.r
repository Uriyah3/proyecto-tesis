cat("=== DAVID BIOLOGICAL EVALUATION TEST ===\n")
cat("Target: Prostate KEGG Table 5.6\n")
cat("  DAVID #groups: 4157\n")
cat("  DAVID max enrichment: 153.918\n")
cat("  DAVID mean enrichment: 3.014\n\n")

source("main.r")
gc()

dataset <- datasets$GSE6919_U95Av2
ev_data <- load.dataset(dataset)
gene_list <- colnames(ev_data)
dmatrix_exp <- expression.matrix(t(ev_data), dataset = dataset$name)
dmatrix_bio <- biological.matrix(gene_list, biological_databases$kegg, dataset = dataset$name)
rm(ev_data); gc()

# Run NSGA-II once with thesis KEGG params
cat("Running NSGA-II + L-MOLS (10000 evals, 1 run)...\n"); flush.console()
params <- best_params$kegg
params$dmatrix_expression <- dmatrix_exp
params$dmatrix_biological <- dmatrix_bio
params$debug <- FALSE

nsga_results <- do.call(nsga2.custom, params)
cat("Pareto front:", nrow(nsga_results$population), "solutions\n"); flush.console()

# Run full evaluation WITH DAVID biological validation
cat("\nRunning evaluator.multiobjective.clustering (with DAVID)...\n")
cat("This calls DAVID for each cluster in the best solution.\n")
cat("May take 15-30 minutes depending on DAVID responsiveness.\n\n"); flush.console()

t1 <- system.time({
  full_eval <- tryCatch(
    evaluator.multiobjective.clustering(nsga_results, dmatrix_exp, debug = TRUE,
      dataset_name = dataset$name, bio = "kegg", iter = 1),
    error = function(e) { cat("ERROR:", e$message, "\n"); NULL }
  )
})

if (!is.null(full_eval)) {
  cat("\n=== RESULTS ===\n")
  cat("Time:", round(t1[3] / 60, 1), "minutes\n\n")

  cat("Silhouette max: ", full_eval$silhouette$max_silhouette, "\n")
  cat("Silhouette mean:", full_eval$silhouette$mean_silhouette, "\n")
  cat("Hypervolume:    ", full_eval$hypervolume$hypervolume, "\n")

  if (!is.null(full_eval$biological)) {
    cat("\nDAVID Results:\n")
    cat("  #groups:         ", full_eval$biological$cluster_count, "\n")
    cat("  max enrichment:  ", full_eval$biological$max_enrichment, "\n")
    cat("  mean enrichment: ", full_eval$biological$mean_enrichment, "\n")

    cat("\n--- Comparison to thesis (Prostate KEGG, Table 5.6) ---\n")
    cat(sprintf("                  Thesis     Ours\n"))
    cat(sprintf("#groups:          4157       %s\n", full_eval$biological$cluster_count))
    cat(sprintf("max enrichment:   153.918    %s\n", round(full_eval$biological$max_enrichment, 3)))
    cat(sprintf("mean enrichment:  3.014      %s\n", round(full_eval$biological$mean_enrichment, 3)))
  } else {
    cat("\nDAVID biological results: NULL (DAVID may have failed)\n")
  }
} else {
  cat("\nFull evaluation failed.\n")
  quit(status = 1, save = "no")
}
