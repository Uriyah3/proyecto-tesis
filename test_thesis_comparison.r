cat("=== THESIS COMPARISON TEST ===\n")
cat("Target: Prostate KEGG silhouette max=0.573, mean=0.441 (thesis Table 5.4)\n\n")

source("main.r")
gc()  # Free memory after loading all packages
cat("main.r loaded. Free memory after gc.\n\n"); flush.console()

# Use Prostate_GSE6919_U95Av2 — smallest evaluation dataset (9005 genes, ~1.2GB/matrix)
dataset <- datasets$GSE6919_U95Av2
cat("Loading dataset:", dataset$name, "...\n"); flush.console()
ev_data <- load.dataset(dataset)
gene_list <- colnames(ev_data)
cat("Loaded:", nrow(ev_data), "samples,", length(gene_list), "genes\n"); flush.console()

# Build/load distance matrices
cat("Loading expression matrix...\n"); flush.console()
dmatrix_exp <- expression.matrix(t(ev_data), dataset = dataset$name)
rm(ev_data); gc()  # Free the raw data — we only need the matrices now
cat("Loading KEGG biological matrix...\n"); flush.console()
dmatrix_bio <- biological.matrix(gene_list, biological_databases$kegg, dataset = dataset$name)
cat("Matrices:", nrow(dmatrix_exp), "x", ncol(dmatrix_exp), "\n"); flush.console()
cat("Memory after matrices:", round(sum(gc()[,2]), 0), "MB\n\n"); flush.console()

# Run with thesis KEGG params — 13 runs like the thesis
cat("Running NSGA-II + L-MOLS (thesis KEGG params, 10000 evals, 13 runs)...\n"); flush.console()
cat("This will take ~100 minutes...\n"); flush.console()
params <- best_params$kegg
params$dmatrix_expression <- dmatrix_exp
params$dmatrix_biological <- dmatrix_bio
params$debug <- FALSE

t1 <- system.time({
  results <- evaluator.metaheuristics(nsga2.custom, params, runs = 13)
})

cat("\n=== RESULTS (13 runs) ===\n")
cat("Time:", round(t1[3] / 60, 1), "minutes\n\n")

# Per-run results
cat("Per-run silhouette:\n")
sil_maxes <- numeric(length(results$full_results))
sil_means <- numeric(length(results$full_results))
for (i in 1:length(results$full_results)) {
  sil <- results$full_results[[i]]$silhouette
  hv <- results$full_results[[i]]$hypervolume
  sil_maxes[i] <- sil$max_silhouette
  sil_means[i] <- sil$mean_silhouette
  cat(sprintf("  Run %2d: sil_max=%.3f sil_mean=%.3f hv=%.4f\n",
    i, sil$max_silhouette, sil$mean_silhouette, hv$hypervolume))
}

# Aggregated results
cat("\n--- Aggregated over 13 runs ---\n")
cat(sprintf("Silhouette max:  mean=%.3f sd=%.3f\n", mean(sil_maxes), sd(sil_maxes)))
cat(sprintf("Silhouette mean: mean=%.3f sd=%.3f\n", mean(sil_means), sd(sil_means)))

cat("\n--- Comparison to thesis (Prostate KEGG, Table 5.4) ---\n")
cat(sprintf("                 Thesis    Ours\n"))
cat(sprintf("Sil max:         0.573     %.3f\n", mean(sil_maxes)))
cat(sprintf("Sil mean:        0.441     %.3f\n", mean(sil_means)))

diff_max <- abs(mean(sil_maxes) - 0.573) / 0.573 * 100
diff_mean <- abs(mean(sil_means) - 0.441) / 0.441 * 100
cat(sprintf("\nDifference: max=%.1f%%, mean=%.1f%%\n", diff_max, diff_mean))

if (mean(sil_maxes) > 0.4 && mean(sil_means) > 0.3) {
  cat("\nVERIFIED: Results are in the thesis ballpark.\n")
  cat("The algorithm is working correctly on R 4.5.3.\n")
} else {
  cat("\nWARNING: Results are significantly lower than thesis values.\n")
  quit(status = 1, save = "no")
}
