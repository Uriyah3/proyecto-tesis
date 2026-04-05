cat("=== PIPELINE TEST ===\n")
cat("Tests: calculate.results -> evaluate.results -> show.results -> friedman\n\n")

source("main.r")
gc()

# Override params for speed (2000 evals, not 10000)
best_params$kegg$evaluations <- 2000
# Only test KEGG (skip go, string, disgenet_dis, base)
original_best_params <- best_params
best_params <- list(kegg = best_params$kegg)

cat("Step 1: calculate.results (3 runs, 2000 evals, KEGG only)...\n"); flush.console()
save.metaheuristic.results <- function(dataset.name, type, metaheuristic, params, runs = 13, debug = FALSE) {
  # Override to use 3 runs regardless of what calculate.results passes
  runs <- 3
  for (n in 1:runs) {
    if (debug) message(paste("Run", n, "/", runs))
    results <- list()
    results$nsga <- do.call(metaheuristic, params)
    results$time <- 0
    filename <- build.saved.results.filename(dataset.name, type, n)
    saveRDS(results, filename)
    gc()
  }
}

calculate.results(list('GSE6919_U95Av2'), debug = TRUE)

cat("\nStep 2: evaluate.results...\n"); flush.console()
evaluate.results(list('GSE6919_U95Av2'), runs = 3, debug = TRUE)

cat("\nStep 3: show.results...\n"); flush.console()
show.results(list('GSE6919_U95Av2'), runs = 3, debug = FALSE)

# Step 4: Load results and verify
cat("\nStep 4: Verify cached results...\n"); flush.console()
dataset <- datasets$GSE6919_U95Av2
sil_values <- sapply(1:3, function(i) {
  sil <- load.evaluation.from.cache(dataset$name, "kegg", i, "silhouette")
  sil$max_silhouette
})
cat("Silhouette max per run:", sil_values, "\n")
avg_sil <- mean(sil_values, na.rm = TRUE)
cat("Average silhouette max:", avg_sil, "\n")

# Thesis value for Prostate KEGG: 0.573
if (is.numeric(avg_sil) && !is.na(avg_sil) && avg_sil > 0.3) {
  cat("\nPIPELINE TEST PASSED\n")
  cat("Results are valid numeric values in acceptable range.\n")
} else {
  cat("\nPIPELINE TEST FAILED\n")
  cat("Silhouette values are invalid or too low.\n")
  quit(status = 1, save = "no")
}
