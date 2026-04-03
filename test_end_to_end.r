options(warn = 1)
PASS <- 0
FAIL <- 0

check <- function(name, expr) {
  cat(sprintf("  %-50s", name)); flush.console()
  result <- tryCatch({ expr; TRUE }, error = function(e) { cat("ERROR:", e$message, "\n"); FALSE })
  if (result) { cat("OK\n"); PASS <<- PASS + 1 } else { FAIL <<- FAIL + 1 }
}

cat("=== END-TO-END TEST SUITE ===\n")
cat("R version:", R.version.string, "\n\n")

# -------------------------------------------------------
cat("--- Phase 1: File loading ---\n")
# -------------------------------------------------------

check("source nsga2.r", source("nsga2.r"))
check("source matrices_evaluation.r", source("matrices_evaluation.r"))
check("source evaluator.r", source("evaluator.r"))
check("source globals.r", source("globals.r"))
check("source file_utils.r", source("file_utils.r"))
check("source david_client.r", source("david_client.r"))
check("library ComplexHeatmap", library(ComplexHeatmap))
check("library circlize", library(circlize))
check("library magick", library(magick))
check("library future", library(future))

# -------------------------------------------------------
cat("\n--- Phase 2: C++ compilation ---\n")
# -------------------------------------------------------

check("sourceCpp performance.cpp", {
  Rcpp::sourceCpp("performance.cpp")
  stopifnot(exists("Rmcalc"))
  stopifnot(exists("clusteringCalc"))
})

# -------------------------------------------------------
cat("\n--- Phase 3: Data loading & gene translation ---\n")
# -------------------------------------------------------

check("read sample dataset", {
  test_data <<- read.dataset("data/training/samples/Leukemia_GSE28497.csv-1-all-sample-5-percent.csv")
  stopifnot(nrow(test_data) > 0, ncol(test_data) > 0)
})

check("translate gene IDs (GPL96)", {
  test_data <<- clean.and.translate.entrez.id(test_data, "GPL96")
  test_data <<- test_data[, order(as.numeric(colnames(test_data)))]
  stopifnot(ncol(test_data) > 500)
})

# -------------------------------------------------------
cat("\n--- Phase 4: Distance matrix computation ---\n")
# -------------------------------------------------------

test_file_name <- "Leukemia_GSE28497.csv-1-all-sample-5-percent"
gene_list <- colnames(test_data)

check("expression matrix (Pearson)", {
  dmatrix_exp <<- expression.matrix(t(test_data), dataset = test_file_name)
  stopifnot(nrow(dmatrix_exp) > 0)
})

check("biological matrix (GO)", {
  dmatrix_bio <<- biological.matrix(gene_list, biological_databases$go, dataset = test_file_name)
  stopifnot(nrow(dmatrix_bio) > 0)
})

# -------------------------------------------------------
cat("\n--- Phase 5: NSGA-II algorithm ---\n")
# -------------------------------------------------------

check("NSGA-II + PLS local search (1 run, 1400 evals)", {
  test_params <<- list(
    dmatrix_expression = dmatrix_exp,
    dmatrix_biological = dmatrix_bio,
    population_size = 80,
    evaluations = 1400,
    num_clusters = 5,
    ls_pos = 1,
    local_search = local_search_algorithms$pls,
    neighborhood = 0.7,
    debug = FALSE,
    ls_params = list(rank_cutoff = 3, acceptance_criteria_fn = helper.non.dominated)
  )
  results <<- evaluator.metaheuristics(nsga2.custom, test_params, runs = 1)
  stopifnot(!is.null(results))
})

# -------------------------------------------------------
cat("\n--- Phase 6: Evaluation metrics ---\n")
# -------------------------------------------------------

check("silhouette metric computed", {
  sil <- results$full_results[[1]]$silhouette
  stopifnot(!is.null(sil$max_silhouette))
  cat(sprintf("[max=%.3f] ", as.numeric(sil$max_silhouette)))
})

check("hypervolume metric computed", {
  hv <- results$full_results[[1]]$hypervolume
  stopifnot(!is.null(hv$hypervolume))
  cat(sprintf("[hv=%.4f] ", hv$hypervolume))
})

# -------------------------------------------------------
cat("\n--- Phase 7: DAVID client ---\n")
# -------------------------------------------------------

check("david_client.r functions exist", {
  stopifnot(exists("david.authenticate"))
  stopifnot(exists("david.add.list"))
  stopifnot(exists("david.set.categories"))
  stopifnot(exists("david.get.term.cluster.report"))
  stopifnot(exists("david.annotate.genes"))
  stopifnot(exists("DAVID_ANNOTATION_CATEGORIES"))
  stopifnot(length(DAVID_ANNOTATION_CATEGORIES) == 44)
})

check("david.new.session() creates handle", {
  david.new.session()
  stopifnot(!is.null(david_session))
})

check("DAVID API reachable (authenticate attempt)", {
  # This will return FALSE because the test email isn't registered,
  # but it should NOT throw an HTTP error — that proves the API is up.
  result <- tryCatch({
    david.authenticate("test@example.org")
  }, error = function(e) {
    stop(paste("HTTP connection failed:", e$message))
  })
  # FALSE is expected (unregistered email), error would mean API is down
  cat("[auth=", result, "] ", sep = "")
})

source("credentials.r")
email <- tryCatch(get.david.email(), error = function(e) { cat("No email configured:", e$message, "\n"); NULL })
if (is.null(email)) stop("Configure DAVID_EMAIL in .env first")
check(paste0("DAVID auth with mail_list.txt (", email, ")"), {
  david.new.session()
  result <- david.authenticate(email)
  cat(sprintf("[auth=%s] ", result))
  if (result) {
    # If auth works, try the full pipeline
    test_genes <- c("6122", "481", "7019", "212", "4739")
    mapped <- david.add.list(test_genes, "ENTREZ_GENE_ID", "e2e_test", 0L)
    cat(sprintf("[mapped=%s] ", mapped))
    david.set.categories(c("GOTERM_BP_ALL", "KEGG_PATHWAY"))
    clusters <- david.get.term.cluster.report()
    cat(sprintf("[clusters=%d] ", length(clusters)))
  } else {
    cat("[skipping API tests - email not registered] ")
  }
})

# -------------------------------------------------------
cat("\n\n=== RESULTS ===\n")
cat(sprintf("PASSED: %d  FAILED: %d\n", PASS, FAIL))
if (FAIL == 0) {
  cat("ALL TESTS PASSED!\n")
} else {
  cat("SOME TESTS FAILED.\n")
  quit(status = 1)
}
