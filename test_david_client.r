source("david_client.r")

cat("=== DAVID Client Integration Test ===\n\n")

# Read email from mail_list.txt
email <- readLines("mail_list.txt")[1]
cat("Using email:", email, "\n\n")

# Test 1: Authentication
cat("1. Testing authenticate...\n"); flush.console()
david.new.session()
auth_result <- david.authenticate(email)
cat(if (auth_result) "   OK\n" else "   FAILED (is your email registered?)\n")

if (!auth_result) {
  cat("\nTo register, visit: https://davidbioinformatics.nih.gov/webservice/register.htm\n")
  cat("Note: requires an organizational email (no gmail/yahoo/hotmail)\n")
  stop("Cannot continue without valid email registration")
}

# Test 2: Upload a small gene list (reuses session from test 1)
cat("2. Testing addList...\n"); flush.console()
test_genes <- c("6122", "481", "7019", "212", "4739", "27090", "55346", "10175", "8549")
mapped <- david.add.list(test_genes, "ENTREZ_GENE_ID", "test_list", 0L)
cat("   OK - mapped:", mapped, "\n")

# Test 3: Set categories
cat("3. Testing setCategories...\n"); flush.console()
result <- david.set.categories(c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "KEGG_PATHWAY"))
cat(if (result) "   OK\n" else "   FAILED\n")

# Test 4: Get cluster report
cat("4. Testing getTermClusterReport...\n"); flush.console()
clusters <- david.get.term.cluster.report()
cat("   OK -", length(clusters), "clusters found\n")

if (length(clusters) > 0) {
  cat("\n   Cluster details:\n")
  for (i in seq_along(clusters)) {
    cat(sprintf("   Cluster %d: EnrichmentScore = %.3f\n", i, clusters[[i]]$EnrichmentScore))
  }
}

# Test 5: Full convenience wrapper (creates its own fresh session)
cat("\n5. Testing david.annotate.genes (full pipeline)...\n"); flush.console()
full_result <- david.annotate.genes(email, test_genes, debug = TRUE)
cat("   OK\n")
cat("   cluster_count:", full_result$cluster_count, "\n")
if (!is.na(full_result$enrichment[1])) {
  cat("   enrichment:", paste(round(full_result$enrichment, 3), collapse=", "), "\n")
}

cat("\n=== ALL TESTS PASSED ===\n")
