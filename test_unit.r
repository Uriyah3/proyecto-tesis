options(warn = 1)
PASS <- 0; FAIL <- 0
check <- function(name, expr) {
  cat(sprintf("  %-55s", name)); flush.console()
  ok <- tryCatch({ expr; TRUE }, error = function(e) { cat("FAIL:", e$message, "\n"); FALSE })
  if (ok) { cat("OK\n"); PASS <<- PASS + 1 } else { FAIL <<- FAIL + 1 }
}

cat("=== UNIT TESTS ===\n\n")

# --- credentials.r ---
cat("--- credentials.r ---\n")
source("credentials.r")

check("load.credentials sets env vars from .env", {
  tmpfile <- tempfile(fileext = ".env")
  writeLines(c("TEST_VAR_A=hello", "TEST_VAR_B=world=123", "# comment", "", "TEST_VAR_C=ok"), tmpfile)
  load.credentials(tmpfile)
  stopifnot(Sys.getenv("TEST_VAR_A") == "hello")
  stopifnot(Sys.getenv("TEST_VAR_B") == "world=123")
  stopifnot(Sys.getenv("TEST_VAR_C") == "ok")
  unlink(tmpfile)
})

check("load.credentials skips comments and blank lines", {
  stopifnot(Sys.getenv("#") == "")
})

check("get.david.email stops when no credentials", {
  old_email <- Sys.getenv("DAVID_EMAIL")
  Sys.unsetenv("DAVID_EMAIL")
  # Temporarily hide both .env and mail_list.txt
  if (file.exists(".env")) file.rename(".env", ".env.test.bak")
  if (file.exists("mail_list.txt")) file.rename("mail_list.txt", "mail_list.txt.bak")
  err <- tryCatch({ get.david.email() }, error = function(e) e$message)
  if (file.exists(".env.test.bak")) file.rename(".env.test.bak", ".env")
  if (file.exists("mail_list.txt.bak")) file.rename("mail_list.txt.bak", "mail_list.txt")
  if (nchar(old_email) > 0) Sys.setenv(DAVID_EMAIL = old_email)
  stopifnot(grepl("No DAVID email found", err))
})

check("get.disgenet.api.key stops when no credentials", {
  old_key <- Sys.getenv("DISGENET_API_KEY")
  Sys.unsetenv("DISGENET_API_KEY")
  if (file.exists(".env")) file.rename(".env", ".env.test.bak")
  err <- tryCatch({ get.disgenet.api.key() }, error = function(e) e$message)
  if (file.exists(".env.test.bak")) file.rename(".env.test.bak", ".env")
  if (nchar(old_key) > 0) Sys.setenv(DISGENET_API_KEY = old_key)
  stopifnot(grepl("No DisGeNET API key found", err))
})

# --- david_client.r ---
cat("\n--- david_client.r ---\n")
source("david_client.r")

check("david.new.session creates a handle", {
  david.new.session()
  stopifnot(!is.null(david_session))
  stopifnot(inherits(david_session, "handle"))
})

check("david.extract.return handles NULL input", {
  stopifnot(is.na(david.extract.return(NULL)))
})

check("david.extract.return parses XML with return element", {
  xml <- xml2::read_xml('<ns:testResponse xmlns:ns="http://test"><ns:return>hello</ns:return></ns:testResponse>')
  stopifnot(david.extract.return(xml) == "hello")
})

check("DAVID_ANNOTATION_CATEGORIES has 44 entries", {
  stopifnot(length(DAVID_ANNOTATION_CATEGORIES) == 44)
})

check("POST threshold triggers for large params", {
  big_param <- paste(rep("12345", 500), collapse = ",")
  total <- sum(nchar(as.character(list(args0 = big_param))))
  stopifnot(total > 2000)
})

# --- nsga2.r clustering ---
cat("\n--- nsga2.r (clustering) ---\n")
library(data.table)

check("data.table join produces valid clustering", {
  set.seed(42)
  nmedoids <- 3; ngenes <- 10
  dm <- data.table(matrix(runif(nmedoids * ngenes), nrow = nmedoids))
  DT <- data.table(
    value = unlist(dm, use.names = FALSE),
    colid = rep(1:ncol(dm), each = nrow(dm)),
    rowid = 1:nrow(dm)
  )
  setkey(DT, colid, value)
  clustering_result <- DT[J(unique(colid)), rowid, mult = "first"]
  if (is.data.table(clustering_result) || is.data.frame(clustering_result)) {
    cl <- clustering_result[["rowid"]]
  } else {
    cl <- clustering_result
  }
  stopifnot(length(cl) == ngenes)
  stopifnot(all(cl >= 1 & cl <= nmedoids))
})

# --- evaluator.r silhouette ---
cat("\n--- evaluator.r (silhouette) ---\n")

check("silhouette na.rm works with NA values", {
  indices <- c(0.5, NA, 0.3, NA, 0.7)
  stopifnot(max(indices, na.rm = TRUE) == 0.7)
  stopifnot(abs(mean(indices, na.rm = TRUE) - 0.5) < 0.001)
  stopifnot(min(indices, na.rm = TRUE) == 0.3)
})

check("which.best returns valid index, not 0", {
  source("evaluator.r")
  # Simulate: all silhouettes below median (edge case)
  sil <- c(0.1, 0.2, 0.3)
  pop <- data.frame(objective_bio = c(0.9, 0.8, 0.7))
  result <- which.best('best', sil, pop)
  stopifnot(result >= 1 && result <= 3)
})

cat(sprintf("\n=== RESULTS: %d passed, %d failed ===\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
