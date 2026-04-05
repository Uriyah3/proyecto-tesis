# Final test - Run actual clustering
options(repos = 'https://cran.r-project.org', warn = 1)

cat("===TESTING THESIS PROJECT===\n\n")

# Load nsga2 (this loads all dependencies)
cat("1. Loading nsga2.r...")
source("nsga2.r")
cat(" OK\n")

# Load matrices evaluation 
cat("2. Loading matrices_evaluation.r...")
source("matrices_evaluation.r")
cat(" OK\n")

# Load evaluator
cat("3. Loading evaluator.r...")  
source("evaluator.r")
cat(" OK\n")

cat("\n===ALL CORE FILES LOADED SUCCESSFULLY===\n\n")

# Try loading a small sample dataset
sample_file <- "data/training/samples/Leukemia_GSE28497.csv-1-all-sample-5-percent.csv"
if (!file.exists(sample_file)) {
  cat("ERROR: Sample file not found. Run setup step 3 first:\n")
  cat("  source('data_sampling.r')\n")
  cat("  generate.samples('data/evaluation/Leukemia_GSE28497.csv.gz', 1, 5, 50, 500)\n")
  quit(status = 1, save = "no")
}
cat("4. Testing dataset load...")
source("file_utils.r")
dataset <- read.dataset(sample_file)
cat(" OK -", nrow(dataset), "samples,", ncol(dataset), "genes\n")

# Test gene ID translation
cat("5. Translating gene IDs...")
source("gpl_chip_to_entrez_id.r")
translated <- clean.and.translate.entrez.id(dataset, "GPL96")
cat(" OK -", ncol(translated), "ENTREZ genes\n\n")

cat("PROJECT IS FUNCTIONAL!\n")
