# Create sampled data
# Creates 3 sampled versions of the datasets contained in data/training:
#   (row_sample) 10% or row
#   (col_sample) 10% of cols
#   (all_sample) 10% of rows and cols
# *The 10% can be configured, a minimum sample size can also be defined
# Can create various samples from the same dataset
# Stores the generated files in data/training/samples

# % that will be sampled from rows and columns
PERCENTAGE_SAMPLING <- 10
# Number of samples that are generated for each file
SAMPLES_PER_FILE <- 40
# Minimum amount of biological samples
MINIMUM_ROW_SIZE <- 50
# Minimum amount of genes
MINIMUM_COL_SIZE <- 2500

###### Function definitions ######

# Load comma separated text file named dataFile
readFile <- function(dataFile) {
  
  data <- read.csv(dataFile, header = TRUE, sep = ",", quote= "", dec = ".", row.names = 1)
}

# Generates and saves sample_number samples from data
generateSamples <- function(dataFile, sample_number) {
  data <- readFile(dataFile)
  outputFile <- tools::file_path_sans_ext(basename(dataFile))
  outputFile <- paste(sampleDir, "/", outputFile, sep = "")
  
  sampleDirSanityCheck()
  for (i in 1:sample_number) {
    rows <- nrow(data)
    cols <- ncol(data)
    
    bounded_rows <- round(rows * (PERCENTAGE_SAMPLING / 100))
    bounded_rows <- ifelse(MINIMUM_ROW_SIZE > bounded_rows, MINIMUM_ROW_SIZE, bounded_rows)
    bounded_rows <- ifelse(bounded_rows > rows, rows, bounded_rows)
    bounded_cols <- round(cols * (PERCENTAGE_SAMPLING / 100))
    bounded_cols <- ifelse(MINIMUM_COL_SIZE > bounded_cols, MINIMUM_COL_SIZE, bounded_cols)
    bounded_cols <- ifelse(bounded_cols > cols, cols, bounded_cols)
    
    generateAndSaveSample(data, rows, bounded_cols, paste(outputFile, "-" , i, "-col-sample-", PERCENTAGE_SAMPLING, "-percent", ".csv", sep=""))
    generateAndSaveSample(data, bounded_rows, cols, paste(outputFile, "-" , i, "-row-sample-", PERCENTAGE_SAMPLING, "-percent", ".csv", sep=""))
    generateAndSaveSample(data, bounded_rows, bounded_cols, paste(outputFile, "-" , i, "-all-sample-", PERCENTAGE_SAMPLING, "-percent", ".csv", sep=""))
  }
}

# Checks if dir data/training/samples exists, otherwise create it
sampleDirSanityCheck <- function() {
  if (!dir.exists(sampleDir)) {
    dir.create(sampleDir)
  }
}

# Generate a single bounded sample
generateAndSaveSample <- function(data, rows, cols, outputFile) {
  sampled_data = data[sample(nrow(data), rows), sample(ncol(data), cols)]
  write.csv(sampled_data, outputFile, row.names = TRUE, quote = FALSE)
}

###### Main loop ######
runSampling <- function() {
  dataDir <- "data/training/"
  dataFiles <- list.files(dataDir)
  dataFiles <- paste0(dataDir, dataFiles)
  sampleDir <- paste(dataDir, "samples", sep="")
  
  for (dataFile in dataFiles) {
    generateSamples(dataFile, SAMPLES_PER_FILE)
  }  
}
