source("package_installer.r")
source("nsga2.r")
source("matrices_evaluation.r")
source("globals.r")
source("evaluator.r")
source("file_utils.r")

datasets <- list(
  GSE53757 = list(
    name = "Renal_GSE53757",
    chip = "GPL570",
    type = "evaluation",
    cancer = "Renal"
  ),
  GSE41657 = list(
    name = "Colorectal_GSE41657",
    chip = "GPL6480",
    type = "evaluation",
    cancer = "Colorectal"
  ),
  GSE70947 = list(
    name = "Breast_GSE70947",
    chip = "GPL13607",
    type = "evaluation",
    cancer = "Breast"
  ),
  GSE28497 = list(
    name = "Leukemia_GSE28497",
    chip = "GPL96",
    type = "evaluation",
    cancer = "Leukemia"
  ),
  GSE6919_U95C = list(
    name = "Prostate_GSE6919_U95C",
    chip = "GPL8300",
    type = "evaluation",
    cancer = "Prostate"
  ),
  GSE6919_U95B = list(
    name = "Prostate_GSE6919_U95B",
    chip = "GPL8300",
    type = "training",
    cancer = "Prostate"
  ),
  GSE6919_U95Av2 = list(
    name = "Prostate_GSE6919_U95Av2",
    chip = "GPL8300",
    type = "training",
    cancer = "Prostate"
  ),
  GSE22405 = list(
    name = "Liver_GSE22405",
    chip = "GPL10553",
    type = "training",
    cancer = "Liver"
  ),
  GSE31189 = list(
    name = "Bladder_GSE31189",
    chip = "GPL570",
    type = "training",
    cancer = "Bladder"
  ),
  GSE63459 = list(
    name = "Lung_GSE63459",
    chip = "GPL6883",
    type = "training",
    cancer = "Lung"
  ),
  GSE71449 = list(
    name = "Leukemia_GSE71449",
    chip = "GPL19197",
    type = "training",
    cancer = "Leukemia"
  ),
  GSE6008 = list(
    name = "Ovary_GSE6008",
    chip = "GPL96",
    type = "training",
    cancer = "Ovary"
  ),
  GSE50161 = list(
    name = "Brain_GSE50161",
    chip = "GPL570",
    type = "training",
    cancer = "Brain"
  ),
  GSE44861 = list(
    name = "Colorectal_GSE44861",
    chip = "GPL3921",
    type = "training",
    cancer = "Colorectal"
  ),
  GSE77953 = list(
    name = "Colorectal_GSE77953",
    chip = "GPL96",
    type = "training",
    cancer = "Colorectal"
  ),
  GSE10797 = list(
    name = "Breast_GSE10797",
    chip = "GPL571",
    type = "training",
    cancer = "Breast"
  ),
  GSE26304 = list(
    name = "Breast_GSE26304",
    chip = "GPL6848",
    type = "training",
    cancer = "Breast"
  ),
  GSE45827 = list(
    name = "Breast_GSE45827",
    chip = "GPL570",
    type = "training",
    cancer = "Breast"
  ),
  GSE7904 = list(
    name = "Breast_GSE7904",
    chip = "GPL570",
    type = "training",
    cancer = "Breast"
  ),
  GSE59246 = list(
    name = "Breast_GSE59246",
    chip = "GPL13607",
    type = "training",
    cancer = "Breast"
  )
)


run.default <- function() {
  test_file_name <- 'Renal_GSE53757'
  test_file <- paste("data/evaluation/", test_file_name, '.csv.gz', sep="")
  data <- read.dataset(test_file)
  
  data <- clean.and.translate.entrez.id(data, "GPL570")
  data <- data[ , order( as.numeric(colnames(data)) ) ]
  
  gene_list <- colnames(data)
  dmatrix_expression = expression.matrix(t(data), dataset=test_file_name)
  dmatrix_biological = biological.matrix(gene_list, biological_databases$go, dataset=test_file_name)
  
  #set.seed(1048)
  #Rprof("profile1.out", line.profiling=TRUE, memory.profiling=TRUE)
  results <- nsga2.custom(dmatrix_expression, dmatrix_biological, population_size = 80, generations = 40, num_clusters = 5, ls_pos=NULL, local_search = NULL)
  #Rprof(NULL)
  
  #summaryRprof("profile1.out", lines = "show")
}

#' Used along run.r.sample.test to check if there is any bias in R's sample method.
#' 
#' @param data Array of data to be tested for random selections
#' @param times Integer. How many times is sample() run.
#' @param test_n Integer. Test n° identifier. Used to save png to different file.
#' @param test_type String. What kind of test. Used to save png to different file.
#' 
r.sample.test <- function(data, times, test_n, test_type) {
  appearances <- hash()
  for (i in data) { appearances[[paste(i)]] <- 0 }
  
  for (i in 1:times) {
    n <- sample(data, 1)
    appearances[[paste(n)]] <- appearances[[paste(n)]] + 1
  }
  
  data.string <- str_interp("${min(data)}:${max(data)}")
  png(filename=str_interp("plots/r_${test_type}_sample_test_${min(data)}_to_${max(data)}_executed_${times}_times_${test_n}.png"))
  plot(as.list(appearances), data, xlab=str_interp("N° apariciones en sample(${data.string}, 1)"), ylab=data.string, main=str_interp("Prueba de distribución de sample de R (${times} veces)."))
  dev.off()
}

#' Run multiple sample tests and save results in png files 
#' under the 'plots' directory
#' 
run.r.sample.test <- function() {
  
  set.seed(12521)
  
  for(big_test in 1:3) {
    r.sample.test(1:10000, 10000000, big_test, "big")
  }
  
  for(small_test in 1:10) {
    r.sample.test(1:1000, 1000000, small_test, "small")
  }
}

#' Load a dataset and translate its genes (colnames) to ENTREZ_GENE_ID
#' 
#' @param dataset A list with the attributes: name (file_name without extension),
#'   chip (Name of the platform used, example: "GPL_570"), type ("evaluation" or
#'   "training") and cancer
load.dataset <- function(dataset) {
  
  file_name <- paste(dataset$name, '.csv.gz', sep='')
  file_path <- paste('data/', dataset$type, '/', file_name, sep="")
  data <- read.dataset(file_path)
  
  data <- clean.and.translate.entrez.id(data, dataset$chip)
  data <- data[ , order( as.numeric(colnames(data)) ) ]
}

#' Calculate all biological distance matrices and store them in .rda files under
#' the 'cache' directory to save time in the future.
#' 
precalculate.biological.dmatrix <- function() {
  plan(multisession, gc = TRUE, workers = 22)
  
  for (dataset in datasets) {
    print(dataset)
    data <- load.dataset(dataset)
    gene_list <- colnames(data)
    for (biological_source in biological_databases) {
      future({biological.matrix(gene_list, biological_source, dataset=dataset$name)})
    }
  }
}

