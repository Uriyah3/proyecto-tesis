if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("GOSemSim")
#BiocManager::install("DOSE")
#X BiocManager::install("biomaRt")
#X BiocManager::install("pathifier")
#BiocManager::install("meshes")
#BiocManager::install("MeSH.Hsa.eg.db")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("KEGGREST")
#X BiocManager::install("KEGG.db")
#BiocManager::install("BioCor")
#BiocManager::install("STRINGdb")
#X library(pathifier)
#X install.packages("neo4jshell")
#BiocManager::install("qusage")

source("nsga2.r")
source("matrices_evaluation.r")
source("globals.r")
source("evaluator.r")
source("file_utils.r")

test_file_name <- 'Leukemia_GSE28497.csv-2-all-sample-5-percent.csv'
test_file <- paste("data/training/samples/", test_file_name, sep="")
data <- read.dataset(test_file)

datasets <- list(
  Leukemia = list(
    name = "Leukemia_GSE28497",
    chip = "GPL96",
    type = "training"
  ),
  Liver = list(
    name = "Liver_GSE14520_U133A",
    chip = "GPL3921",
    type = "training"
  ),
  Renal = list(
    name = "Renal_GSE53757",
    chip = "GPL570",
    type = "training"
  ),
  Breast = list(
    name = "Breast_GSE70947",
    chip = "GPL13607",
    type = "evaluation"
  ),
  Colorectal = list(
    name = "Colorectal_GSE44076",
    chip = "GPL13667",
    type = "evaluation"
  ),
  Lung = list(
    name = "Lung_GSE19804",
    chip = "GPL570",
    type = "evaluation"
  ),
  Prostate = list(
    name = "Prostate_GSE46602",
    chip = "GPL570",
    type = "evaluation"
  )
)

data <- clean.and.translate.entrez.id(data, "GPL96")
data <- data[ , order( as.numeric(colnames(data)) ) ]

gene_list <- colnames(data)
dmatrix_expression = expression.matrix(t(data), dataset=test_file_name)
dmatrix_biological = biological.matrix(gene_list, biological_databases$go, dataset=test_file_name)

set.seed(1048)
results <- nsga2.custom(dmatrix_expression, dmatrix_biological, ls_pos = 2, local_search = local_search_algorithms$pls, population_size = 20, generations = 10, neighborhood = 0.80)