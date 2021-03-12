source("package_installer.r")
source("nsga2.r")
source("matrices_evaluation.r")
source("globals.r")
source("evaluator.r")
source("file_utils.r")

test_file_name <- 'Leukemia_GSE28497.csv-2-all-sample-5-percent.csv'
test_file <- paste("data/training/samples/", test_file_name, sep="")
data <- read.dataset(test_file)

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

data <- clean.and.translate.entrez.id(data, "GPL96")
data <- data[ , order( as.numeric(colnames(data)) ) ]

gene_list <- colnames(data)
dmatrix_expression = expression.matrix(t(data), dataset=test_file_name)
dmatrix_biological = biological.matrix(gene_list, biological_databases$go, dataset=test_file_name)

#set.seed(1048)
#Rprof("profile1.out", line.profiling=TRUE, memory.profiling=TRUE)
results <- nsga2.custom(dmatrix_expression, dmatrix_biological, population_size = 80, generations = 40, num_clusters = 5, ls_pos=NULL, local_search = NULL)
#Rprof(NULL)

#summaryRprof("profile1.out", lines = "show")
