if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("GOSemSim")
#BiocManager::install("biomaRt")
#BiocManager::install("pathifier")
#library(pathifier)

#source("nsga2.r")
source("matrices_evaluation.r")
source("data_sampling.r")

testFile = paste(sampleDir, "/Leukemia_GSE28497-3-col-sample-5-percent.csv", sep="")
data <- readFile(testFile)

datasets = list(
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
  ),
)

data <- clean.and.translate.entrez.id(data, "GPL96")

expressionMatrix = expression.matrix(t(data))
biologicalMatrix = biological.matrix.go(colnames(data))