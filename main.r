if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("GOSemSim")
#BiocManager::install("biomaRt")
#BiocManager::install("pathifier")
#library(pathifier)

source("nsga2.r")
source("matrices_evaluation")

samplesDir = "data/training/samples/"
testFile = paste(samplesDir, "Leukemia_GSE28497-29-col-sample-10-percent.csv", sep="")
data <- read.csv(testFile, header = TRUE, sep = ",", quote= "", dec = ".", row.names = 1)

chips = list(
  Leukemia = "[HG-U133A] Affymetrix Human Genome U133A Array", #affy_hg_u133a
  Liver = "[HT_HG-U133A] Affymetrix HT Human Genome U133A Array",
  Renal = "[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array", #affy_hg_u133_plus_2
  Breast = "Agilent-028004 SurePrint G3 Human GE 8x60K Microarray", #agilent_sureprint_g3_ge_8x60k_v2
  Colorectal = "[HG-U219] Affymetrix Human Genome U219 Array",
  Lung = "[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array", #affy_hg_u133_plus_2
  Prostate = "[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array" #affy_hg_u133_plus_2
)

data <- clean.data.without.entrez.id(data)
colnames(data) <- gene.id.translate(colnames(data), 'affy_hg_u133a')

expressionMatrix = expression.matrix.pearson(t(testFile))
biologicalMatrix = biological.matrix.GO(colnames(testFile))