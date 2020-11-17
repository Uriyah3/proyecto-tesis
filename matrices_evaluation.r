source("data_sampling.R")
library("amap")
library("biomaRt")
library(org.Hs.eg.db)

if( !exists("human_mart") ) {
  human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
}
gene.id.translate <- function(gene_list, chip = 'affy_hg_u133_plus_2') {
  # Sanity check
  if( substring(gene_list[1], 1, 1) == "X" ) {
    gene_list = sapply(gene_list, function(gene) return(substring(gene, 2)), USE.NAMES = FALSE)
  }
  
  translated_gene_list <- getBM(attributes = c(chip, 'entrezgene_id'),
                                 filters = chip,
                                 values = gene_list, 
                                 mart = human_mart)
  print(translated_gene_list)
  
  translated_gene_list <- translated_gene_list[!duplicated(translated_gene_list[,1]),]
  return( translated_gene_list$entrezgene_id )
}

expression.matrix <- function(gene_profiles, metric = 'abscorrelation', nbproc = 2) {
  data <- as.matrix(amap::Dist(gene_profiles, metric, nbproc = 2))
  return( data )
}

biological.matrix.go <- function(gene_list) {
  d <- GOSemSim::godata('org.Hs.eg.db', ont="MF", computeIC=FALSE)
  matrix <- GOSemSim::mgeneSim(gene_list, semData=d, measure="Wang")
  return( matrix )
}

biological.matrix.kegg.pathways <- function(gene_list) {
  
}

biological.matrix.string <- function(gene_list) {
  
}

biologiacal.matrix.mesh <- function(gene_list) {
  
}
