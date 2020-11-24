source("data_sampling.r")
source("gpl_chip_to_entrez_id.r")
library("amap")
library(org.Hs.eg.db)

clean.and.translate.entrez.id <- function(data, chip) {
  # Sanity check if R added an X to the column name
  gene_list = colnames(data)
  if( substring(gene_list[1], 1, 1) == "X" ) {
    gene_list = sapply(gene_list, function(gene) return(substring(gene, 2)), USE.NAMES = FALSE)
    
    colnames(data) <- gene_list
  }
  
  gene_translator <- process.gpl(chip)
  genes_to_keep <- gene_translator[gene_translator$ENTREZ_GENE_ID != "", , drop=FALSE]
  
  # Translate to ENTREZ_GENE_ID
  data <- data[, colnames(data) %in% rownames(genes_to_keep), drop=FALSE]
  colnames(data) <- sapply(colnames(data), function(gene_name) {
    unlist(strsplit(gene_translator[gene_name, ], " /// "))[1]
  })
  # remove duplicates
  data <- data[, !duplicated(colnames(data))]
  return( data )
}

expression.matrix <- function(gene_profiles, metric = 'abscorrelation', nbproc = 2) {
  data <- as.matrix(amap::Dist(gene_profiles, metric, nbproc = 2))
  return( data )
}

biological.matrix.fill.missing <- function(gene_list, bmatrix) {
  fill_value <- min(max(bmatrix) * 1.01, 1.00)
  
  bmatrix <- as.data.frame(bmatrix)
  bmatrix[gene_list[!(gene_list %in% rownames(bmatrix))], ] <- NA
  bmatrix[ , gene_list[!(gene_list %in% colnames(bmatrix))]] <- NA
  bmatrix[is.na(bmatrix)] <- fill_value
  
  # Temporal fix for distance = 0 (on non diagonals)
  bmatrix[ bmatrix == 0 ] <- min( bmatrix[bmatrix > 0.001] ) * 0.5
  # Set diagonal to be 0
  bmatrix[ row(bmatrix) == col(bmatrix) ] <- 0
  return( bmatrix )
}

biological.matrix.go <- function(gene_list) {
  d <- GOSemSim::godata('org.Hs.eg.db', ont="MF", computeIC=FALSE)
  matrix <- GOSemSim::mgeneSim(gene_list, semData=d, measure="Wang")
  # Transform the similarity matrix into a distance matrix 1-Wang
  matrix <- as.matrix(as.dist(1-matrix))
  return( biological.matrix.fill.missing(gene_list, matrix) )
}

biological.matrix.kegg.pathways <- function(gene_list) {
  
}

biological.matrix.string <- function(gene_list) {
  
}

biologiacal.matrix.mesh <- function(gene_list) {
  
}
