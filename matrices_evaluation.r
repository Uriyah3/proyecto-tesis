source("data_sampling.r")
source("gpl_chip_to_entrez_id.r")
source("globals.r")
library("amap")
library(org.Hs.eg.db)
library(KEGGREST)
library(STRINGdb)

#' Calculate the expression matrix using the data in \code{gene_profiles} and store
#' the resulting data frame into cache. Loads the data from cache it it was
#' calculated before.
#' 
#' @param gene_profiles Matrix of gene/samples.
#' @param metric Metric used to calculate amap:Dist. Default: 'abscorrelation'
#' @param nbproc Number of processed to use in amap::Dist. Default: 2
#' @param dataset Name of the dataset being processed. Used to store the resulting
#' matrix into cache.
#' @return Data frame with expression matrix. This means, pairwise distance
#' between genes.
#'
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

#' Calculate the expression matrix using the data in \code{gene_profiles} and store
#' the resulting data frame into cache. Loads the data from cache it it was
#' calculated before.
#' 
#' @param gene_profiles Matrix of gene/samples.
#' @param metric Metric used to calculate amap:Dist. Default: 'abscorrelation'
#' @param nbproc Number of processed to use in amap::Dist. Default: 2
#' @param dataset Name of the dataset being processed. Used to store the resulting
#' matrix into cache.
#' @return Data frame with expression matrix. This means, pairwise distance
#' between genes.
#'
expression.matrix <- function(gene_profiles, metric = 'abscorrelation', nbproc = 2, dataset = '') {
  if( dataset != '' ) {
    datafile <- tools::file_path_sans_ext(basename(dataset))
    cached_filename <- paste("cache/", datafile, "-expression.rda", sep="")
    
    if( file.exists(cached_filename) ) {
      return( readRDS(cached_filename) )
    }
  }
  
  data <- as.data.frame( as.matrix(amap::Dist(gene_profiles, metric, nbproc = 2)) )
  if( dataset != '' ) {
    saveRDS(data, file = cached_filename)
  }
  
  return( data )
}

#' Check if the pair of dataset and biological database have been calculated before
#' (check in cache) and return it. Otherwise, calculate them and store into cache.
#' 
#' @param gene_list Vector of genes using their entrezgene_id.
#' @param database Name of the biological data source (check globals.r)
#' @param dataset Name of the dataset being processed
#' @return Data frame with biological matrix.
#'
biological.matrix <- function(gene_list, database, dataset = '') {
  # Create the unique filename of the cached data
  datafile <- tools::file_path_sans_ext(basename(dataset))
  cached_filename <- paste("cache/", datafile, "-", database, ".rda", sep="")
  
  # Check the cache
  if( file.exists(cached_filename) ) {
    dmatrix <- readRDS(cached_filename)
  }
  # If there isn't anything, create it and save it into the file cache
  else {
    if( database == biological_databases$go ) {
      dmatrix <- biological.matrix.go(gene_list)
    } else if( database == biological_databases$string ) {
      dmatrix <- biological.matrix.string(gene_list)
    } else if( database == biological_databases$kegg ) {
      dmatrix <- biological.matrix.kegg.pathway(gene_list)
    } else if( database == biological_databases$mesh ) {
      dmatrix <- biological.matrix.mesh(gene_list)
    } else if( database == biological_databases$do ) {
      dmatrix <- biological.matrix.disease.ontology(gene_list)
    } else {
      stop("Wrong database name given")
    }
    
    saveRDS(dmatrix, file = cached_filename)
  }
  
  return( dmatrix )
}

#' Fill genes that are missing from the \code{bmatrix} with a distance greater
#' than the maximum distance in bmatrix but no greater than 1.00. These values 
#' may be missing because the biological source used doesn't have information
#' about said missing genes.
#' 
#' @param gene_list Vector of genes using their entrezgene_id.
#' @param bmatrix Biological distance matrix (using the genes in gene_list).
#' @return Data frame with biological distance between all pair of genes in \code{gene_list}
#' where NA values are filled with a distance greater than max() but no more than 1.00
#' @note To avoid dividing by 0 in other places, all distances=0 are replaced by 
#' a distance that is 0.5 times smaller than the minimum distance in \code{bmatrix}
#'
biological.matrix.fill.missing <- function(gene_list, bmatrix) {
  fill_value <- min(max(bmatrix[!is.na(bmatrix)]) * 1.01, 1.00)
  
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

#' Calculate a biological distance matrix between the genes in \code{gene_list}
#' using Gene Ontology as the biological data source.
#' 
#' @param gene_list Vector of genes using their entrezgene_id.
#' @return Data frame with biological distance between all pair of genes in \code{gene_list}
#'
biological.matrix.go <- function(gene_list) {
  d <- GOSemSim::godata('org.Hs.eg.db', ont="MF", computeIC=FALSE)
  matrix <- GOSemSim::mgeneSim(gene_list, semData=d, measure="Wang")
  # Transform the similarity matrix into a distance matrix 1-Wang
  matrix <- as.matrix(as.dist(1-matrix))
  return( biological.matrix.fill.missing(gene_list, matrix) )
}

#' Download the KEGG database and store the gene data needed by the method
#' BioCor::mgeneSim to calculate similarities between genes in cache.
#' 
#' @return A list of lists where the names are the geneentrez_ids and the values
#' are the paths each gene participates in.
#' @note The cache created by this method will be overwritten when its more than
#' 10 days old. This is necessary since this database is constantly being updated.
#' 
helper.download.kegg <- function() {
  # Save the file to avoid too many requests to the API
  cached_time_in_days <- 10
  cached_filename <- paste("cache/","kegg-hsa-data-by-gene.rda", sep="")
  
  if( file.exists(cached_filename) ) {
    oldness <- as.numeric( difftime(Sys.Date(), file.info(cached_filename)$ctime, units=c('days')) )
    if( oldness < cached_time_in_days ) {
      return( as.list(readRDS(cached_filename)) )
    }
  }
  hsa_links <- keggLink("pathway", "hsa")
  
  names(hsa_links) <- sapply(names(hsa_links), function(hsa) gsub("hsa:", "", hsa))
  hsa_data_by_gene <- split(unlist(hsa_links, use.names = FALSE), rep(names(hsa_links), lengths(hsa_links)))
  
  saveRDS(hsa_data_by_gene, file = cached_filename)
  
  return(as.list(hsa_data_by_gene))
}

#' Calculate a biological distance matrix between the genes in \code{gene_list}
#' using KEGG as the biological data source.
#' 
#' @param gene_list Vector of genes using their entrezgene_id.
#' @return Data frame with biological distance between all pair of genes in \code{gene_list}
#'
biological.matrix.kegg.pathway <- function(gene_list) {
  d <- helper.download.kegg()
  matrix <- BioCor::mgeneSim(gene_list, d, method="BMA")
  # Transform the similarity matrix into a distance matrix 1-sim
  matrix <- as.matrix(as.dist(1-matrix))
  return( biological.matrix.fill.missing(gene_list, matrix) )
}

#' Download the string database and store its adjacency matrix into cache.
#' 
#' Note: This method follows what was done by Jonathan Ronen and Altuna Akalin in
#' https://bioconductor.org/packages/devel/bioc/vignettes/netSmooth/inst/doc/buildingPPIsFromStringDB.html
#' 
#' @param version The version of the database to download/use from stringdb
#' @param species The identifier of a species from string. Default: 9606 (human)
#' @param score_threshold Minimum combined_score acceptable. Default: 0
#' @param input_directory Directory where stringdb data is stored. Default: "stringdb/"
#' @return An adjacency (sparse) matrix with the combined scores between all pairs 
#' of genes that exist in the human stringdb. The entrezgene_id is used for
#' colnames and rownames.
#' 
helper.download.string <- function(version="11.0", species=9606, score_threshold = 0, input_directory = "stringdb/") {
  cached_filename <- paste("cache/","stringdb-score-data.rda", sep="")
  if( file.exists(cached_filename) ) {
    return( readRDS(cached_filename) )
  }
  
  string_db <- STRINGdb$new( version=version, species=species,
                             score_threshold=score_threshold, input_directory=input_directory)
  
  human_graph <- string_db$get_graph()
  adj_matrix <- as_adjacency_matrix(human_graph, attr="combined_score")
  
  mart=useMart(host = 'grch37.ensembl.org',
               biomart='ENSEMBL_MART_ENSEMBL',
               dataset='hsapiens_gene_ensembl')
  
  protein_ids <- sapply(strsplit(rownames(adj_matrix), '\\.'), function(x) x[2])
  
  mart_results <- getBM(attributes = c("entrezgene_id", "ensembl_peptide_id"),
                        filters = "ensembl_peptide_id", values = protein_ids,
                        mart = mart)
  
  ix <- match(protein_ids, mart_results$ensembl_peptide_id)
  ix <- ix[!is.na(ix)]
  
  newnames <- protein_ids
  newnames[match(mart_results[ix,'ensembl_peptide_id'], newnames)] <-
    mart_results[ix, 'entrezgene_id']
  rownames(adj_matrix) <- newnames
  colnames(adj_matrix) <- newnames
  
  ppi <- adj_matrix[!duplicated(newnames), !duplicated(newnames)]
  nullrows <- Matrix::rowSums(ppi)==0
  ppi <- ppi[!nullrows,!nullrows]
  
  saveRDS(ppi, file = cached_filename)
  
  return(ppi)
}

#' Calculate a biological distance matrix between the genes in \code{gene_list}
#' using STRING as the biological data source.
#' 
#' @param gene_list Vector of genes using their entrezgene_id.
#' @return Data frame with biological distance between all pair of genes in \code{gene_list}
#'
biological.matrix.string <- function(gene_list) {
  string_db <- helper.download.string()
  
  sim <- string_db[ rownames(string_db) %in% gene_list, colnames(string_db) %in% gene_list ]
  matrix <- as.matrix(sim)
  # Make the difference between no score and score>=1 a bit bigger.
  matrix[matrix == 0] <- - (max(matrix) / 100)
  matrix <- apply(matrix, 1, function(x)(x-min(x))/(max(x)-min(x)))
  matrix <- as.matrix(as.dist(1-matrix))
  # Skip filling values that have no data inside the method biological.matrix.fill.missing
  matrix[matrix == 1] <- NA
  matrix[is.na(matrix)] <- min(max(matrix[!is.na(matrix)]) * 1.12, 1.00)
  
  return( biological.matrix.fill.missing(gene_list, matrix) )
} 

#' Calculate a biological distance matrix between the genes in \code{gene_list}
#' using MeSH as the biological data source.
#' 
#' @param gene_list Vector of genes using their entrezgene_id.
#' @return Data frame with biological distance between all pair of genes in \code{gene_list}
#'
biological.matrix.mesh <- function(gene_list) {
  hsamd_cache_filename <- "cache/hsamd.rda"
  if( file.exists(hsamd_cache_filename) ) {
    hsamd <- readRDS(hsamd_cache_filename)
  } else {
    hsamd <- meshes::meshdata("MeSH.Hsa.eg.db", category='A', computeIC=T, database="gendoo")
    saveRDS(hsamd, file = hsamd_cache_filename)
  }
  
  meshes::geneSim(gene_list, gene_list, semData=hsamd, measure="Wang")
  # Transform the similarity matrix into a distance matrix 1-Wang
  matrix <- as.matrix(as.dist(1-matrix))
  return( biological.matrix.fill.missing(gene_list, matrix) )
}

#' Calculate a biological distance matrix between the genes in \code{gene_list}
#' using Disease Ontology as the biological data source.
#' 
#' @param gene_list Vector of genes using their entrezgene_id.
#' @return Data frame with biological distance between all pair of genes in \code{gene_list}
#'
biological.matrix.disease.ontology <- function(gene_list) {
  matrix <- DOSE::geneSim(gene_list, gene_list, measure="Wang")
  # Transform the similarity matrix into a distance matrix 1-Wang
  matrix <- as.matrix(as.dist(1-matrix))
  return( biological.matrix.fill.missing(gene_list, matrix) )
}
