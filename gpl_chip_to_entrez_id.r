#' Reads a GPL data file using read.csv
#' Note: The GPL files were downloaded from CuMiDa
#' 
#' @param filename The path to the file
#' @return A matrix with the information stored in \code{filename}
#' @examples
#' gpl_data = read.gpl("data/chip_to_entrez_id/GPL96-57554.txt.gz")
#' 
read.gpl <- function(filename) {
  gpl = read.csv(filename, header = TRUE, sep = "\t", quote= "", dec = ".", row.names = 1, check.names = FALSE)
}

# Convert the respective gpl tool to the file that stores its data
gpl_to_file <- list(
  GPL96 = "data/chip_to_entrez_id/GPL96_limpo.txt.gz",
  GPL3921 = "data/chip_to_entrez_id/GPL3921_limpo.txt.gz",
  GPL570 = "data/chip_to_entrez_id/GPL570_limpo.txt.gz",
  GPL13607 = "data/chip_to_entrez_id/GPL13607_limpo.txt.gz",
  GPL13667 = "data/chip_to_entrez_id/GPL13667_limpo.txt.gz"
)

#' Returns the ENTREZ_GENE_ID for all genes in the GPL tool
#' 
#' @param chip Code of the GPL chip used in the experiment
#' @return A matrix with gene to ENTREZ_GENE_ID info for the GPL \code{chip}
#' @examples
#' gene_to_id = process.gpl("GPL96")
#'
process.gpl <- function(chip) {
  gpl_data = read.gpl(gpl_to_file[[chip]])
  gpl_data = gpl_data[ ,"ENTREZ_GENE_ID", drop=FALSE]
}

