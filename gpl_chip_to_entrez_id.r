
read.gpl <- function(filename) {
  gpl = read.csv(filename, header = TRUE, sep = "\t", quote= "", dec = ".", row.names = 1)
}

GPL_CHIPS <- list(
  gpl96 = "data/chip_to_entrez_id/GPL96-57554.txt"
  
)

process.gpl <- function(chip) {
  gpl_data = read.gpl(GPL_CHIPS[[chip]])
  gpl_data = gpl_data[ ,"ENTREZ_GENE_ID", drop=FALSE]
}

