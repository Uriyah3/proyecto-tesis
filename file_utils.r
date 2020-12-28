library(utils)

#' Load a comma separated dataset file
#' 
#' @param data_file String with relative path to the dataset
#' @return data frame with data from the file \code{data_file}
read.dataset <- function(data_file) {
  read.csv(data_file, header = TRUE, sep = ",", quote= "", dec = ".", row.names = 1, check.names = FALSE)
}
