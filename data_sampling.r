library(stringr)
library(tools)
source("file_utils.r")

#' Create sampled data
#' Creates 3 sampled versions of the datasets contained in data/training:
#'   (row_sample) 10% or row
#'   (col_sample) 10% of cols
#'   (all_sample) 10% of rows and cols
#' Can create various samples from the same dataset
#' Stores the generated files in data/training/samples
#' 
#' @param samples_per_file Amount of sampled dataset to generate from each 
#' training dataset
#' @param sampling_percent Percentaje to sample from rows and columns
#' @param minimum_row_size The minimum amount of rows that can be left after sampling
#' @param minimum_col_size The minimum amount of columns that can be left after sampling
#' @param data_dir String with path to the data to be sampled
#' @return nothing
run.sampling <- function(samples_per_file = 10, sampling_percent = 5, minimum_row_size = 50, minimum_col_size = 500, data_dir="data/training/") {
  
  data_files <- list.files(data_dir)
  data_files <- paste0(data_dir, data_files)
  for (data_file in data_files) {
    if(file_ext(data_file) == 'gz' || file_ext(data_file) == 'csv') {
      message(str_interp("Generating samples for: ${data_file} ..."))
      generate.samples(data_file, samples_per_file, sampling_percent, minimum_row_size, minimum_col_size)
      gc()
    }
  }
  return()
}

#' Generates multiples samples using the data from \code{dataset_file}
#' @param dataset_file Path to the file to be sampled.
generate.samples <- function(dataset_file, number_of_samples, sampling_percent, minimum_row_size, minimum_col_size, sample_dir = "data/training/samples") {
  data <- read.dataset(dataset_file)
  output_file <- tools::file_path_sans_ext(basename(dataset_file))
  output_file <- paste(sample_dir, "/", output_file, sep = "")

  for (i in 1:number_of_samples) {
    message(str_interp("Generating sample ${i}/${number_of_samples}"))
    rows <- nrow(data)
    cols <- ncol(data)
    
    bounded_rows <- round(rows * (sampling_percent / 100))
    bounded_rows <- ifelse(minimum_row_size > bounded_rows, minimum_row_size, bounded_rows)
    bounded_rows <- ifelse(bounded_rows > rows, rows, bounded_rows)
    bounded_cols <- round(cols * (sampling_percent / 100))
    bounded_cols <- ifelse(minimum_col_size > bounded_cols, minimum_col_size, bounded_cols)
    bounded_cols <- ifelse(bounded_cols > cols, cols, bounded_cols)
    
    generate.plus.save.sample(data, rows, bounded_cols, str_interp("${output_file}-${i}-col-sample-${sampling_percent}-percent.csv"))
    generate.plus.save.sample(data, bounded_rows, cols, str_interp("${output_file}-${i}-row-sample-${sampling_percent}-percent.csv"))
    generate.plus.save.sample(data, bounded_rows, bounded_cols, str_interp("${output_file}-${i}-all-sample-${sampling_percent}-percent.csv"))
  }
  return()
}

#' Generate and save a single bounded sample
#' @param data Data to be sampled
#' @param rows Amount of rows to sample from \code{data}
#' @param cols Amount of columns to sample from \code{data}
#' @param output_file Where the sampled data is stored
generate.plus.save.sample <- function(data, rows, cols, output_file) {
  sampled_data <- data[sample(nrow(data), rows), sample(ncol(data), cols)]
  write.csv(sampled_data, output_file, row.names = TRUE, quote = FALSE)
  return()
}

