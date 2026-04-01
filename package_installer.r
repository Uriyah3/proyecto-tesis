options(repos = c(CRAN = "https://cran.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

list.of.packages <- c("ggplot2", "hash", "nsga2R", "amap", "lobstr", "profvis",
        "stringr", "tools", "mco", "cluster", "utils", "reshape2", "irace", "future",
        "igraph", "dplyr", "optparse", "future.apply", "inline", "magick",
        "httr", "xml2")
list.of.biocmanager.packages <- c("GOSemSim", "DOSE",
        "org.Hs.eg.db", "KEGGREST", "BioCor", "STRINGdb", "qusage",
        "biomaRt", "ComplexHeatmap")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

new.biocmanager.packages <- list.of.biocmanager.packages[!(list.of.biocmanager.packages %in% installed.packages()[,"Package"])]
if(length(new.biocmanager.packages)) BiocManager::install(new.biocmanager.packages)
