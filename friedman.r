# https://gist.github.com/vv111y/aa6a67fb080814687a3acfe327fca614 
friedman.test.with.post.hoc <- function(data, alpha = 0.05, metrica.nombre = 'Silueta máxima')
{
  print("Check if you missing the packages 'graph' and 'Rgraphviz'. Try to install them using bioconductor")
  
  #source("http://bioconductor.org/biocLite.R")
  #biocLite(c("graph","Rgraphviz"))
  
  # Loading needed packages
  if(!require(ggplot2))
  {
    print("You are missing the package 'ggplot2', we will now try to install it...")
    install.packages("ggplot2")
    library(ggplot2)
  }
  
  if(!require(scmamp))
  {
    print("You are missing the package 'scmamp', we will now try to install it...")
    devtools::install_github("b0rxa/scmamp")
    library(scmamp)
  }
  
  if(!require(pgirmess))
  {
    print("You are missing the package 'pgirmess', we will now try to install it...")
    install.packages("pgirmess")
    library(pgirmess)
  }
  
  pre.results <- friedmanTest(data)
  imanDavenport.result <- imanDavenportTest(data)
  
  if(pre.results$p.value < alpha){
    post.results <- NULL
    
    if(length(colnames(data)) > 9){
      post.results <- postHocTest(data=data, test="friedman", correct="shaffer")
    }else{
      post.results <- postHocTest(data=data, test="friedman", correct="bergmann")
    }
    
    ## LaTeX formated: Significances highlighted in bold
    avg.val <- post.results$summary
    
    if(metrica.nombre == 'DAVID #grupos') {
      best.res <- avg.val == min(avg.val)
    } else {
      best.res <- avg.val == max(avg.val)
    }
    stat.diff <- post.results$corrected.pval < alpha
    stat.diff[is.na(stat.diff)] <- FALSE
    stat.diff <- colSums(stat.diff) == 4
    
    writeTabular(table = avg.val, format = 'f', bold = best.res, mark = stat.diff, digits = 3)
    
    bold <- post.results$corrected.pval < alpha
    bold[is.na(bold)] <- FALSE
    writeTabular(table=post.results$corrected.pval, format='f', bold=bold, hrule=0, vrule=0, digits = 5)
    
    friedman.mc <- friedmanmc(data.matrix(data))
    
    plt <- plotPvalues(post.results$corrected.pval, alg.order=order(post.results$summary)) + labs(title=paste("p-values corregidos utilizando procedimiento \nde Bergmann y Hommel para ", metrica.nombre,sep="")) + xlab("Metaheurística") + ylab("Metaheurística") + scale_fill_gradientn("p-values corregidos" , colours = c("grey15" , "grey30")) + theme(plot.title = element_text(hjust = 0.5))
    file.metric <- name.mapper[[metrica.nombre]]
    ggsave(str_interp("friedman-bergmann-${file.metric}.png"), device="png", path="plots", width=6.62, height = 3.47, dpi=100) # in pixels width=662, height=347
    list.to.return <- list(Friedman = pre.results, ImanDavenport = imanDavenport.result, PostHoc = post.results, FriedmanMC = friedman.mc, Plt = plt)
    print(list.to.return$Friedman)
    return(list.to.return)
  }
  else{
    print("The results where not significant. There is no need for a post-hoc test.")
    list.to.return <- list(Friedman = pre.results, ImanDavenport = imanDavenport.result, PostHoc = NULL, FriedmanMC = NULL, Plt = NULL)
    return(list.to.return)
  }
}

name.mapper <- list(
  "Silueta máxima" = 'silueta-maxima',
  "Silueta promedio" = 'silueta-promedio',
  "DAVID #grupos" = 'david-grupos',
  "DAVID enrichment máximo" = 'david-enrichment-maximo',
  "DAVID enrichment promedio" = 'david-enrichment-promedio'
)


statistic.test <- function() {
  ev.datasets <- list('GSE89116', 'GSE53757', 'GSE31189', 'GSE50161', 'GSE6919_U95Av2')
  metaheuristics <-  c("go", "string", "kegg", "disgenet_dis", "base")
  bio.columns <-  c("GO", "STRING", "KEGG", "DisGeNET", "Base")
  metrics <- list('mean.silhouette', 'max.silhouette', 'david.groups', 'david.max.enrichment', 'david.mean.enrichment')
  
  mean.silhouette <- sapply(metaheuristics, function(metaheuristic) {
    sapply(ev.datasets, function(dataset.key) {
      dataset <- datasets[[dataset.key]]
      sapply(1:13, function(i) {
        silhouette_results <- load.evaluation.from.cache(dataset$name, metaheuristic, i, 'silhouette')
        silhouette_results$mean_silhouette
      })
    })
  })
  
  max.silhouette <- sapply(metaheuristics, function(metaheuristic) {
    sapply(ev.datasets, function(dataset.key) {
      dataset <- datasets[[dataset.key]]
      sapply(1:13, function(i) {
        silhouette_results <- load.evaluation.from.cache(dataset$name, metaheuristic, i, 'silhouette')
        silhouette_results$max_silhouette
      })
    })
  })
  
  colnames(mean.silhouette) <- c("GO", "STRING", "KEGG", "DisGeNET", "Base")
  colnames(max.silhouette) <- c("GO", "STRING", "KEGG", "DisGeNET", "Base")
  # load certain results from all datasets (silhouette avg y max, david 3 metrics)
  # apply statistic test considering each iteration result with same label
}
# 
# max.silhouette <- structure(c(0.692, 0.639, 0.572, 0.593, 0.569,
#                                0.66, 0.619, 0.519, 0.535, 0.478,
#                                0.694, 0.64, 0.588, 0.557, 0.573,
#                                0.696, 0.633, 0.545, 0.585, 0.505,
#                                0.494, 0.424, 0.304, 0.259, 0.371), .Dim = c(5L, 5L), 
#                             .Dimnames = list(c("Mama GSE89116", "Riñón GSE53757", "Vejiga GSE31189", "Cerebro GSE50161", "Próstata GSE6919_U95Av2"), c("GO", "STRING", "KEGG", "DisGeNET", "Base")))
# 
# mean.silhouette <- structure(c(0.541, 0.508, 0.422, 0.357, 0.432,
#                                0.519, 0.399, 0.344, 0.332, 0.34,
#                                0.553, 0.473, 0.476, 0.381, 0.441,
#                                0.56, 0.439, 0.372, 0.327, 0.379,
#                                0.327, 0.284, 0.207, 0.184, 0.242), .Dim = c(5L, 5L), 
#                             .Dimnames = list(c("Mama GSE89116", "Riñón GSE53757", "Vejiga GSE31189", "Cerebro GSE50161", "Próstata GSE6919_U95Av2"), c("GO", "STRING", "KEGG", "DisGeNET", "Base")))

david.groups <- structure(c(6219, 5971, 5939, 4237, 6428, 
                            6558, 6695, 4665, 5754, 5676, 
                            5590, 4157, 6385, 6421, 6150, 
                            4077, 6430, 5866, 6152, 3964), .Dim = c(4L, 5L), 
                             .Dimnames = list(c("Riñón GSE53757", "Vejiga GSE31189", "Cerebro GSE50161", "Próstata GSE6919_U95Av2"), c("GO", "STRING", "KEGG", "DisGeNET", "Base")))

david.max.enrichment <- structure(c(97.24, 121.559, 88.23, 48.157, 43.512, 
                                    163.725, 76.481, 67.359, 150.392, 206.956, 
                                    169.315, 153.918, 58.471, 108.748, 94.667, 
                                    77.844, 68.25, 170.936, 82.301, 164.975), .Dim = c(4L, 5L), 
                          .Dimnames = list(c("Riñón GSE53757", "Vejiga GSE31189", "Cerebro GSE50161", "Próstata GSE6919_U95Av2"), c("GO", "STRING", "KEGG", "DisGeNET", "Base")))

david.mean.enrichment <- structure(c(0.914, 1, 0.843, 1.27, 0.715, 
                                     0.856, 0.718, 1.196, 2.337, 2.161, 
                                     2.56, 3.014, 0.769, 0.831, 0.874, 
                                     1.418, 0.903, 1.011, 0.869, 1.602), .Dim = c(4L, 5L), 
                                  .Dimnames = list(c("Riñón GSE53757", "Vejiga GSE31189", "Cerebro GSE50161", "Próstata GSE6919_U95Av2"), c("GO", "STRING", "KEGG", "DisGeNET", "Base")))

overall.ranking <- t(structure(c(2, 4, 1.6, 2.4, 5, 
                                2, 3.8, 1.4, 2.8, 5, 
                                2.75, 4.75, 1.5, 3, 3, 
                                3.5, 4.25, 1.25, 3.5, 2.5,
                                3.25, 4.75, 1, 3.5, 2.5), .Dim = c(5L, 5L), 
                              .Dimnames = list(c("GO", "STRING", "KEGG", "DisGeNET", "Base"), c("Max sil", "Mean sil", "groups", "Max enrichment", "Mean enrichment"))))


all.mean.ranking <- t(structure(c(2, 4, 1, 3, 5, 
                                 2, 4, 1, 3, 5, 
                                 2.75, 4.75, 1.5, 3, 3, 
                                 3.5, 4.25, 1.25, 3.5, 2.5,
                                 3.25, 4.75, 1, 3.5, 2.5), .Dim = c(5L, 5L), 
                               .Dimnames = list(c("GO", "STRING", "KEGG", "DisGeNET", "Base"), c("Max sil", "Mean sil", "groups", "Max enrichment", "Mean enrichment"))))

friedman.test.with.post.hoc(max.silhouette)
friedman.test.with.post.hoc(mean.silhouette, metrica.nombre = 'Silueta promedio')
friedman.test.with.post.hoc(david.groups, metrica.nombre = 'DAVID #grupos', alpha=0.1)
friedman.test.with.post.hoc(david.max.enrichment, metrica.nombre = 'DAVID enrichment máximo', alpha=0.1)
friedman.test.with.post.hoc(david.mean.enrichment, metrica.nombre = 'DAVID enrichment promedio', alpha=0.1)