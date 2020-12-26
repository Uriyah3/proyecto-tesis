library(KEGGREST)
library(neo4jshell)
library(stringr)

neo4j.query(query_string) <- function {
  shell_path = "C:\\Users\\BlackOnyx\\AppData\\Local\\Neo4j\\Relate\\Data\\dbmss\\dbms-61d09d6c-dcf1-47d2-a2e0-38f53e229a21\\bin\\cypher-shell"
  neo_conn_kegg <- list(address = "bolt://localhost:7687", uid = "neo4j", pwd = "KEGG")
  neo4j_query(con = neo_conn_kegg, qry = query_string, shell_path = shell_path)
}

kegg.build.neo4j.database <- function() {
  # Download link/pathway data for human genome and insert to neo4j
  hsa_links <- keggLink("pathway", "hsa")
  hsa_paths <- keggList('pathway', 'hsa')
  hsa_genes <- keggList('hsa')
  
  # Group of cliques that belong to a pathway
  # hsa_data <- split(rep(names(hsa_links), lengths(hsa_links)), unlist(hsa_links))
  
  create_hsa_paths <- sapply(1:length(hsa_paths), function(i) {
    path_id <- attr(hsa_paths[i], 'name')
    path_name <- hsa_paths[[i]]
    str_interp("CREATE (`${path_id}`:Pathway {name: '${path_name}'})")
  })
  create_hsa_paths <- paste( c(paste(create_hsa_paths, collapse="\n"), ';'), collapse='')
  
  for(i in 1:length(hsa_paths)) {
    path_id <- attr(hsa_paths[i], 'name')
    path_name <- hsa_paths[[i]]
    insert <- str_interp("CREATE (`${path_id}`:Pathway {name: '${path_name}'})")
  }
  
  for(i in 1:length(hsa_genes)) {
    gene_id <- attr(hsa_genes[i], 'name')
    gene_name <- hsa_genes[[i]]
    insert <- str_interp("CREATE (`${gene_id}`:Gene {name: '${gene_name}'})")
  }
  
  for(i in 1:length(hsa_links)) {
    gene_id <- attr(hsa_links[i], 'name')
    path_id <- hsa_links[[i]]
    insert <- str_interp("CREATE (`${gene_id}`)-[:PARTICIPATES_IN]->(`${path_id}`)")
  }
 
  calculate_similarities <- "CALL algo.nodeSimilarity('Gene | Pathway', 'PARTICIPATES_IN', {
  direction: 'OUTGOING',
  write: true
})
YIELD nodesCompared, relationships, write, writeProperty, writeRelationshipType, loadMillis, computeMillis, writeMillis, postProcessingMillis, min, max, mean, stdDev, p1, p5, p10, p25, p50, p75, p90, p95, p99, p100;"
  results <- neo4j.query(calculate_similarities);
  print(results)
}