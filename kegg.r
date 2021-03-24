library(KEGGREST)
library(neo4jshell)
library(stringr)
library(igraph)

# This part of the code is now abandoned, while it was fun to try and use Neo4j
# to approach this problem, a more solid and robust alternative was found inside
# the BioCor package.
# NOTE: Some of the final methods are incomplete or haven't been tested, i.e.
# the KEGGSim and KEGGSim.matrix methods specifically.

neo4j.query <- function(query_string) {
  shell_path = "C:\\Users\\BlackOnyx\\AppData\\Local\\Neo4j\\Relate\\Data\\dbmss\\dbms-61d09d6c-dcf1-47d2-a2e0-38f53e229a21\\bin\\cypher-shell.bat"
  neo_conn_kegg <- list(address = "bolt://localhost:7687", uid = "neo4j", pwd = "KEGG")
  neo4j_query(con = neo_conn_kegg, qry = query_string, shell_path = shell_path)
}

neo4j.create.graph <- function() {
  create_graph_query <- "CALL gds.graph.create(
  'KEGG_Graph',
  ['Gene', 'Pathway'],
  'PARTICIPATES_IN'
);"
  neo4j.query(create_graph_query)
}

neo4j.calculate.similarities <- function() {

  calculate_similarities <- "CALL gds.nodeSimilarity.write('KEGG_Graph', {
    writeRelationshipType: 'SIMILAR',
    writeProperty: 'score'
})
YIELD nodesCompared, relationshipsWritten, createMillis, computeMillis, postProcessingMillis, similarityDistribution, configuration"
  neo4j.query(calculate_similarities);
}

# Esta función utiliza la Neo4j Graph Data Science (GDS) library.
kegg.build.neo4j.database <- function() {
  # Download link/pathway data for human genome and insert to neo4j
  hsa_links <- keggLink("pathway", "hsa")
  hsa_paths <- keggList('pathway', 'hsa')
  hsa_genes <- keggList('hsa')
  
  # Group genes by pathway
  # hsa_data <- split(rep(names(hsa_links), lengths(hsa_links)), unlist(hsa_links))
  
  message("Creating Pathways (hsa_paths)...")
  insert <- ""
  group_by = 50
  for(i in 1:length(hsa_paths)) {
    path_id <- attr(hsa_paths[i], 'name')
    path_name <- hsa_paths[[i]]
    insert <- paste( insert, str_interp("CREATE (p${i}:Pathway {id: \"${path_id}\", name: \"${path_name}\"})"), sep="\n" )
    if(i %% group_by == 0 || i == length(hsa_paths)) {
      neo4j.query(paste(insert, ';', sep='')); 
      insert <- ""
    }
  }
  
  message("Creating Genes (hsa_genes)...")
  insert <- ""
  group_by = 40
  for(i in 1:length(hsa_genes)) {
    gene_id <- attr(hsa_genes[i], 'name')
    gene_name <- hsa_genes[[i]]
    insert <- paste( insert, str_interp("CREATE (g${i}:Gene {id: \"${gene_id}\", name: \"${gene_name}\"})"), sep="\n" )
    if(i %% group_by == 0 || i == length(hsa_genes)) {
      neo4j.query(paste(insert, ';', sep='')); 
      insert <- ""
    }
  }
  
  message("Creating Links (hsa_links)...")
  insert <- ""
  group_by = 50
  for(i in 30431:length(hsa_links)) {
    gene_id <- attr(hsa_links[i], 'name')
    path_id <- hsa_links[[i]]
    #insert <- paste( insert, str_interp("CREATE (`${gene_id}`)-[:PARTICIPATES_IN]->(`${path_id}`)"), collapse="\n")
    insert <- paste( insert, str_interp("MATCH (g:Gene), (p:Pathway) WHERE g.id=\"${gene_id}\" AND p.id=\"${path_id}\" CREATE (g)-[:PARTICIPATES_IN]->(p)"), sep="\n") 
    if(i %% group_by == 0 || i == length(hsa_links)) {
      neo4j.query(paste(insert, ';', sep='')); 
      insert <- ""
    } else {
      # Este comando se agrega para poder enviar varios MATCH en una misma query
      insert <- paste(insert, "\nWITH count(*) as dummy", sep="")
    }
  }
  
  neo4j.create.graph()
  message("Calculating similarity between genes...")
  results <- neo4j.calculate.similarities()
}

kegg.build.sim.matrix <- function() {
  qry <- "MATCH (g1:Gene)-[r:SIMILAR]->(g2:Gene) RETURN g1.id, g2.id, r.score;"
  sims <- neo4j.query(qry)
  
  G <- graph.data.frame(sims,directed=TRUE);
  sims_matrix <- as_adjacency_matrix(G,names=TRUE,sparse=FALSE,attr="r.score");
}

KEGGSim <- function(gene1, gene2) {
  kegg_datafile <- 'kegg-sim'
  cached_filename <- paste("cache/", kegg_datafile, ".rda", sep="")
  
  if( file.exists(cached_filename) ) {
    kegg_sims <- readRDS(cached_filename)
  } else {
    return(NULL)
  }
  
  return( kegg_sims[gene1, gene2] )
}

KEGGSim.matrix <- function(gene_list) {
  kegg_datafile <- 'kegg-sim'
  cached_filename <- paste("cache/", kegg_datafile, ".rda", sep="")
  
  if( file.exists(cached_filename) ) {
    kegg_sims <- readRDS(cached_filename)
  } else {
    kegg.build.neo4j.database()
    kegg_sims <- kegg.build.sim.matrix()
    saveRDS(kegg_sims, file = cached_filename)
  }
  
}
