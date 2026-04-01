#' DAVID Web Service Client
#'
#' A lightweight R client for the DAVID Bioinformatics web service API.
#' Replaces the abandoned RDAVIDWebService Bioconductor package (which required
#' Java/rJava) with simple HTTP GET requests via httr.
#'
#' DAVID is a stateful web service -- all calls within a workflow must share
#' the same HTTP session (cookies). This client uses an httr handle as a
#' cookie jar to maintain session state.
#'
#' @note DAVID moved from david.ncifcrf.gov to davidbioinformatics.nih.gov
#'   in March 2024.

library(httr)
library(xml2)

# Base URL for the DAVID web service API.
DAVID_BASE_URL <- "https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService"

# Session handle (cookie jar) -- shared across all calls in a workflow.
david_session <<- NULL

# ---------------------------------------------------------------------------
# Session management
# ---------------------------------------------------------------------------

#' Start a new DAVID session (resets cookies)
#'
#' Creates a fresh httr handle bound to the DAVID base URL so that all
#' subsequent requests share the same cookie jar.
#'
#' @return Invisibly returns the new handle.
david.new.session <- function() {
  david_session <<- handle(DAVID_BASE_URL)
}

#' Make a GET request to a DAVID endpoint within the current session
#'
#' Low-level helper that constructs the full URL, issues the request with
#' the shared session handle, and parses the XML response.
#'
#' @param operation Character string naming the DAVID SOAP operation
#'   (e.g. \code{"authenticate"}, \code{"addList"}).
#' @param params Named list of query parameters to pass.
#' @param timeout_secs Numeric timeout in seconds (default 60).
#' @return Parsed \code{xml2::xml_document}, or \code{NULL} on error.
david.request <- function(operation, params = list(), timeout_secs = 60) {
  if (is.null(david_session)) david.new.session()
  url <- paste0(DAVID_BASE_URL, "/", operation)

  # Use POST for operations with large payloads (addList with many genes)
  # to avoid HTTP 414 URI Too Long errors
  total_param_length <- sum(nchar(as.character(params)))
  if (total_param_length > 2000) {
    response <- POST(url, body = params, encode = "form", handle = david_session, timeout(timeout_secs))
  } else {
    response <- GET(url, query = params, handle = david_session, timeout(timeout_secs))
  }

  if (http_error(response)) {
    warning(paste("DAVID", operation, "HTTP error:", status_code(response)))
    return(NULL)
  }
  xml <- tryCatch(
    read_xml(content(response, as = "text", encoding = "UTF-8")),
    error = function(e) { warning(paste("DAVID XML parse error:", e$message)); NULL }
  )
  return(xml)
}

#' Extract the text content of the <return> element from DAVID XML
#'
#' Most DAVID SOAP responses wrap their payload in a \code{<return>} element.
#' This helper extracts the text content of the first such element using a
#' local-name match (namespace-agnostic).
#'
#' @param xml An \code{xml2::xml_document} or \code{NULL}.
#' @return Character string with the return value, or \code{NA} if not found.
david.extract.return <- function(xml) {
  if (is.null(xml)) return(NA)
  xml_text(xml_find_first(xml, ".//*[local-name()='return']"))
}

# ---------------------------------------------------------------------------
# API operations
# ---------------------------------------------------------------------------

#' Authenticate with the DAVID web service
#'
#' @param email Character string with a registered DAVID account email.
#' @return Logical \code{TRUE} if authentication succeeded, \code{FALSE}
#'   otherwise.
david.authenticate <- function(email) {
  xml <- david.request("authenticate", list(args0 = email), timeout_secs = 30)
  result <- david.extract.return(xml)
  return(!is.na(result) && grepl("true", result, ignore.case = TRUE))
}

#' Upload a gene list to the current DAVID session
#'
#' @param gene_list Character or numeric vector of gene identifiers.
#' @param id_type Character string specifying the identifier type
#'   (default \code{"ENTREZ_GENE_ID"}).
#' @param list_name Character string label for the list in DAVID.
#' @param list_type Integer flag: 0 for gene list, 1 for background
#'   (default 0).
#' @return Numeric count of successfully mapped genes, or \code{NA} on
#'   failure.
david.add.list <- function(gene_list, id_type = "ENTREZ_GENE_ID",
                           list_name = "gene_list", list_type = 0L) {
  ids_string <- paste(gene_list, collapse = ",")
  xml <- david.request("addList", list(
    args0 = ids_string, args1 = id_type, args2 = list_name, args3 = as.integer(list_type)
  ), timeout_secs = 60)
  result <- david.extract.return(xml)
  return(as.numeric(result))
}

#' Set the annotation categories for subsequent queries
#'
#' @param categories Character vector of DAVID annotation category names
#'   (see \code{DAVID_ANNOTATION_CATEGORIES} for the full list).
#' @return Logical \code{TRUE} if the request succeeded (non-NA response).
david.set.categories <- function(categories) {
  categories_string <- paste(categories, collapse = ",")
  xml <- david.request("setCategories", list(args0 = categories_string), timeout_secs = 30)
  result <- david.extract.return(xml)
  return(!is.na(result))
}

#' Retrieve the term cluster report from DAVID
#'
#' Runs DAVID's functional annotation clustering on the currently uploaded
#' gene list and selected categories.
#'
#' @param overlap Integer minimum gene overlap (default 4).
#' @param initial_seed Integer initial group membership (default 4).
#' @param final_seed Integer final group membership (default 4).
#' @param linkage Numeric linkage threshold (default 0.5).
#' @param kappa Integer kappa similarity threshold as percentage (default 35).
#' @return A list of clusters, each with \code{name} and
#'   \code{EnrichmentScore} elements. Returns an empty list if no clusters
#'   are found.
david.get.term.cluster.report <- function(overlap = 4L, initial_seed = 4L,
                                          final_seed = 4L, linkage = 0.5,
                                          kappa = 35L) {
  xml <- david.request("getTermClusterReport", list(
    args0 = as.integer(overlap), args1 = as.integer(initial_seed),
    args2 = as.integer(final_seed), args3 = as.numeric(linkage),
    args4 = as.integer(kappa)
  ), timeout_secs = 300)
  if (is.null(xml)) return(list())
  cluster_nodes <- xml_find_all(xml, ".//*[local-name()='return']")
  if (length(cluster_nodes) == 0) return(list())
  clusters <- lapply(cluster_nodes, function(node) {
    score_node <- xml_find_first(node, ".//*[local-name()='score']")
    name_node <- xml_find_first(node, ".//*[local-name()='name']")
    list(
      name = if (!is.na(name_node)) xml_text(name_node) else NA,
      EnrichmentScore = if (!is.na(score_node)) as.numeric(xml_text(score_node)) else 0
    )
  })
  return(clusters)
}

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#' All available DAVID annotation categories
DAVID_ANNOTATION_CATEGORIES <- c(
  "ENTREZ_GENE_ID", "BIOCARTA", "BBID", "BIOGRID_INTERACTION",
  "CGAP_EST_QUARTILE", "CGAP_SAGE_QUARTILE", "CHROMOSOME",
  "ENSEMBL_GENE_ID", "ENTREZ_GENE_SUMMARY", "GAD_DISEASE",
  "GAD_DISEASE_CLASS", "GENERIF_SUMMARY", "GNF_U133A_QUARTILE",
  "GOTERM_BP_ALL", "GOTERM_BP_DIRECT", "GOTERM_CC_ALL", "GOTERM_CC_DIRECT",
  "GOTERM_MF_ALL", "GOTERM_MF_DIRECT", "HIV_INTERACTION",
  "HIV_INTERACTION_CATEGORY", "HIV_INTERACTION_PUBMED_ID", "KEGG_PATHWAY",
  "MINT", "OMIM_DISEASE", "PFAM", "PIR_SEQ_FEATURE", "PIR_SUMMARY",
  "PIR_SUPERFAMILY", "PRINTS", "PRODOM", "PROSITE", "PUBMED_ID",
  "REACTOME_PATHWAY", "SMART", "SP_COMMENT", "SP_COMMENT_TYPE", "SUPFAM",
  "TIGRFAMS", "UCSC_TFBS", "UNIGENE_EST_QUARTILE", "UP_KEYWORDS",
  "UP_SEQ_FEATURE", "UP_TISSUE"
)

# ---------------------------------------------------------------------------
# Convenience wrapper
# ---------------------------------------------------------------------------

#' Run a complete DAVID annotation clustering workflow
#'
#' Starts a fresh session, authenticates, uploads a gene list, sets
#' annotation categories, and retrieves term cluster results.
#'
#' @param email Character string with a registered DAVID account email.
#' @param gene_list Character or numeric vector of Entrez Gene IDs.
#' @param categories Character vector of annotation categories to use
#'   (default: all categories in \code{DAVID_ANNOTATION_CATEGORIES}).
#' @param debug Logical; if \code{TRUE}, prints progress messages.
#' @return A list with two elements:
#'   \describe{
#'     \item{cluster_count}{Integer number of enrichment clusters found,
#'       or \code{NA} on failure.}
#'     \item{enrichment}{Numeric vector of enrichment scores per cluster,
#'       or \code{NA} if none found.}
#'   }
david.annotate.genes <- function(email, gene_list,
                                 categories = DAVID_ANNOTATION_CATEGORIES,
                                 debug = FALSE) {
  david.new.session()
  if (!david.authenticate(email)) {
    if (debug) message(paste("DAVID: Authentication failed for", email))
    return(list(cluster_count = NA, enrichment = NA))
  }
  list_name <- paste("gene_list", sample(1:10000, 1))
  mapped <- david.add.list(gene_list, "ENTREZ_GENE_ID", list_name, 0L)
  if (is.na(mapped)) {
    if (debug) message("DAVID: Failed to upload gene list")
    return(list(cluster_count = NA, enrichment = NA))
  }
  if (debug) message(paste("DAVID: Mapped", mapped, "of", length(gene_list), "genes"))
  david.set.categories(categories)
  clusters <- david.get.term.cluster.report()
  enrichment <- sapply(clusters, `[[`, "EnrichmentScore")
  if (length(enrichment) == 0) {
    if (debug) message("DAVID: No enrichment clusters found")
    return(list(cluster_count = min(100, length(gene_list)), enrichment = NA))
  }
  return(list(cluster_count = length(clusters), enrichment = enrichment))
}
