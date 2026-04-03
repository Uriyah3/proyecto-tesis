#' Load credentials from .env file
#'
#' Reads key=value pairs from .env and sets them as environment variables.
#' Falls back to mail_list.txt for DAVID email (legacy support).
#'
load.credentials <- function(env_file = ".env") {
  if (file.exists(env_file)) {
    lines <- readLines(env_file, warn = FALSE)
    for (line in lines) {
      line <- trimws(line)
      if (nchar(line) == 0 || startsWith(line, "#")) next
      parts <- strsplit(line, "=", fixed = TRUE)[[1]]
      if (length(parts) >= 2) {
        key <- trimws(parts[1])
        value <- trimws(paste(parts[-1], collapse = "="))
        do.call(Sys.setenv, setNames(list(value), key))
      }
    }
  }
}

#' Get DAVID registered email
#'
#' Checks .env first (DAVID_EMAIL), then falls back to mail_list.txt
#'
get.david.email <- function() {
  load.credentials()
  email <- Sys.getenv("DAVID_EMAIL", "")
  if (nchar(email) > 0) return(email)

  # Legacy fallback
  if (file.exists("mail_list.txt")) {
    emails <- readLines("mail_list.txt", warn = FALSE)
    emails <- emails[nchar(trimws(emails)) > 0]
    if (length(emails) > 0) return(emails[1])
  }

  stop("No DAVID email found. Set DAVID_EMAIL in .env or add email to mail_list.txt.\n",
       "Register at: https://davidbioinformatics.nih.gov/webservice/register.htm")
}

#' Get DisGeNET API key
#'
get.disgenet.api.key <- function() {
  load.credentials()
  key <- Sys.getenv("DISGENET_API_KEY", "")
  if (nchar(key) > 0) return(key)

  stop("No DisGeNET API key found. Set DISGENET_API_KEY in .env\n",
       "Register at: https://disgenet.com/academic-apply")
}
