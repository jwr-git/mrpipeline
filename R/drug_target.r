#' Generate drug target evidence
#'
#' Uses the Drug Genome Interaction DB's API to search for drug target-related
#' evidence, including on: Druggable Genome, Clinically Actionable and
#' Drug Resistant ontologies.
#'
#' The lookup MUST be ENSG IDs.
#'
#' @param dat A data.frame or named list
#' @param ensg_col Column, or name, to be accessed in `dat`
#'
#' @return data.frame of results
#' @export
drug_target_evidence <- function(dat,
                                 ensg_col = "exposure")
{
  dgi <- TRUE

  if (!(ensg_col %in% names(dat))) {
    warning("Could not find column for ENSG IDs named \"", ensg_col)
    dgi <- FALSE
  }

  if (dgi == TRUE)
  {
    # Split ENSGS if too large for API
    ensgs <- unique(dat[[ensg_col]])
    if (length(ensgs) > 100) {
      search_ensgs <- .chunk(ensgs, length(ensgs) / 50)

      dgidb_res <- lapply(1:length(search_ensgs), function(i)
      {
        .dgidb_linkage(search_ensgs[[i]])
      }) %>% dplyr::bind_rows()
    } else {
      dgidb_res <- .dgidb_linkage(ensgs)
    }
  }

  return(dgidb_res)
}

#' Link ENSG IDs with DGIdb
#'
#' Lookup ENSGs using the Drug Genome Interaction DB API.
#'
#' @param ensgs Vector of ENSG IDs
#'
#' @return data.frame of results
#' @keywords Internal
.dgidb_linkage <- function(ensgs)
{
  interactions <- c("DRUGGABLE GENOME", "CLINICALLY ACTIONABLE", "DRUG RESISTANT")
  gene_names <- paste(ensgs, collapse = ",")

  terms <- tibble::as_tibble(ensgs) %>%
    tibble::add_column(Druggable.Genome = FALSE,
                       Clinically.Actionable = FALSE,
                       Drug.Resistant = FALSE) %>%
    dplyr::rename(Gene = value)

  url <- utils::URLencode(paste0("https://dgidb.org/api/v2/interactions.json?genes=", gene_names))

  response <- try(r <- jsonlite::fromJSON(url), silent = T)
  if (class(response) == "try-error") {
    warning("DGIdb API is under heavy load or unavailable and so evidence will not be retreived. Please try again later.")
    return(NULL)
  }

  if (length(r$matchedTerms)) {
    for (i in 1:length(r$matchedTerms$geneCategories)) {
      terms[terms$Gene == r$matchedTerms$searchTerm[i], ][c("Druggable.Genome", "Clinically.Actionable", "Drug.Resistant")] <- as.list(interactions %in% r$matchedTerms$geneCategories[[i]]$name)
    }
  }

  return(terms)
}

#' Splits vector into chunks
#'
#' @param x Vector
#' @param n Number of chunks to create
#'
#' @return list of chunks
#' @keywords Internal
.chunk <- function(x, n)
{
  split(x, cut(seq_along(x), n, labels = FALSE))
}
