get_ensg <- function(target, id1)
{
  # First check if ensg is already given
  if (id1$info[id1$info$trait == target, ]$ontology != target) {
    return(target)
  }

  # Sanitise string for URL
  target <- iconv(target, from = 'UTF-8', to = 'ASCII//TRANSLIT')
  search_url <- "https://api.opentargets.io/v3/platform/public/search?q="
  url <- paste0(search_url, target)
  url <- URLencode(url)
  r <- jsonlite::fromJSON(url)

  if (r$size == 0) {
    ensg <- ""
  } else {
    ensg <- r$data$data$id[1]
  }
  return(ensg)
}

get_efo <- function(disease)
{
  # Sanitise string for URL
  disease <- iconv(disease, from = 'UTF-8', to = 'ASCII//TRANSLIT')
  r <- try(epigraphdb::ontology_gwas_efo(disease, fuzzy = T, mode = "table"), silent = T)

  if (class(r) == "try-error") {
    efo <- ""
  } else {
    efo <- r[r$r.score > 0.95,][1, ]$efo.id
    efo <- tail(strsplit(efo, "/", fixed = T)[[1]], n = 1)
  }
  return(efo)
}

ot_linkage <- function(id1, id2, report)
{
  ensgs_efos <- tidyr::crossing(id1$info$ontology, id2$info$ontology)

  #target_str <- paste0("\"",paste(ensgs$ensg, collapse='", "'), "\"")
  #disease_str <- paste0("\"",paste(efos$efo, collapse='", "'), "\"")
  # Batch query for OT scores
  #result <- httr::POST("https://platform-api.opentargets.io/v3/platform/public/evidence/filter",
  #                     httr::content_type_json(),
  #                     body = paste0('{"target":[', target_str, '], "disease":[', disease_str, ']}'))

  res <- lapply(1:nrow(ensgs_efos), function(i)
  {
    ensg <- ensgs_efos[i, ]$ensg
    efo <- ensgs_efos[i, ]$efo
    url <- paste0("https://api.opentargets.io/v3/platform/public/association/filter?target=", ensg, "&disease=", efo, "&direct=true&facets=false")

    response <- try(r <- jsonlite::fromJSON(url[1]))
    if (class(response) != "try-error" & r$size != 0) {
      res <- tibble::tibble(trait1 = ensgs_efos[i, ]$trait,
                            trait2 = ensgs_efos[i, ]$trait2,
                            overall.ot = r$data$association_score$overall[1],
                            literature.ot = r$data$association_score$datatypes$literature[1],
                            rna_expression.ot = r$data$association_score$datatypes$rna_expression[1],
                            genetic_assoc.ot = r$data$association_score$datatypes$genetic_association[1],
                            somatic_mute.ot = r$data$association_score$datatypes$somatic_mutation[1],
                            known_drug.ot = r$data$association_score$datatypes$known_drug[1],
                            animal_model.ot = r$data$association_score$datatypes$animal_model[1],
                            affected_pathway.ot = r$data$association_score$datatypes$affected_pathway[1])
    } else {
      res <- tibble::tibble(trait1 = ensgs_efos[i, ]$trait,
                            trait2 = ensgs_efos[i, ]$trait2,
                            overall.ot = 0,
                            literature.ot = 0,
                            rna_expression.ot = 0,
                            genetic_assoc.ot = 0,
                            somatic_mute.ot = 0,
                            known_drug.ot = 0,
                            animal_model.ot = 0,
                            affected_pathway.ot = 0)
    }
    return(res)
  }) %>%
    dplyr::bind_rows()

  # Merge with IDs
  traits_merge <- id1$info[c("id", "trait")]
  res <- base::merge(res, traits_merge, by.x = "trait1", by.y = "trait", all.x = T) %>%
    dplyr::rename(id1 = id, name1 = trait1)

  traits_merge <- id2$info[c("id", "trait")]
  res <- base::merge(res, traits_merge, by.x = "trait2", by.y = "trait", all.x = T) %>%
    dplyr::rename(id2 = id, name2 = trait2)

  report$add_otresults(res)
  return(res)
}

dgidb_linkage <- function(id1, id2, report)
{
  interactions <- c("DRUGGABLE GENOME", "CLINICALLY ACTIONABLE", "DRUG RESISTANT")
  gene_names <- id1$info$trait %>%
    paste(., collapse = ",")

  # Set up results first
  res <- id1$info[c("id", "trait")] %>%
    tibble::add_column(Druggable.Genome = F,
                       Clinically.Actionable = F,
                       Drug.Resistant = F)

  drugs <- tibble::tibble(Interaction.Types = character(),
                          Drug.Name = character(),
                          CHEMBL.ID = character(),
                          Sources = character(),
                          PMIDs = character(),
                          Score = numeric(),
                          trait = character(),
                          id = character())

  url <- utils::URLencode(paste0("https://dgidb.org/api/v2/interactions.json?genes=", gene_names))

  response <- try(r <- jsonlite::fromJSON(url), silent = T)
  if (class(response) == "try-error") {
    warning("DGIdb API is under heavy load or unavailable and so evidence will not be retreived. Please try again later.")
    report$add_dgidbresults(NA)
    report$add_dgidbdrugs(NA)
    return()
  }

  if (length(r$matchedTerms)) {
    for (i in 1:length(r$matchedTerms$geneCategories)) {
      res[res$trait == r$matchedTerms$geneName[i], ][c("Druggable.Genome", "Clinically.Actionable", "Drug.Resistant")] <- as.list(interactions %in% r$matchedTerms$geneCategories[[i]]$name)

      temp <- r$matchedTerms$interactions[[i]]
      if (nrow(temp) == 0 | length(temp) == 0) {
        next
      }
      temp$trait <- r$matchedTerms$geneName[i]
      temp$id <- id1$info[id1$info$trait == r$matchedTerms$geneName[i], ]$id

      # Collapse vectors/lists to strings
      temp$sources <- sapply(temp$sources, paste, collapse=", ")
      #temp$pmids <- sapply(temp$pmids, paste, collapse=", ")

      drugs <- temp %>%
        dplyr::rename(Interaction.Types = interactionTypes,
                      Drug.Name = drugName,
                      CHEMBL.ID = drugConceptId,
                      Sources = sources,
                      PMIDs = pmids,
                      Score = score) %>%
        dplyr::select(-interactionId) %>%
        rbind(drugs, .) %>%
        dplyr::relocate(c("trait", "id"), .before = "Interaction.Types")
    }
  }

  report$add_dgidbresults(res)
  report$add_dgidbdrugs(drugs)
}

clinvar_linkage <- function(id1, id2, report)
{
  ensgs_efos <- tidyr::crossing(id1$info$trait, id2$info$trait)
  colnames(ensgs_efos) <- c("ensgs", "efos")


  res <- lapply(1:nrow(ensgs_efos), function(i)
  {
    ensg <- ensgs_efos[i, ]$ensgs
    efo <- ensgs_efos[i, ]$efos

    url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=(",
                  utils::URLencode(ensg), "[Gene+Name])+AND+(",
                  utils::URLencode(gsub("[[:punct:]]", "", efo)), "[Disease%2FPhenotype])&retmode=json")

    r <- jsonlite::fromJSON(url)
  }
  )
}

epigdb_linkage <- function(id1, id2, report)
{
  lapply(1:nrow(id2$info), function(i)
  {
    outcome_trait <- id2$info$trait[i]
    gene_list <- id1$info$trait

    extract_literature <- function(outcome_trait, gene_list) {
      per_gene <- function(gene_name) {
        endpoint <- "/gene/literature"
        params <- list(
          gene_name = gene_name,
          object_name = outcome_trait %>% stringr::str_to_lower()
        )
        df <- epigraphdb::query_epigraphdb(route = endpoint, params = params, mode = "table")
      }
      res_df <- gene_list %>% purrr::map_df(per_gene)

      if (length(res_df) == 0) {
        report$add_epidbresults(NULL)
        return()
      }

      res_df %>%
        mutate(literature_count = purrr::map_int(pubmed_id, function(x) length(x)))
    }

    literature_df <- extract_literature(
      outcome_trait = outcome_trait,
      gene_list = gene_list
    )

    if (is.null(literature_df)) {
      report$add_epidbresults(NULL)
      return()
    }

    # Append our identifiers to the df
    identifiers <- id1$info %>%
      dplyr::select(id, trait)
    literature_df <- base::merge(literature_df, identifiers, by.x = "gene.name", by.y = "trait", all.x = T, all.y = F) %>%
      dplyr::rename(id1 = id)
    literature_df$id2 <- id2$info$id[i]

    # Add to report
    report$add_epidbresults(literature_df)
  }
  )
}

pmid_to_link <- function(id, hyperlink = F)
{
  url <- paste0("https://pubmed.ncbi.nlm.nih.gov/", id, "/")

  if (hyperlink == T) {
    url <- paste0("<a href='", url, "'>", id, "</a>")
  }

  return(url)
}

intermine <- function(id1, id2, report)
{
  mgis <- .ensg_to_mgi(id1$info$ontology)

  lapply(1:nrow(mgis), function(i)
  {
    response <- try(r <- RCurl::getURL(paste0('https://www.ebi.ac.uk/mi/impc/solr/phenodigm/select?rows=250000&q=type:disease_model_summary%20AND%20marker_id:"',
                                mgis[i, 3],
                                '"&wt=csv&fl=marker_id,model_id,model_description,model_genetic_background,disease_id,disease_term,association_curated,disease_model_avg_norm,disease_model_max_norm')))

    if (class(response) == "try-error") {
      return()
    }

    res <- read.csv(text = r)

    report$add_mouseresults(res)
  })
}
