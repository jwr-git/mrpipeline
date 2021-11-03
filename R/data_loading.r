#' Reads exposure and outcome datasets from local .vcf files
#' or OpenGWAS DB
#'
#' @param id1 DatasetsID class of exposure IDs
#' @param id2 DatasetsID class of outcome IDs
#' @param report QCReport class of results, etc. for reporting
#' @param conf config::config file of parameters
#'
#' @return Data.frame of harmonised exposure-outcome SNPs
read_datasets <- function(id1, id2,
                          report,
                          conf)
{
  # Bad practice to loop over an object which may be altered by that loop
  rows_id1 <- nrow(id1$info)

  exposure_dat <- parallel::mclapply(1:rows_id1, function(i)
  {
    f <- id1$info$filename[i]

    # Load locally from vcf or gwas OpenGWAS
    # and convert to TwoSampleMR format ready for
    # the rest of the analyses
    if (file.exists(f) & tools::file_ext(f) == "vcf")
    {
      dat <- gwasvcf::query_gwas(vcf = VariantAnnotation::readVcf(f), pval = -log10(conf$p_cutoff))

      if (nrow(dat) < 1) {
        warning("No SNPs selected for exposure: ", id1$info$trait[i])
        report$add_dataset(report$get_exposures_name(),
                           c("",
                             id1$info$id[i],
                             id1$info$trait[i],
                             0,
                             F))
        id1$bad_id(i)
        return(NULL)
      }

      dat <- dat %>% gwasglue::gwasvcf_to_TwoSampleMR("exposure")
    } else
    {
      # From OpenGWAS
      dat <- TwoSampleMR::extract_instruments(id1$info$id[i],
                                              p1 = conf$p_cutoff,
                                              clump = ifelse(is.null(conf$bfile_path) & is.null(conf$plink_path),
                                                             T,
                                                             F),
                                              r2 = conf$clump_r2,
                                              kb = conf$clump_kb)

      if (nrow(dat) < 1) {
        warning("No SNPs selected for exposure: ", id1$info$trait[i])
        report$add_dataset(report$get_exposures_name(),
                           c(unique(dat$exposure),
                             id1$info$id[i],
                             id1$info$trait[i],
                             0,
                             F))
        id1$bad_id(i)
        return(NULL)
      }

      # Some cleaning due to OpenGWAS format
      dat <- dat %>%
        tidyr::separate(exposure, "exposure", sep = "\\|\\|", extra = "drop")
      dat$exposure <- trimws(dat$exposure)
    }

    # "Our" exposure name
    dat$info.exposure <- id1$info$trait[i]

    ############################################################################
    # F-statistic test
    ############################################################################
    ret <- calc_f_stat(dat, conf$f_cutoff)
    dat <- ret[[1]]
    rem.f.stat <- ret[[2]]

    if (nrow(dat) < 1) {
      warning("No SNPs passed the F-statistic threshold for: ", id1$info$trait[i])
      report$add_dataset(report$get_exposures_name(),
                         c("",
                           id1$info$id[i],
                           id1$info$trait[i],
                           0,
                           F))
      id1$bad_id(i)
      return(NULL)
    }

    ############################################################################
    # Local clumping, if specified
    # Both datasets will be in TwoSampleMR format by this stage
    ############################################################################
    pre_clump <- nrow(dat)
    if (!is.null(conf$bfile_path) & !is.null(conf$plink_path))
    {
      attempt <- try(dat <- dplyr::mutate(dat,
                                          rsid = SNP,
                                          pval = pval.exposure) %>%
                       ieugwasr::ld_clump(clump_kb = conf$clump_kb,
                                          clump_r2 = conf$clump_r2,
                                          pop = conf$pop,
                                          bfile = conf$bfile_path,
                                          plink_bin = conf$plink_path), silent = T)

      # Clumping failed?
      if (class(attempt) == "try-error" | nrow(dat) < 1) {
        warning("Error with clumping, or no SNPs retained thereafter for exposure: ", id1$info$trait[i])
        report$add_dataset(report$get_exposures_name(),
                           c("",
                             id1$info$id[i],
                             id1$info$trait[i],
                             0,
                             F))
        id1$bad_id(i)
        return(NULL)
      }
    }

    report$add_dataset(report$get_exposures_name(),
                       c(unique(dat$exposure),
                         id1$info$id[i],
                         id1$info$trait[i],
                         nrow(dat),
                         ifelse(file.exists(f) & tools::file_ext(f) == "vcf", F,
                           ifelse(is.null(conf$bfile_path) & is.null(conf$plink_path), T, F)),
                         rem.f.stat,
                         pre_clump - nrow(dat),
                         nrow(dat),
                         !(file.exists(f) & tools::file_ext(f) == "vcf")))

    return (dat)
  }, mc.cores = conf$cores) %>%
    dplyr::bind_rows()

  # Some tidying up again!
  # Remove those IDs with no data from future analyses
  id1$remove_bad_ids()

  rsids <- unique(exposure_dat$SNP)

  rows_id2 <- nrow(id2$info)
  outcome_dat <- parallel::mclapply(1:rows_id2, function(i)
  {
    f <- id2$info$filename[i]

    # Is local vcf file? Load from there
    if (file.exists(f))
    {
      # Set data as null so that we can tell if the local proxy search
      # failed or did not run
      dat <- NULL

      if (conf$proxies == TRUE) {
        if (!is.null(conf$bfile_path)) {
          dat <- gwasvcf::query_gwas(vcf = VariantAnnotation::readVcf(f),
                                     proxies = "yes", bfile = conf$bfile_path,
                                     rsid = rsids)
        }
        else if (!is.null(conf$dbfile_path)) {
          dat <- gwasvcf::query_gwas(vcf = VariantAnnotation::readVcf(f),
                                     proxies = "yes", dbfile = conf$dbfile_path,
                                     rsid = rsids)
        }
      }

      # If this is null, then no local search for proxies
      if (is.null(dat)) {
        dat <- gwasvcf::query_gwas(vcf = VariantAnnotation::readVcf(f), rsid = rsids)
      }

      if (nrow(dat)) {
        dat <- dat %>%
          gwasglue::gwasvcf_to_TwoSampleMR("outcome")
      }
    } else
    {
      # From OpenGWAS
      dat <- TwoSampleMR::extract_outcome_data(rsids, id2$info$id[i], proxies = conf$proxies)

      if (nrow(dat)) {
        dat <- dat %>%
          tidyr::separate(outcome, "outcome", sep = "\\|\\|", extra = "drop")

        dat$outcome <- trimws(dat$outcome)
      }
    }

    # No SNPs
    if (nrow(dat) < 1) {
      dat <- NULL

      warning("No SNPs selected for outcome: ", id2$info$trait[i])

      report$add_dataset(report$get_outcomes_name(),
                         c(unique(dat$outcome),
                           id2$info$id[i],
                           id2$info$trait[i],
                           0,
                           F,
                           0,
                           0,
                           0))
    }


    # Standardise outcome data
    if (!"proxy.outcome" %in% names(dat)) {
      dat$proxy.outcome <- F
      dat$proxy_snp.outcome <- ""
    }

    # EFO linkage - do it here as we may need the exposure name from
    # IEUGWAS DB if data are extracted from there
    id2$annotate_efo(id2$info$trait[i], unique(dat$outcome))

    report$add_dataset(report$get_outcomes_name(),
                       c(unique(dat$outcome),
                         id2$info$id[i],
                         id2$info$trait[i],
                         nrow(dat),
                         nrow(dat[!is.na(dat$proxy.outcome) & dat$proxy.outcome == T,]),
                         nrow(dat),
                         1))

    return(dat)
  }, mc.cores = conf$cores) %>%
    dplyr::bind_rows()

  dat <- TwoSampleMR::harmonise_data(exposure_dat, outcome_dat, action = conf$harmonise_action)
  return(dat)
}

#' A Reference Class for the dataset IDs
#'
#' @field info A data.frame of details of each dataset, including filename,
#'             annotated trait name and ID if from OpenGWAS DB
DatasetIDs <- setRefClass("DatasetIDs",
                        fields = list(info = "data.frame"),
                        methods = list(
                          initialise = function(ids) {
                            info <<- dplyr::tibble(id = ids, trait = ids, filename = ids, ontology = ids, good = T)
                            info <<- info %>% dplyr::filter(duplicated(id) == F)

                            loc <- file.exists(ids)
                            if (any(loc))
                            {
                              info[loc,]$filename <<- file.path(info[loc,]$filename)
                              info[loc,]$trait <<- tools::file_path_sans_ext(basename(info[loc,]$filename))
                              info[loc,]$id <<- tools::file_path_sans_ext(basename(info[loc,]$filename))
                              info[loc,]$ontology <<- tools::file_path_sans_ext(basename(info[loc,]$filename))
                            }
                          },

                          add_id = function(id, trait, filename) {
                            "Adds new ID with trait and filename"
                            info <<- info %>%
                              tibble::add_row(id = id, trait = trait, filename = filename, ontology = trait, good = rep(T, length(id)))
                          },

                          bad_id = function(id) {
                            "Notes that an ID is 'bad', i.e. no data loaded"
                            info[info$id == id, ]$good <<- F
                          },

                          remove_bad_ids = function() {
                            "Removes those IDs from the data.frame which are 'bad'"
                            bad <- info[info$good == F, ]
                            info <<- info[info$good == T, ]

                            if (nrow(bad)) {
                              warning("The following datasets were removed due to there being no SNPs present for analysis: ", paste(unique(bad$id), collapse=", "))
                            }
                          },

                          annotate_ensg = function() {
                            "Annotates IDs, currently using biomaRt for ENSGs"
                            lookup <- info[grepl("eqtl-a-ENSG[0-9]+", info$id), ]$id

                            if (length(lookup) < 1) {
                              return()
                            }

                            lookup_ensg <- stringr::str_extract(lookup, "ENSG[0-9]+")

                            if (length(lookup_ensg)) {
                              hgnc <- ensg_to_name(lookup_ensg)
                              for (i in 1:nrow(hgnc)) {
                                info[grepl(paste0("eqtl-a-", hgnc[i, "ensembl_gene_id"]), info$id), ]$trait <<- hgnc[i, "hgnc_symbol"]
                              }
                            }

                            info <<- info %>% dplyr::mutate(trait = ifelse(trait %in% "", id, trait))
                          },

                          annotate_efo = function(id, name) {
                            "Annotates EFO"
                            disease <- iconv(name, from = 'UTF-8', to = 'ASCII//TRANSLIT')
                            attempt <- try(r <- epigraphdb::ontology_gwas_efo(disease, fuzzy = T, mode = "table"), silent = T)

                            if (inherits(attempt, "try-error")) {
                              efo <- name
                            } else {
                              efo <- r[r$r.score > 0.95,][1, ]$efo.id
                              efo <- tail(strsplit(efo, "/", fixed = T)[[1]], n = 1)
                            }

                            info[info$id == id, ]$trait <<- name
                            info[info$id == id, ]$ontology <<- efo
                          }
                        )
)

#' Pairwise combination of two DatasetIDs classes
#'
#' @param id1 DatasetsIDs class of exposure IDs
#' @param id2 DatasetsIDs class of outcome IDs
#'
#' @return data.frame
.combine_ids <- function(id1, id2)
{
  # Not the most efficient or R-like but since these should be relatively small,
  # I think it's fine to do this. tidyr::crossing does something similar but
  # renaming the cols is weird

  df = data.frame(id.x = character(), trait.x = character(), filename.x = character(), ontology.x = character(), good.x = logical(),
                  id.y = character(), trait.y = character(), filename.y = character(), ontology.y = character(), good.y = logical())
  for (r1 in 1:nrow(id1$info))
  {
    for (r2 in 1:nrow(id2$info))
    {
      df[nrow(df) + 1, ] <- c(id1$info[r1, ], id2$info[r2, ])
    }
  }
  #return(df %<>% dplyr::rename())
  return(df)
}

#' Converts ENSG to gene name using biomaRt
#'
#' @param ensg Vector of or single character, ENSG IDs
#'
#' @return Vector of or single character, gene names
ensg_to_name <- function(ensg)
{
  if (!require("biomaRt"))
  {
    warning("biomaRt not found! To annotate ENSG IDs, the biomaRt package must be installed.")
    return(NA)
  }

  mart.gene <- NA

  # Just in case biomaRt is down, this would stop the entire pipeline
  attempt <- tryCatch({
    mart.gene <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                                  host="grch37.ensembl.org",
                                  path="/biomart/martservice",
                                  dataset="hsapiens_gene_ensembl")
  }, error = function(e_code) {
    warning("biomaRt is currently unavailable and so ENSG IDs will not be annotated.")
  })

  if (inherits(attempt, "error")) {
    return(NA)
  }

  results <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'),
                   filters = 'ensembl_gene_id',
                   values = ensg,
                   mart = mart.gene)

  return(results)
}

.ensg_to_mgi <- function(ensg)
{
  if (!require("biomaRt"))
  {
    warning("biomaRt not found! To annotate ENSG IDs, the biomaRt package must be installed.")
    return(NA)
  }

  mart.mouse <- NA
  mart.human <- NA

  # Just in case biomaRt is down, this would stop the entire pipeline
  attempt <- tryCatch({
    mart.mouse <- biomaRt::useMart(biomart="ensembl",
                          dataset="mmusculus_gene_ensembl")
    mart.human <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                          dataset="hsapiens_gene_ensembl")

  }, error = function(e) {
    warning("biomaRt is currently unavailable and so ENSG IDs will not be converted to mouse gene IDs.")
  })

  if (inherits(attempt, "error")) {
    return(NA)
  }

  results <- biomaRt::getLDS(attributes = c("ensembl_gene_id"),
                             filters = "ensembl_gene_id",
                             values = ensg,
                             mart = mart.human,
                             attributesL = c("mgi_symbol", "mgi_id"),
                             martL = mart.mouse,
                             uniqueRows=T)

  return(results)
}
