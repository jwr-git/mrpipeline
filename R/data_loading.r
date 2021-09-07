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

    # Load locally from vcf
    if (file.exists(f) & tools::file_ext(f) == "vcf")
    {
      dat <- gwasvcf::query_gwas(vcf = VariantAnnotation::readVcf(f), pval = -log10(conf$p_cutoff)) %>%
        gwasglue::gwasvcf_to_TwoSampleMR("exposure") %>%
        dplyr::mutate(exposure = id1$info$trait[i],
                      trait.exposure = id1$info$trait[i],
                      id.exposure = id1$info$trait[i])

      # Reporting
      if (length(dat)) {
        # F-statistic
        ret <- .calc_f_stat(dat, conf$f_cutoff)
        dat <- ret[[1]]
        rem.f.stat <- ret[[2]]

        # Annotate
        id1$annotate_ensg(id1$info$id[i], unique(dat$exposure))

        pre_clump <- nrow(dat)

        dat <- dplyr::mutate(dat,
                             rsid = SNP,
                             pval = pval.exposure) %>%
          ieugwasr::ld_clump(clump_kb = conf$clump_kb,
                             clump_r2 = conf$clump_r2,
                             pop = conf$pop,
                             bfile = conf$bfile_path,
                             plink_bin = conf$plink_path)
        #dat %<>% TwoSampleMR::clump_data(clump_kb = kb, clump_r2 = r2, pop = pop)

        #f_plot(dat, report)

        report$add_dataset(report$get_exposures_name(),
                           c(unique(dat$exposure),
                             id1$info$id[i],
                             id1$info$trait[i],
                             pre_clump,
                             F,
                             rem.f.stat,
                             pre_clump - nrow(dat),
                             nrow(dat),
                             0))
      } else {
        dat <- NULL

        warning("No SNPs selected for exposure: ", id1$info$trait[i])

        report$add_dataset(report$get_exposures_name(),
                           c(unique(dat$exposure),
                             id1$info$id[i],
                             id1$info$trait[i],
                             0,
                             F))
      }
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

      if (length(dat)) {
        dat <- dat %>%
          tidyr::separate(exposure, "exposure", sep = "\\|\\|", extra = "drop") %>%
          dplyr::mutate(id.exposure = id1$info$trait[i])

        dat$exposure <- trimws(dat$exposure)

        if (!is.null(conf$bfile_path) & !is.null(conf$plink_path))
        {
          dat <- dplyr::mutate(dat,
                               rsid = SNP,
                               pval = pval.exposure) %>%
            ieugwasr::ld_clump(clump_kb = conf$clump_kb,
                               clump_r2 = conf$clump_r2,
                               pop = conf$pop,
                               bfile = conf$bfile_path,
                               plink_bin = conf$plink_path)
        }

        # F-statistic
        ret <- .calc_f_stat(dat, conf$f_cutoff)
        dat <- ret[[1]]
        rem.f.stat <- ret[[2]]

        # Annotate
        id1$annotate_ensg(id1$info$id[i], unique(dat$exposure))

        #f_plot(dat, report)

        report$add_dataset(report$get_exposures_name(),
                           c(unique(dat$exposure),
                             id1$info$id[i],
                             id1$info$trait[i],
                             nrow(dat),
                             ifelse(is.null(conf$bfile_path) & is.null(conf$plink_path), T, F),
                             rem.f.stat,
                             0,
                             nrow(dat),
                             1))
      } else {
        dat <- NULL

        warning("No SNPs selected for exposure: ", id1$info$trait[i])

        report$add_dataset(report$get_exposures_name(),
                           c(unique(dat$exposure),
                             id1$info$id[i],
                             id1$info$trait[i],
                             0,
                             F))
      }
    }

    return (dat)
  }, mc.cores = conf$cores) %>%
    dplyr::bind_rows()

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
      dat <- gwasvcf::query_gwas(vcf = VariantAnnotation::readVcf(f), rsid = rsids) %>%
        gwasglue::gwasvcf_to_TwoSampleMR("outcome") %>%
        dplyr::mutate(id.outcome = id2$info$trait[i])

      if (length(dat)) {
        # Annotate EFO
        id2$annotate_efo(id2$info$trait[i], unique(dat$outcome))

        report$add_dataset(report$get_outcomes_name(),
                           c(unique(dat$outcome),
                             id2$info$id[i],
                             id2$info$trait[i],
                             nrow(dat),
                             0,
                             nrow(dat)))
      } else {
        dat <- NULL

        warning("No SNPs selected for outcome: ", id2$info$trait[i])

        report$add_dataset(report$get_outcomes_name(),
                           c(unique(dat$outcome),
                             id2$info$id[i],
                             id2$info$trait[i],
                             0,
                             0,
                             0))
      }
    } else
    {
      # From OpenGWAS
      dat <- TwoSampleMR::extract_outcome_data(rsids, id2$info$id[i], proxies = conf$proxies) %>%
        tidyr::separate(outcome, "outcome", sep = "\\|\\|", extra = "drop") %>%
        dplyr::mutate(id.outcome = id2$info$trait[i])

      if (length(dat)) {
        dat$outcome <- trimws(dat$outcome)

        # Annotate EFO
        id2$annotate_efo(id2$info$trait[i], unique(dat$outcome))

        report$add_dataset(report$get_outcomes_name(),
                           c(unique(dat$outcome),
                             id2$info$id[i],
                             id2$info$trait[i],
                             nrow(dat),
                             nrow(dat[!is.na(dat$proxy.outcome) & dat$proxy.outcome == T,]),
                             nrow(dat),
                             1))
      } else {
        dat <- NULL

        warning("No SNPs selected for outcome: ", id2$info$trait[i])

        report$add_dataset(report$get_exposures_name(),
                           c(unique(dat$outcome),
                             id2$info$id[i],
                             id2$info$trait[i],
                             0,
                             F,
                             0,
                             0,
                             0))
      }
    }

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

                          annotate_ensg = function(id, name) {
                            "Annotates IDs, currently using biomaRt for ENSGs"
                            name <- fs::path_sanitize(name)
                            name_san <- stringr::str_extract(name, "ENSG[0-9]+")

                            # If trait name is ENSG, annotate these using biomaRt to hgnc_symbol
                            if (!is.na(name_san)) {
                              hgnc <- .ensg_to_name(name_san)

                              if (length(hgnc)) {
                                info[info$id == id, ]$trait <<- hgnc$hgnc_symbol
                              } else {
                                info[info$id == id, ]$trait <<- name
                              }

                              info[info$id == id, ]$ontology <<- name_san
                            } else {
                              info[info$id == id, ]$trait <<- name
                              info[info$id == id, ]$ontology <<- name
                            }
                          },

                          annotate_efo = function(id, name) {
                            "Annotates EFO"
                            disease <- iconv(name, from = 'UTF-8', to = 'ASCII//TRANSLIT')
                            r <- try(epigraphdb::ontology_gwas_efo(disease, fuzzy = T, mode = "table"), silent = T)

                            if (class(r) == "try-error") {
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
.ensg_to_name <- function(ensg)
{
  if (!require("biomaRt"))
  {
    warning("biomaRt not found! To annotate ENSG IDs, the biomaRt package must be installed.")
    return(NULL)
  }

  library(biomaRt)

  mart.gene <- NA

  # Just in case biomaRt is down, this would stop the entire pipeline
  tryCatch({
    mart.gene <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                                  host="grch37.ensembl.org",
                                  path="/biomart/martservice",
                                  dataset="hsapiens_gene_ensembl")
  }, error = function(e_code) {
    warning("biomaRt is currently unavailable and so ENSG IDs will not be annotated.")
  })

  if (exists("mart.gene")) {
    if (class(mart.gene) != "Mart") {
      return(NULL)
    }
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
    return(NULL)
  }

  library(biomaRt)

  mart.mouse <- NA
  mart.human <- NA

  # Just in case biomaRt is down, this would stop the entire pipeline
  tryCatch({
    mart.mouse <- useMart(biomart="ensembl",
                          dataset="mmusculus_gene_ensembl")
    mart.human <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                          dataset="hsapiens_gene_ensembl")

  }, error = function(e_code) {
    warning("biomaRt is currently unavailable and so ENSG IDs will not be converted to mouse gene IDs.")
  })

  results <- biomaRt::getLDS(attributes = c("ensembl_gene_id"),
                             filters = "ensembl_gene_id",
                             values = ensg,
                             mart = mart.human,
                             attributesL = c("mgi_symbol", "mgi_id"),
                             martL = mart.mouse,
                             uniqueRows=T)

  return(results)
}
