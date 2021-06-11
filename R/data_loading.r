#' Reads exposure and outcome datasets from local .vcf files
#' or OpenGWAS DB
#'
#' @param id1 DatasetsID class of exposure IDs
#' @param id2 DatasetsID class of outcome IDs
#' @param report QCReport class of results, etc. for reporting
#' @param p_cutoff P value cutoff for "tophits" from summary statistics
#' @param chrompos_query Select SNPs within this region only,
#'                       format: "chrom:start-end"
#' @param gene_query Select SNPs nearby this gene only, NOT YET IMPLEMENTED
#' @param f_cutoff F-statistic cutoff for exposure data, default: 10
#' @param r2 Clumping threshold, default: 0.001
#' @param kb Clumping window, default: 10000
#' @param pop Super-population for clumping, default: EUR
#' @param proxies Whether to search for and use proxy SNPs in the outcome
#'                dataset, default: TRUE
#' @param action See: TwoSampleMR::harmonise_data, default: 1
#' @param bfile Plink 1.x format reference files for use in local clumping,
#'              etc, default: NULL
#' @param plink_bin Plink 1.x location, default: NULL
#' @param nthreads Number of threads for use in multithreaded functions,
#'                 default: 1
#'
#' @return Data.frame of harmonised exposure-outcome SNPs
read_datasets <- function(id1, id2,
                          report,
                          p_cutoff, chrompos_query, gene_query,
                          f_cutoff = 10,
                          r2 = 0.001, kb = 10000, pop = "EUR",
                          proxies = T,
                          action = 1,
                          bfile = NULL, plink_bin = NULL,
                          nthreads = 1)
{
  # Bad practice to loop over an object which may be altered by that loop
  rows_id1 <- nrow(id1$info)

  exposure_dat <- parallel::mclapply(1:rows_id1, function(i)
  {
    f <- id1$info$filename[i]

    # Load locally from vcf
    if (file.exists(f) & tools::file_ext(f) == "vcf")
    {
      dat <- gwasvcf::query_gwas(vcf = VariantAnnotation::readVcf(f), pval = -log10(p_cutoff)) %>%
        gwasglue::gwasvcf_to_TwoSampleMR("exposure") %>%
        dplyr::mutate(exposure = id1$info$trait[i],
                      trait.exposure = id1$info$trait[i],
                      id.exposure = id1$info$trait[i])

      # F-statistic
      ret <- .calc_f_stat(dat, f_cutoff)
      dat <- ret[[1]]
      rem.f.stat <- ret[[2]]

      # Annotate
      id1$annotate(id1$info$id[i], unique(dat$exposure))

      # Reporting
      if (length(dat)) {
        pre_clump <- nrow(dat)

        dat <- dplyr::mutate(dat,
                             rsid = SNP,
                             pval = pval.exposure) %>%
          ieugwasr::ld_clump(clump_kb = kb, clump_r2 = r2, pop = pop,
                             bfile = bfile, plink_bin = plink_bin)
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
                                              p1 = p_cutoff,
                                              clump = ifelse(is.null(bfile) & is.null(plink_bin),
                                                             T,
                                                             F),
                                              r2 = r2,
                                              kb = kb)

      # Nothing extracted - remove ID and return
      if (length(dat) < 1) {
        id1$bad_id(id1$info$id[i])
        return(NULL)
      }

      dat <- dat %>%
        tidyr::separate(exposure, "exposure", sep = "\\|\\|", extra = "drop") %>%
        dplyr::mutate(id.exposure = id1$info$trait[i])

      dat$exposure <- trimws(dat$exposure)

      if (!is.null(bfile) & !is.null(plink_bin))
      {
        dat <- dplyr::mutate(dat,
                             rsid = SNP,
                             pval = pval.exposure) %>%
          ieugwasr::ld_clump(clump_kb = kb, clump_r2 = r2, pop = pop,
                             bfile = bfile, plink_bin = plink_bin)
      }

      # F-statistic
      ret <- .calc_f_stat(dat, f_cutoff)
      dat <- ret[[1]]
      rem.f.stat <- ret[[2]]

      # Annotate
      id1$annotate(id1$info$id[i], unique(dat$exposure))

      if (length(dat)) {
        #f_plot(dat, report)

        report$add_dataset(report$get_exposures_name(),
                           c(unique(dat$exposure),
                             id1$info$id[i],
                             id1$info$trait[i],
                             nrow(dat),
                             ifelse(is.null(bfile) & is.null(plink_bin), T, F),
                             rem.f.stat,
                             0,
                             nrow(dat),
                             1))
      } else {
        report$add_dataset(report$get_exposures_name(),
                           c(unique(dat$exposure),
                             id1$info$id[i],
                             id1$info$trait[i],
                             0,
                             F))
      }
    }

    return (dat)
  }, mc.cores = nthreads) %>%
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
        report$add_dataset(report$get_outcomes_name(),
                           c(unique(dat$outcome),
                             id2$info$id[i],
                             id2$info$trait[i],
                             nrow(dat),
                             0,
                             nrow(dat)))
      } else {
        dat <- NULL

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
      dat <- TwoSampleMR::extract_outcome_data(rsids, id2$info$id[i], proxies = proxies) %>%
        tidyr::separate(outcome, "outcome", sep = "\\|\\|", extra = "drop") %>%
        dplyr::mutate(id.outcome = id2$info$trait[i])

      dat$outcome <- trimws(dat$outcome)

      # Annotate
      id2$annotate(id2$info$id[i], unique(dat$outcome))

      if (length(dat)) {
        report$add_dataset(report$get_outcomes_name(),
                           c(unique(dat$outcome),
                             id2$info$id[i],
                             id2$info$trait[i],
                             nrow(dat),
                             nrow(dat[!is.na(dat$proxy.outcome) & dat$proxy.outcome == T,]),
                             nrow(dat),
                             1))
      } else {
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
  }, mc.cores = nthreads) %>%
    dplyr::bind_rows()

  dat <- TwoSampleMR::harmonise_data(exposure_dat, outcome_dat, action = action)
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
                            info <<- dplyr::tibble(id = ids, trait = ids, filename = ids, good = T)
                            info <<- info %>% dplyr::filter(duplicated(id) == F)

                            loc <- file.exists(ids)
                            if (any(loc))
                            {
                              info[loc,]$filename <<- file.path(info[loc,]$filename)
                              info[loc,]$trait <<- tools::file_path_sans_ext(basename(info[loc,]$filename))
                              info[loc,]$id <<- tools::file_path_sans_ext(basename(info[loc,]$filename))
                            }
                          },

                          add_id = function(id, trait, filename) {
                            "Adds new ID with trait and filename"
                            info <<- info %>%
                              tibble::add_row(id = id, trait = trait, filename = filename, good = rep(T, length(id)))
                          },

                          bad_id = function(id) {
                            "Notes that an ID is 'bad', i.e. no data loaded"
                            info[info$id == id, ]$good <<- F
                          },

                          remove_bad_ids = function() {
                            "Removes those IDs from the data.frame which are 'bad'"
                            info <<- info[info$good == T, ]
                          },

                          annotate = function(id, name) {
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
                            } else {
                              info[info$id == id, ]$trait <<- name
                            }
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

  df = data.frame(id.x = character(), trait.x = character(), filename.x = character(), good.x = logical(),
                  id.y = character(), trait.y = character(), filename.y = character(), good.y = logical())
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
    mart.gene <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                        host="grch37.ensembl.org",
                        path="/biomart/martservice",
                        dataset="hsapiens_gene_ensembl")
  }, error = function(e_code) {
    warning("biomaRt is currently unavailable and so ENSG IDs will not be annotated.")
  })

  if (!exists("mart.gene")) {
    return(NULL)
  }

  results <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'),
                   filters = 'ensembl_gene_id',
                   values = ensg,
                   mart = mart.gene)

  return(results)
}
