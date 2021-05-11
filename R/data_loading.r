#' Reads all datasets
#'
#'
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
  rows_id1 <- nrow(id1$info)

  exposure_dat <- parallel::mclapply(1:rows_id1, function(i)
  {
    f <- id1$info$filename[i]

    # Load locally from vcf
    if (file.exists(f))
    {
      dat <- gwasvcf::query_gwas(vcf = VariantAnnotation::readVcf(f), pval = p_cutoff) %>%
        gwasglue::gwasvcf_to_TwoSampleMR("exposure") %>%
        dplyr::mutate(exposure = id1$info$trait[i],
                      trait.exposure = id1$info$trait[i],
                      id.exposure = id1$info$trait[i])

      # PVE / F-statistic
      dat$maf.exposure <- ifelse(dat$eaf.exposure < 0.5, dat$eaf.exposure, 1 - dat$eaf.exposure)
      dat$pve.exposure <- mapply(.calc_pve,
                                 dat$beta.exposure,
                                 dat$maf.exposure,
                                 dat$se.exposure,
                                 dat$samplesize.exposure)
      dat$f.stat.exposure <- ((dat$samplesize.exposure - 1 - 1) / 1) * (dat$pve.exposure / (1 - dat$pve.exposure))

      #dat$f.stat.exposure <- dat$beta.exposure ** 2 / dat$se.exposure ** 2
      rem.f.stat <- nrow(dat[dat$f.stat.exposure < f_cutoff, ])
      dat <- dat[dat$f.stat.exposure >= f_cutoff & !is.na(dat$f.stat.exposure), ]

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
                                              clump = ifelse(is.null(bfile) && is.null(plink_bin),
                                                             T,
                                                             F),
                                              r2 = r2,
                                              kb = kb) %>%
        tidyr::separate(exposure, "exposure", sep = "\\|\\|", extra = "drop") %>%
        dplyr::mutate(id.exposure = id1$info$trait[i])

      dat$exposure <- trimws(dat$exposure)

      if (!is.null(bfile) && !is.null(plink_bin))
      {
        dat <- dplyr::mutate(dat,
                             rsid = SNP,
                             pval = pval.exposure) %>%
          ieugwasr::ld_clump(clump_kb = kb, clump_r2 = r2, pop = pop,
                             bfile = bfile, plink_bin = plink_bin)
      }

      # PVE / F-statistic
      dat$maf.exposure <- ifelse(dat$eaf.exposure < 0.5, dat$eaf.exposure, 1 - dat$eaf.exposure)
      dat$pve.exposure <- mapply(.calc_pve,
                                 dat$beta.exposure,
                                 dat$maf.exposure,
                                 dat$se.exposure,
                                 dat$samplesize.exposure)
      dat$f.stat.exposure <- ((dat$samplesize.exposure - 1 - 1) / 1) * (dat$pve.exposure / (1 - dat$pve.exposure))

      #dat$f.stat.exposure <- dat$beta.exposure ** 2 / dat$se.exposure ** 2
      rem.f.stat <- nrow(dat[dat$f.stat.exposure < f_cutoff, ])
      dat <- dat[dat$f.stat.exposure >= f_cutoff, ]

      # Annotate
      id1$annotate(id1$info$id[i], unique(dat$exposure))

      if (length(dat)) {
        #f_plot(dat, report)

        report$add_dataset(report$get_exposures_name(),
                           c(unique(dat$exposure),
                             id1$info$id[i],
                             id1$info$trait[i],
                             nrow(dat),
                             ifelse(is.null(bfile) && is.null(plink_bin), T, F),
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

DatasetIDs <- setRefClass("DatasetIDs",
                        fields = list(info = "data.frame"),
                        methods = list(
                          initialise = function(ids) {
                            info <<- dplyr::tibble(id = ids, trait = ids, filename = ids)

                            loc <- file.exists(ids)
                            if (any(loc))
                            {
                              info[loc,]$filename <<- file.path(info[loc,]$filename)
                              info[loc,]$trait <<- tools::file_path_sans_ext(basename(info[loc,]$filename))
                              info[loc,]$id <<- tools::file_path_sans_ext(basename(info[loc,]$filename))
                            }
                          },

                          annotate = function(id, name) {
                            name <- fs::path_sanitize(name)

                            # If trait name is ENSG, annotate these using biomaRt to hgnc_symbol
                            if (substring(name, 1, 4) == "ENSG") {
                              hgnc <- .ensg_to_name(name)

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


# Not the most efficient or R-like but since these should be relatively small,
# I think it's fine to do this. tidyr::crossing does something similar but
# renaming the cols is weird
.combine_ids <- function(id1, id2)
{
  df = data.frame(id.x = character(), trait.x = character(), filename.x = character(),
                  id.y = character(), trait.y = character(), filename.y = character())
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

.calc_pve <- function(b, maf, se, n)
{
  tryCatch(
    expr = {
      pve <- (2 * (b^2) * maf * (1 - maf)) /
        ((2 * (b^2) * maf * (1 - maf)) + ((se^2) * 2 * n * maf * (1 - maf)))
    },
    error = function(x) {
      pve <- -1
    }
  )
  return(pve)
}

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

  if (is.na(mart.gene)) {
    return(NULL)
  }

  results <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'),
                   filters = 'ensembl_gene_id',
                   values = ensg,
                   mart = mart.gene)

  return(results)
}
