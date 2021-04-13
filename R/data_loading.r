read_datasets <- function(id1, id2, report, p_cutoff, chrompos_query, gene_query,
                          f_cutoff = 10,
                          r2 = 0.001, kb = 10000, pop = "EUR",
                          proxies = T,
                          action = 1,
                          nthreads = 1)
{
  exposure_dat <- parallel::mclapply(1:nrow(id1), function(i)
  {
    f <- id1$filename[i]

    # Load locally from vcf
    if (file.exists(f))
    {
      dat <- gwasvcf::query_gwas(vcf = VariantAnnotation::readVcf(f), pval = p_cutoff) %>%
        gwasglue::gwasvcf_to_TwoSampleMR("exposure") %>%
        dplyr::mutate(exposure = id1$trait[i], trait.exposure = id1$trait[i], id.exposure = id1$trait[i])

      dat$f.stat.exposure <- dat$beta.exposure ** 2 / dat$se.exposure ** 2
      rem.f.stat <- nrow(dat[dat$f.stat.exposure < f_cutoff, ])
      dat <- dat[dat$f.stat.exposure >= f_cutoff, ]

      if (length(dat)) {
        pre_clump <- nrow(dat)

        dat %<>% TwoSampleMR::clump_data(clump_kb = kb, clump_r2 = r2, pop = pop)

        pipeline_f_plot(dat, report)

        report$add_dataset(report$get_exposures_name(), c(unique(dat$exposure), id1$id[i], pre_clump, F, rem.f.stat, pre_clump - nrow(dat), nrow(dat)))
      } else {
        dat <- NULL

        report$add_dataset(report$get_exposures_name(), c(unique(dat$exposure), id1$id[i], 0, F))
      }
    } else
    {
      # From OpenGWAS
      dat <- TwoSampleMR::extract_instruments(id1$id[i], p1 = p_cutoff, clump = T, r2 = r2, kb = kb) %>%
        tidyr::separate(exposure, "exposure", sep = "\\|\\|", extra = "drop") %>%
        dplyr::mutate(id.exposure = id1$trait[i])

      dat$f.stat.exposure <- dat$beta.exposure ** 2 / dat$se.exposure ** 2
      rem.f.stat <- nrow(dat[dat$f.stat.exposure < as.numeric(f_cutoff), ])
      dat <- dat[dat$f.stat.exposure >= as.numeric(f_cutoff), ]

      if (length(dat)) {
        pipeline_f_plot(dat, report)

        report$add_dataset(report$get_exposures_name(), c(unique(dat$exposure), id1$id[i], nrow(dat), T, rem.f.stat, 0, nrow(dat)))
      } else {
        report$add_dataset(report$get_exposures_name(), c(unique(dat$exposure), id1$id[i], 0, F, 0, 0, 0))
      }
    }

    return (dat)
  }, mc.cores = nthreads) %>%
    dplyr::bind_rows()

  rsids <- unique(exposure_dat$SNP)

  outcome_dat <- parallel::mclapply(1:nrow(id2), function(i)
  {
    f <- id2$filename[i]

    # Is local vcf file? Load from there
    if (file.exists(f))
    {
      dat <- gwasvcf::query_gwas(vcf = VariantAnnotation::readVcf(f), rsid = rsids) %>%
        gwasglue::gwasvcf_to_TwoSampleMR("outcome") %>%
        dplyr::mutate(id.outcome = id2$trait[i])

      if (length(dat)) {
        report$add_dataset(report$get_outcomes_name(), c(unique(dat$outcome), id2$id[i], nrow(dat), 0, nrow(dat)))
      } else {
        dat <- NULL

        report$add_dataset(report$get_outcomes_name(), c(unique(dat$outcome), id2$id[i], 0, 0, 0))
      }
    } else
    {
      # From OpenGWAS
      dat <- TwoSampleMR::extract_outcome_data(rsids, id2$id[i], proxies = proxies) %>%
        tidyr::separate(outcome, "outcome", sep = "\\|\\|", extra = "drop") %>%
        dplyr::mutate(id.outcome = id2$trait[i])

      if (length(dat)) {
        report$add_dataset(report$get_outcomes_name(), c(unique(dat$outcome), id2$id[i], nrow(dat), nrow(dat[!is.na(dat$proxy.outcome) & dat$proxy.outcome == T,]), nrow(dat)))
      } else {
        report$add_dataset(report$get_exposures_name(), c(unique(dat$outcome), id2$id[i], 0, F, 0, 0, 0))
      }
    }

    return(dat)
  }, mc.cores = nthreads) %>%
    dplyr::bind_rows()

  dat <- TwoSampleMR::harmonise_data(exposure_dat, outcome_dat, action = action)

  return(dat)
}

# Not the most efficient or R-like but since these should be relatively small,
# I think it's fine to do this. tidyr::crossing does something similar but
# renaming the cols is weird
.combine_ids <- function(id1, id2)
{
  df = data.frame(id.x = character(), trait.x = character(), filename.x = character(),
                  id.y = character(), trait.y = character(), filename.y = character())
  for (r1 in 1:nrow(id1))
  {
    for (r2 in 1:nrow(id2))
    {
      df[nrow(df) + 1, ] <- c(id1[r1, ], id2[r2, ])
    }
  }
  #return(df %<>% dplyr::rename())
  return(df)
}

