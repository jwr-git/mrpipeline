do_mr <- function(dat, report)
{
  res <- TwoSampleMR::mr(dat) %>% TwoSampleMR::generate_odds_ratios()

  res %>%
    dplyr::group_by(id.exposure, id.outcome) %>%
    dplyr::group_map(~ scatter_plot(.x, .y, dat, report), .keep = T)

  # Volcano plot of MR results
  res %>%
    dplyr::group_by(id.exposure) %>%
    dplyr::group_map(~ volcano_plot(.x, report), .keep = T)

  # PheWAS plot of MR results
  res %>%
    dplyr::group_by(id.exposure) %>%
    dplyr::group_map(~ phewas_plot(.x, report), .keep = T)

  main <- res[res$method %in% c("Wald ratio", "Inverse variance weighted"),]
  sensitivity <- res[!(res$method %in% c("Wald ratio", "Inverse variance weighted")),]

  report$add_results(main[c("id.exposure", "id.outcome", "exposure", "outcome",
                            "method", "nsnp", "pval",
                            "or", "or_lci95", "or_uci95")])

  report$add_sresults(sensitivity[c("id.exposure", "id.outcome", "exposure", "outcome",
                                    "method", "nsnp", "pval",
                                    "or", "or_lci95", "or_uci95")])

  report$bonferroni <- 0.05 / nrow(main)

  return(res)
}

do_coloc <- function(id1, id2, dat, report,
                     window = 400, # Kb
                     chrpos = "",
                     nthreads = 1)
{
  # Each unique pair of traits need to be colocalised
  pairs <- .combine_ids(id1, id2)

  cres <- parallel::mclapply(1:nrow(pairs), function(i)
  {
    subdat <- dat[dat$id.exposure == pairs[i, "id.x"] & dat$id.outcome == pairs[i, "id.y"],]

    # Select region for which to do coloc
    # atm very simply the lowest P-value region
    if (chrpos == "")
    {
      idx <- which.min(subdat$pval.exposure)
      chrpos <- paste0(subdat$chr[idx], ":", max(as.numeric(subdat$pos[idx]) - window * 100, 0), "-", as.numeric(subdat$pos[idx]) + window * 100)
    }

    f1 <- pairs[i, "filename.x"]
    f2 <- pairs[i, "filename.y"]

    if (file.exists(f1) && file.exists(f2))
    {
      cdat <- gwasglue::gwasvcf_to_coloc(f1, f2, chrpos)
      regional_plot(cdat, report)

      return(coloc_sub(cdat[[1]], cdat[[2]], pairs[i, "id.x"], pairs[i, "id.y"], chrpos, report))
    } else if (!file.exists(f1) && !file.exists(f2))
    {
      cdat <- gwasglue::ieugwasr_to_coloc(pairs[i, "id.x"], pairs[i, "id.y"], chrpos)
      regional_plot(cdat, report)

      return(coloc_sub(cdat[[1]], cdat[[2]], pairs[i, "id.x"], pairs[i, "id.y"], chrpos, report))
    } else {
      # One is file, one is not
      ## TODO me :(
    }
  }, mc.cores = nthreads)

  return(cres)
}

coloc_sub <- function(dat1, dat2, id1, id2, chrpos, report)
{
  # MAFs should be similar and, as some GWAS are missing MAFs,
  # these should be mergeable between the two datasets
  if (any(is.na(dat1$MAF)) || any(is.na(dat2$MAF)))
  {
    dat1$MAF <- dat2$MAF <- dplyr::coalesce(dat1$MAF, dat2$MAF)
  }

  # Even after merge, if any MAF are missing, it's probs
  # best to ignore this coloc
  # TODO - better way of dealing with this?
  if (any(is.na(dat1$MAF)) || any(is.na(dat2$MAF)))
  {
    report$add_cresults(list(id1 = id1,
                             id2 = id2,
                             nsps = 0,
                             H0 = 0,
                             H1 = 0,
                             H2 = 0,
                             H3 = 0,
                             H4 = 0,
                             chrpos = chrpos))
    return(NA)
  }

  cres <- coloc::coloc.abf(dat1, dat2)

  report$add_cresults(list(id1 = id1,
                           id2 = id2,
                           nsnps = cres$summary[[1]],
                           H0 = cres$summary[[2]],
                           H1 = cres$summary[[3]],
                           H2 = cres$summary[[4]],
                           H3 = cres$summary[[5]],
                           H4 = cres$summary[[6]],
                           chrpos = chrpos))
  return(cres)
}
