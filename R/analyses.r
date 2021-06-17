#' Calculates F-statistic and related
#'
#' @param dat A data.frame of data
#'
#' @return List of data.frame and integer
.calc_f_stat <- function(dat, f_cutoff)
{
  full.f.stat <- T

  if (any(is.na(dat$eaf.exposure))) {
    warning("Using approximate F-statistic as some allele frequencies are missing.")
    full.f.stat <- F
  }
  else if (any(is.na(dat$samplesize.exposure))) {
    warning("Using approximate F-statistic as some sample sizes are missing.")
    full.f.stat <- F
  }
  else if (length(unique(dat$SNP)) == 1) {
    full.f.stat <- F
  }

  # PVE / F-statistic
  if (full.f.stat) {
    dat$maf.exposure <- ifelse(dat$eaf.exposure < 0.5, dat$eaf.exposure, 1 - dat$eaf.exposure)
    dat$pve.exposure <- mapply(.calc_pve,
                               dat$beta.exposure,
                               dat$maf.exposure,
                               dat$se.exposure,
                               dat$samplesize.exposure)
    # From https://doi.org/10.1093/ije/dyr036
    dat$f.stat.exposure <- ifelse(dat$pve.exposure == -1, 0,
                                  ((dat$samplesize.exposure - length(unique(dat$SNP)) - 1) / length(unique(dat$SNP))) * (dat$pve.exposure / (1 - dat$pve.exposure)))
  } else {
    message("Using approximate F-statistic as no allele frequency has been provided.")
    dat$f.stat.exposure <- dat$beta.exposure ** 2 / dat$se.exposure ** 2
  }

  rem.f.stat <- nrow(dat[dat$f.stat.exposure < f_cutoff, ])
  dat <- dat[dat$f.stat.exposure >= f_cutoff & !is.na(dat$f.stat.exposure), ]

  return(list(dat, rem.f.stat))
}

#' Calculates proportion of variance explained
#' From https://doi.org/10.1371/journal.pone.0120758 S1 Text
#'
#' @param b Vector or number, beta
#' @param maf Vector or number, minor allele frequency
#' @param se Vector or number, standard error of beta
#' @param n Vector or number, sample size
#'
#' @return Vector or number, proportion of variance explained
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

#' Calculates the second term Taylor approximation for standard error of
#' the Wald ratio method.
#' From https://doi.org/10.1101/2021.03.01.433439 supplementary
#'
#' @param object Data.frame or vector
#'
#' @return Appended result
.wr_taylor_approx <- function(dat)
{
  b <- dat$beta.outcome / dat$beta.exposure
  se <- (dat$se.outcome ** 2 / dat$beta.exposure ** 2) + ((dat$beta.outcome ** 2 * dat$se.exposure ** 2) / (dat$beta.exposure ** 4))
  se <- sqrt(se)
  pval <- pnorm(abs(b) / se, lower.tail = F) * 2

  res <- data.frame(
    id.exposure = dat$id.exposure[1],
    id.outcome = dat$id.outcome[1],
    outcome = dat$outcome[1],
    exposure = dat$exposure[1],
    method = "Wald ratio",
    nsnp = 1,
    b = b,
    se = se,
    pval = pval
  )

  return(res)
}

#' Runs Mendelian randomisation and related analyses, including heterogeneity,
#' Steiger filtering and pleiotropy analyses.
#' Also generates a number of plots, e.g. volcano plots, for MR results.
#'
#' @param dat A data.frame of harmonised data
#' @param report QCReport class of results, etc. for reporting
#' @param conf config::config file of parameter
#'
#' @return A data.frame of MR results
do_mr <- function(dat, report, conf)
{
  res <- plyr::ddply(dat, c("id.exposure", "id.outcome"), function(x1)
  {
    x <- subset(x1, mr_keep.exposure)

    if (nrow(x) == 1 & conf$wr_taylor_approx == T) {
      res <- .wr_taylor_approx(x)
    }
    else {
      res <- TwoSampleMR::mr(x)
    }

    res <- subset(res, !(is.na(b) & is.na(se) & is.na(pval)))
  }) %>%
    TwoSampleMR::generate_odds_ratios()

  #res %>%
  #  dplyr::group_by(id.exposure, id.outcome) %>%
  #  dplyr::group_map(~ mr_scatter_plot(.x, .y, dat, report), .keep = T)

  # Volcano plot of each SNPs' WR results
  TwoSampleMR::mr_singlesnp(dat) %>%
    TwoSampleMR::generate_odds_ratios() %>%
    dplyr::filter(!startsWith(getElement(., "SNP"), "All")) %>%
    dplyr::group_by(id.exposure) %>%
    dplyr::group_map(~ volcano_plot(.x, report), .keep = T)

  # PheWAS plot of MR results
  res %>%
    dplyr::group_by(id.exposure) %>%
    dplyr::group_map(~ phewas_plot(.x, report), .keep = T)

  # Report results separately
  main <- res[res$method %in% c("Wald ratio", "Inverse variance weighted"),]
  sensitivity <- res[!(res$method %in% c("Wald ratio", "Inverse variance weighted")),]

  report$add_results(main[c("id.exposure", "id.outcome", "exposure", "outcome",
                            "method", "nsnp", "pval",
                            "or", "or_lci95", "or_uci95")])

  report$add_sresults("sresults", sensitivity[c("id.exposure", "id.outcome", "exposure", "outcome",
                                    "method", "nsnp", "pval",
                                    "or", "or_lci95", "or_uci95")])

  # Heterogeneity
  hetres <- TwoSampleMR::mr_heterogeneity(dat)
  report$add_sresults("hetresults", hetres)

  # Pleiotropy
  pleiores <- TwoSampleMR::mr_pleiotropy_test(dat)
  report$add_sresults("pleioresults", pleiores)

  # Steiger
  steigerres <- TwoSampleMR::directionality_test(dat)
  steigerres$flag <- ifelse(steigerres$correct_causal_direction == T & steigerres$steiger_pval < 0.05,
                            "True",
                            ifelse(steigerres$correct_causal_direction == F & steigerres$steiger_pval < 0.05,
                                   "False",
                                   "Unknown"))
  report$add_sresults("steigerresults", steigerres)

  report$bonferroni <- 0.05 / nrow(main)

  return(res)
}

#' Runs colocalisation using the coloc R package's coloc.abf function
#'
#' @param dat A data.frame of harmonised data
#' @param id1 DatasetsID class of exposure IDs
#' @param id2 DatasetsID class of outcome IDs
#' @param report QCReport class of results, etc. for reporting
#' @param conf config::config file of parameter
#'
#' @return A data.frame of colocalistion results
do_coloc <- function(dat, id1, id2, report, conf)
{
  # Each unique pair of traits need to be colocalised
  pairs <- .combine_ids(id1, id2)

  cres <- parallel::mclapply(1:nrow(pairs), function(i)
  {
    subdat <- dat[dat$id.exposure == pairs[i, "id.x"] & dat$id.outcome == pairs[i, "id.y"],]

    if (length(subdat) < 1 || nrow(subdat) < 1) {
      return(NULL)
    }

    # Select region for which to do coloc
    # atm very simply the lowest P-value region
    idx <- which.min(subdat$pval.exposure)
    chrpos <- paste0(subdat$chr[idx],
                     ":",
                     max(as.numeric(subdat$pos[idx]) - conf$coloc_window * 100, 0),
                     "-",
                     as.numeric(subdat$pos[idx]) + conf$coloc_window * 100)

    f1 <- pairs[i, "filename.x"]
    f2 <- pairs[i, "filename.y"]

    if (file.exists(f1) & file.exists(f2))
    {
      cdat <- gwasglue::gwasvcf_to_coloc(f1, f2, chrpos)
      regional_plot(cdat, pairs, i, report)

      return(coloc_sub(cdat[[1]], cdat[[2]], pairs, i, chrpos, report, conf))
    } else if (!file.exists(f1) & !file.exists(f2))
    {
      cdat <- gwasglue::ieugwasr_to_coloc(pairs[i, "id.x"], pairs[i, "id.y"], chrpos)
      regional_plot(cdat, pairs, i, report)

      return(coloc_sub(cdat[[1]], cdat[[2]], pairs, i, chrpos, report, conf))
    } else {
      # One is file, one is not
      ## TODO me :(
    }
  }, mc.cores = conf$cores)
  return(cres)
}

#' Sub-function for the colocalisation analyses
#'
#' @param dat1 SNPs, etc. from first dataset
#' @param dat2 SNPs, etc. from second dataset
#' @param pairs A data.frame of pairwise combined IDs
#' @param i Location in pairs being colocalised
#' @param chrpos Chromosome position to search for SNPs, default uses region
#'               around lead SNP but can be overwritten if provided
#' @param report QCReport class of results, etc. for reporting
#' @param conf config::config file of parameter
#'
#' @return Results, or empty if cannot run the colocalisation
coloc_sub <- function(dat1, dat2, pairs, i, chrpos, report, conf)
{
  if (!length(dat1) || !length(dat2))
  {
    report$add_cresults(list(id1 = pairs[i, "id.x"],
                             id2 = pairs[i, "id.y"],
                             name1 = pairs[i, "trait.x"],
                             name2 = pairs[i, "trait.y"],
                             nsnps = 0,
                             H0 = 0,
                             H1 = 0,
                             H2 = 0,
                             H3 = 0,
                             H4 = 0,
                             chrpos = chrpos,
                             plot = ""))
    return()
  }

  # MAFs should be similar and, as some GWAS are missing MAFs,
  # these should be mergeable between the two datasets
  if (any(is.na(dat1$MAF)) || any(is.na(dat2$MAF)))
  {
    dat1$MAF <- dat2$MAF <- dplyr::coalesce(dat1$MAF, dat2$MAF)
  }

  # Even after merge, if any MAF are missing, it's probs
  # best to ignore this coloc
  # TODO - better way of dealing with this?
  if (any(is.na(dat1$MAF)) || any(is.na(dat2$MAF)) ||
      length(dat1$snp) < conf$coloc_min_snps || length(dat2$snp) < conf$coloc_min_snps)
  {
    report$add_cresults(list(id1 = pairs[i, "id.x"],
                             id2 = pairs[i, "id.y"],
                             name1 = pairs[i, "trait.x"],
                             name2 = pairs[i, "trait.y"],
                             nsnps = length(dat1$snp),
                             H0 = -1,
                             H1 = -1,
                             H2 = -1,
                             H3 = -1,
                             H4 = -1,
                             chrpos = chrpos,
                             plot = ""))
    return()
  }

  cres <- coloc::coloc.abf(dat1, dat2,
                           p1 = conf$coloc_p1,
                           p2 = conf$coloc_p2,
                           p12 = conf$coloc_p12)

  report$add_cresults(list(id1 = pairs[i, "id.x"],
                           id2 = pairs[i, "id.y"],
                           name1 = pairs[i, "trait.x"],
                           name2 = pairs[i, "trait.y"],
                           nsnps = cres$summary[[1]],
                           H0 = cres$summary[[2]],
                           H1 = cres$summary[[3]],
                           H2 = cres$summary[[4]],
                           H3 = cres$summary[[5]],
                           H4 = cres$summary[[6]],
                           chrpos = chrpos,
                           plot = ""))
  return(cres)
}

#' Sub-function to run heterogeneity analyses
#' UNIMPLEMENTED
do_heterogeneity <- function(dat, report)
{
  single <- TwoSampleMR::mr_singlesnp(dat)

  for (exposure in single$exposure)
  {
    sres <- single[single$exposure == exposure, ]
    sres <- sres[!startsWith(sres$SNP, "All"), ] # Exclude all analyses

    if (length(sres) == 0 | nrow(sres) < 2)
    {
      report$add_hresults(list(id.exposure = sres$id.exposure[1],
                               id.outcome = sres$id.outcome[1],
                               exposure = sres$exposure[1],
                               outcome = sres$outcome[1],
                               method = "",
                               Q = 0,
                               Q_df = 0,
                               Q_pval = 0))
      next
    }

    # Multiplicative random effects
    res <- summary(lm(b_out ~ -1 + b_exp, weights = 1/se_out^2))
    Q_df = length(b_exp) - 1
    Q = res$sigma^2 * Q_df
    Q_pval = pchisq(Q, Q_df, low=FALSE)
  }

  #res_nosens <- res[res$method %in% c("Wald ratio", "Inverse variance weighted"), ]

  # Heterogeneity if more than 1 SNP
  #if (length(dat) && nrow(dat) > 1) {
  #  hetero <- TwoSampleMR::mr_ivw(res_nosens$beta.exposure,
  #                                res_nosens$beta.outcome,
  #                                res_nosens$se.exposure,
  #                                res_nosens$se.outcome)
  #}
}
