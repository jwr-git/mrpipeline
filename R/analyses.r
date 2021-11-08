#' Calculates F-statistic and related
#'
#' @param dat A data.frame of data
#'
#' @return List of data.frame and integer
calc_f_stat <- function(dat, f_cutoff = 10, verbose = TRUE)
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
    .print_msg("Using approximate F-statistic as no allele frequency has been provided.", verbose)
    dat$f.stat.exposure <- dat$beta.exposure ** 2 / dat$se.exposure ** 2
  }

  rem.f.stat <- nrow(dat[dat$f.stat.exposure < f_cutoff, ])
  .print_msg(paste0("Amount of SNPs which did not meet threshold: ", rem.f.stat), verbose)

  #dat <- dat[dat$f.stat.exposure >= f_cutoff & !is.na(dat$f.stat.exposure), ]

  return(dat)
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
    snp = unique(dat$SNP),
    b = b,
    se = se,
    pval = pval
  )

  return(res)
}

.ivw_delta <- function(dat)
{
  nsnps <- nrow(dat)
  psi <- 0

  summary <- summary(lm(dat$beta.outcome ~ dat$beta.exposure - 1, weights = (dat$se.outcome^2 + dat$beta.outcome^2*dat$se.exposure^2/dat$beta.exposure^2-2*psi*dat$beta.outcome*dat$se.exposure*dat$se.outcome/dat$beta.exposure)^-1))
  IVWbeta <- summary$coef[1]
  IVWse <- summary$coef[1,2] / min(summary$sigma, 1)
  pval <- 2 * pnorm(-abs(IVWbeta / IVWse))

  #omega <- sqrt(dat$se.outcome ** 2 + dat$beta.outcome ** 2 * dat$se.exposure ** 2 / dat$beta.exposure ** 2) %o%
  #  sqrt(dat$se.outcome ** 2 + dat$beta.outcome ** 2 * dat$se.exposure ** 2 / dat$beta.exposure ** 2)

  #omega_s <- solve(omega)

  #IVWbeta <- as.numeric(solve(t(dat$beta.exposure) %*% omega_s %*% dat$beta.exposure)
  #                      * t(dat$beta.exposure) %*% omega_s %*% dat$beta.outcome)

  # Fixed effect error
  #IVWse <- sqrt(solve(t(dat$beta.exposure) %*% omega_s %*% dat$beta.exposure))

  # Random effect error
  #rse <- dat$beta.outcome - IVWbeta * dat$beta.exposure
  #IVWse <- sqrt(solve(t(dat$beta.exposure) %*% omega_s %*% dat$beta.exposure)) *
  #  max(sqrt(t(rse) %*% omega_s %*% rse / (nsnps - 1)), 1)

  #pval <- 2 * pnorm(-abs(IVWbeta / IVWse))

  res <- data.frame(
    id.exposure = dat$id.exposure[1],
    id.outcome = dat$id.outcome[1],
    outcome = dat$outcome[1],
    exposure = dat$exposure[1],
    method = "Inverse variance weighted",
    nsnp = nsnps,
    snp = paste(dat$SNP, collapse = ", "),
    b = IVWbeta,
    se = IVWse,
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
do_mr <- function(dat, f_cutoff = 10, verbose = TRUE)
{
  if (!is.null(f_cutoff)) {
    if ("f.stat.exposure" %in% names(dat)) {
      dat <- dat[dat$f.stat.exposure >= f_cutoff & !is.na(dat$f.stat.exposure), ]
    } else {
      .print_msg("F-statistic cut-off given but could not find column \"f.stat.exposure\". Calculating these now.", verbose = verbose)
      dat <- calc_f_stat(dat, f_cutoff = f_cutoff, verbose = verbose)
    }
  }

  res <- plyr::ddply(dat, c("id.exposure", "id.outcome"), function(x1)
  {
    x <- subset(x1, mr_keep.exposure)

    nsnps <- nrow(x)

    # WR on all single SNPs
    res <- lapply(1:nsnps, function(i)
    {
      with(x, .wr_taylor_approx(x[i, ]))
    })
    res <- do.call(rbind, res)

    # IVW if applicable
    if (nsnps > 1)
    {
      tryCatch(
        expr = {
          res_ <- .ivw_delta(x)
          res <- rbind(res, res_)
        },
        error = function(e) {
          message("Error encounted in IVW, no result for this will be given: ")
          message(e)
        }
      )
    }

    if (nrow(res) > 0) {
      res <- subset(res, !(is.na(b) & is.na(se) & is.na(pval)))
    }
  }) %>%
    TwoSampleMR::generate_odds_ratios()

  # Steiger
  if (any(is.na(dat$samplesize.exposure)) || any(is.na(dat$samplesize.outcome))) {
    warning("Samplesizes are required for Steiger filtering.")
  } else {
    steigerres <- TwoSampleMR::directionality_test(dat)
    steigerres$steigerflag <- ifelse(steigerres$correct_causal_direction == T & steigerres$steiger_pval < 0.05,
                              "True",
                              ifelse(steigerres$correct_causal_direction == F & steigerres$steiger_pval < 0.05,
                                     "False",
                                     "Unknown"))

    res <- base::merge(res, steigerres, by = c("exposure", "outcome", "id.exposure", "id.outcome"), all.x = TRUE)
  }

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
do_coloc <- function(dat,
                     method = "coloc.abf",
                     coloc_window = 500000,
                     bfile = NULL,
                     pwcoco = NULL,
                     workdir = tempdir(),
                     cores = 1,
                     verbose = TRUE)
{
  if (!all(c("file.exposure", "file.outcome") %in% names(dat))) {
    warning("Requires harmonised dataset for colocalisation; \"file.exposure\" and \"file.outcome\" columns not found.")
    return(NULL)
  }

  if (!(method %in% c("coloc.abf", "pwcoco", "coloc.susie"))) {
    warning("Colocalisation method not understood, defaulting to \"coloc.abf\". Must be one of the following: coloc.abf, pwcoco, coloc.susie.")
    method = "coloc.abf"
  }

  if (method == "pwcoco" && Sys.info()['sysname'] == "Windows") {
    .print_msg("PWCoCo is not an allowed method when using this pipeline in Windows, defaulting to \"coloc.abf\".", verbose = verbose)
    method = "coloc.abf"
  }
  else if (method == "pwcoco" && (is.null(bfile) || is.null(pwcoco))) {
    warning("PWCoCo requires the bfile and pwcoco arguments for paths to the Plink reference data and PWCoCo executible, respectively.")
    return(NULL)
  }

  # Each unique pair of traits need to be colocalised
  pairs <- tidyr::crossing(dat$file.exposure, dat$file.outcome)
  names(pairs) <- c("file.exposure", "file.outcome")

  coloc_res <- parallel::mclapply(1:nrow(pairs), function(i)
  {
    subdat <- dat[which(dat$file.exposure == pairs[i, "file.exposure"][[1]] & dat$file.outcome == pairs[i, "file.outcome"][[1]]), ]

    if (length(subdat) < 1 || nrow(subdat) < 1) {
      return(NULL)
    }

    # Select region for which to do coloc
    # atm very simply the lowest P-value region
    idx <- which.min(subdat$pval.exposure)
    chrpos <- paste0(subdat$chr.exposure[idx],
                     ":",
                     max(as.numeric(subdat$position.exposure[idx]) - coloc_window, 0),
                     "-",
                     as.numeric(subdat$position.exposure[idx]) + coloc_window)

    f1 <- pairs[i, "file.exposure"][[1]]
    f2 <- pairs[i, "file.outcome"][[1]]

    if (method == "coloc.abf") {
      if (file.exists(f1) && file.exists(f2))
      {
        cdat <- gwasglue::gwasvcf_to_coloc(f1, f2, chrpos)
      } else if (!file.exists(f1) && !file.exists(f2))
      {
        cdat <- gwasglue::ieugwasr_to_coloc(f1, f2, chrpos)
      } else {
        # One is file, one is not
        ## TODO me :(
      }

      cres <- .coloc_sub(cdat[[1]], cdat[[2]], verbose = verbose)
    } else if (method == "pwcoco") {
      if (file.exists(f1) && file.exists(f2)) {
        .gwasvcf_to_pwcoco(f1, f2, chrpos, outfile = workdir)
      } else if (!file.exists(f1) && !file.exists(f2)) {
        .ieugwasr_to_pwcoco(f1, f2, chrpos, outfile = workdir)
      }

      cres <- .pwcoco_sub(chrpos, bfile, pwcoco, workdir)
    }



    if (method == "coloc.abf") {
      cres <- .coloc_sub(cdat[[1]], cdat[[2]], verbose = verbose)
    } else if (method == "pwcoco") {
      cres <- .pwcoco_sub(cdat[[1]], cdat[[2]], verbose = verbose)
    }

    #regional_plot(cdat, pairs, i, report)

    if (length(cres)) {
      cres[length(cres) + 1] <- f1
      names(cres)[length(cres)] <- "file.exposure"
      cres[length(cres) + 1] <- f2
      names(cres)[length(cres)] <- "file.outcome"
    }

    return(cres)
  }, mc.cores = cores) %>%
    dplyr::bind_rows()

  return(coloc_res)
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
.coloc_sub <- function(dat1, dat2,
                       min_snps = 100,
                       p1 = 1e-4,
                       p2 = 1e-4,
                       p12 = 1e-5,
                       verbose = TRUE)
{
  if (!length(dat1) || !length(dat2)) {
    .print_msg("No SNPs extracted for colocalisation analysis.", verbose)
    return(NULL)
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
      length(dat1$snp) < min_snps || length(dat2$snp) < min_snps)
  {
    .print_msg("Minor allele frequenices are required for coloc but were not found in these datasets.", verbose)
    return(NULL)
  }

  cres <- coloc::coloc.abf(dat1, dat2,
                           p1 = p1,
                           p2 = p2,
                           p12 = p12)

  return(cres$summary)
}

#' Write files for PWCoCo where data are read from two VCF objects or files.
#'
#' @param vcf1 VCF object or path to VCF file
#' @param vcf2 VCF object or path to VCF file
#' @param chrompos Character of the format chr:pos1-pos2
#' @param type1 How to treat vcffile1 for coloc, either "quant" or "cc"
#' @param type2 How to treat vcffile2 for coloc, either "quant" or "cc"
#' @param outfile Path to output files, without file ending
#'
#' return 0 if success, 1 if there was a problem
.gwasvcf_to_pwcoco <- function(vcf1, vcf2, chrompos, type1=NULL, type2=NULL, outfile)
{
  overlap <- gwasvcf::vcflist_overlaps(list(vcf1, vcf2), chrompos)
  vcf1 <- overlap[[1]]
  vcf2 <- overlap[[2]]

  if (length(vcf1) == 0 || length(vcf2) == 0)
  {
    message("No overlaps for the given chrompos in ", ifelse(length(vcf1) == 0, "vcf1", "vcf2"), ".")
    return(1)
  }

  # vcf1
  tib1 <- vcf1 %>% gwasvcf::vcf_to_granges() %>% dplyr::as_tibble() %>%
    dplyr::select(rsid, ALT, REF, AF, ES, SE, LP, SS, NC) %>%
    dplyr::rename(
      SNP = rsid,
      A1 = ALT,
      A2 = REF,
      freq = AF,
      b = ES,
      se = SE,
      p = LP,
      N = ss,
      N_case = NC
    )
  tib1$p <- 10^(-tib1$p)

  # Coloc type -- if study type is continuous then do not need the case column
  if (type1 == "quant" || VariantAnnotation::header(vcf1) %>% VariantAnnotation::meta() %>% {.[["SAMPLE"]][["StudyType"]]} == "Continuous")
  {
    tib1 <- tib1[c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")]
  }

  # vcf2
  tib2 <- vcf2 %>% gwasvcf::vcf_to_granges() %>% dplyr::as_tibble() %>%
    dplyr::select(rsid, ALT, REF, AF, ES, SE, LP, SS, NC) %>%
    dplyr::rename(
      SNP = rsid,
      A1 = ALT,
      A2 = REF,
      freq = AF,
      b = ES,
      se = SE,
      p = LP,
      N = ss,
      N_case = NC
    )
  tib2$p <- 10^(-tib2$p)

  # Coloc type -- if study type is continuous then do not need the case column
  if (type2 == "quant" || VariantAnnotation::header(vcf2) %>% VariantAnnotation::meta() %>% {.[["SAMPLE"]][["StudyType"]]} == "Continuous")
  {
    tib2 <- tib2[c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")]
  }

  write.table(tib1, file=paste0(outfile, "1.txt"), row=F, col=T, qu=F)
  write.table(tib1, file=paste0(outfile, "2.txt"), row=F, col=T, qu=F)
  return(0)
}

#' Write files for PWCoCo where data are read from the OpenGWAS DB.
#'
#' @param id1 ID for trait 1
#' @param id2 ID for trait 2
#' @param chrompos Character of the format chr:pos1-pos2
#' @param type1 How to treat vcffile1 for coloc, either "quant" or "cc"
#' @param type2 How to treat vcffile2 for coloc, either "quant" or "cc"
#' @param outfile Path to output files, without file ending
#'
#' return 0 if success, 1 if there was a problem
.ieugwasr_to_pwcoco <- function(id1, id2, chrompos, type1=NULL, type2=NULL, outfile)
{
  tib1 <- ieugwasr::associations(id=id1, variants=chrompos) %>% subset(., !duplicated(rsid))
  tib2 <- ieugwasr::associations(id=id2, variants=chrompos) %>% subset(., !duplicated(rsid))

  if (length(tib1) < 1 || length(tib2) < 1)
  {
    message("Data could not be read using the ieugwasr package for id1 = ", id1, " and id2 = ", id2, ".")
    return(1)
  }

  # Matching the files is quicker for PWCoCo, so best to off-load to that?
  # Save data -- PWCoCo handles the matching and cleaning mostly by itself
  tib1 %<>% dplyr::select(rsid, ea, nea, eaf, beta, se, p, n) %>%
    dplyr::rename(
      SNP = rsid,
      A1 = ea,
      A2 = nea,
      freq = eaf,
      b = beta,
      se = se,
      p = p,
      N = n
    )
  # Need to determine whether there are cases
  info1 <- ieugwasr::gwasinfo(id1)
  if ("ncase" %in% colnames(info1))
  {
    tib1$N_case <- info1$ncase
  }

  tib2 %<>% dplyr::select(rsid, ea, nea, eaf, beta, se, p, n) %>%
    dplyr::rename(
      SNP = rsid,
      A1 = ea,
      A2 = nea,
      freq = eaf,
      b = beta,
      se = se,
      p = p,
      N = n
    )
  info2 <- ieugwasr::gwasinfo(id2)
  if ("ncase" %in% colnames(info2))
  {
    tib2$N_case <- info2$ncase
  }

  write.table(tib1, file=paste0(outfile, "1.txt"), row=F, col=T, qu=F)
  write.table(tib2, file=paste0(outfile, "2.txt"), row=F, col=T, qu=F)
  return(0)
}

.pwcoco_sub <- function(bfile,
                        chrpos,
                        pwcoco,
                        maf = 0.01,
                        p1 = 1e-4,
                        p2 = 1e-4,
                        p12 = 1e-5,
                        workdir = tempdir(),
                        verbose = TRUE)
{
  chr <- as.integer(strsplit(chrpos, ":")[[1]][1])

  cmd <- glue::glue("{pwcoco} --bfile {bfile} --sum_stats1 {file.path(workdir, '1.txt')} --sum_stats2 {file.path(workdir, '2.txt')} --out {file.path(workdir, 'out')} --chr {chr} --maf {maf}")
  system(cmd)

  res <- data.table::fread(file.path(workdir, 'out.coloc'))

  return(res)
}
