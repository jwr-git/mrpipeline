#' Run colocalisation analyses
#'
#' Runs colocalisation using any of the following methods:
#' \enumerate{
#' \item Coloc.abf, see coloc::coloc.abf()
#' \item Coloc.susie, see coloc::coloc.susie()
#' \item PWCoCo, see \href{https://github.com/jwr-git/pwcoco}{PWCoCo}
#' }
#' NB: PWCoCo is not available on Windows.
#'
#' @seealso [coloc::coloc.abf()], [coloc::coloc.susie()]
#'
#' @param dat A data.frame of harmonised data
#' @param method Which method of colocalisation to use: coloc.abf, coloc.susie, pwcoco (Optional)
#' @param coloc_window Size (+/-) of region to extract for colocalisation analyses (Optional)
#' @param plot_region Whether to plot the regions or not
#' @param bfile Path to Plink bed/bim/fam files (Optional)
#' @param plink Path to Plink binary (Optional)
#' @param pwcoco If PWCoCo is the selected coloc method, path to PWCoCo binary (Optional)
#' @param workdir Path to save temporary files (Optional)
#' @param cores Number of cores for multi-threaded tasks (Optional)
#'              NB: Unavailable on Windows machines
#' @param verbose Display verbose information (Optional, boolean)
#'
#' @return A data.frame of colocalistion results
#' @export
#' @importFrom tidyr crossing
#' @importFrom parallel mclapply
#' @importFrom gwasglue ieugwasr_to_coloc
#' @importFrom tibble tibble
do_coloc <- function(dat,
                     method = "coloc.abf",
                     coloc_window = 500000,
                     plot_region = F,
                     bfile = NULL,
                     plink = NULL,
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

    if ("pos.exposure" %in% names(subdat)) {
      pos_var <- "pos.exposure"
    } else if ("position.exposure" %in% names(subdat)) {
      pos_var <- "position.exposure"
    } else {
      pos_var <- "bp.exposure"
    }

    # Select region for which to do coloc
    # atm very simply the lowest P-value region
    idx <- which.min(subdat$pval.exposure)
    chrpos <- paste0(as.character(subdat$chr.exposure[idx]),
                     ":",
                     max(as.numeric(subdat[idx, pos_var]) - coloc_window, 0),
                     "-",
                     as.numeric(subdat[idx, pos_var]) + coloc_window)

    f1 <- pairs[i, "file.exposure"][[1]]
    f2 <- pairs[i, "file.outcome"][[1]]

    if (method == "coloc.abf" || method == "coloc.susie") {
      if (file.exists(f1) && file.exists(f2))
      {
        cdat <- .gwasvcf_to_coloc_rsid(f1, f2, chrpos)
      } else if (!file.exists(f1) && !file.exists(f2))
      {
        cdat <- gwasglue::ieugwasr_to_coloc(f1, f2, chrpos)
      } else {
        # One is file, one is not
        cdat <- .cdat_from_mixed(f1, f2, chrpos, verbose)
      }

      if (all(is.na(cdat))) {
        return(list(data.frame = list(file.exposure = f1,
                                      file.outcome = f2,
                                      nsnps = NA,
                                      PP.H0.abf = NA,
                                      PP.H1.abf = NA,
                                      PP.H2.abf = NA,
                                      PP.H3.abf = NA,
                                      PP.H4.abf = NA),
                    regional = NA, zscore = NA))
      }

      cres <- .coloc_sub(cdat[[1]], cdat[[2]],
                         susie = ifelse(method == "coloc.abf", FALSE, TRUE),
                         verbose = verbose)
    } else if (method == "pwcoco") {
      if (file.exists(f1) && file.exists(f2)) {
        .gwasvcf_to_pwcoco(f1, f2, chrpos, outfile = workdir)
      } else if (!file.exists(f1) && !file.exists(f2)) {
        .ieugwasr_to_pwcoco(f1, f2, chrpos, outfile = workdir)
      }

      cres <- .pwcoco_sub(chrpos, bfile, pwcoco, workdir)
    }

    if (plot_region) {
      if (length(cdat[[1]]$snp) > 500 && (is.null(plink) || is.null(bfile))) {
        warning("Cannot generate regional plot as the number of SNPs in the region is above 500 and no bfile/Plink arguments have been given.")
        rp <- NA
      } else if (length(cdat[[1]]$snp) <= 500) {
        rp <- regional_plot(cdat, subdat$exposure[1], subdat$outcome[1], verbose = verbose)
      } else {
        rp <- regional_plot(cdat, subdat$exposure[1], subdat$outcome[1], bfile = bfile, plink = plink, verbose = verbose)
      }
    } else {
      rp <- NA
    }

    zp <- z_comparison_plot(cdat[[1]], cdat[[2]], verbose = verbose)

    if (all(is.na(cres))) {
      return(list(data.frame = list(file.exposure = f1,
                                    file.outcome = f2,
                                    nsnps = NA,
                                    PP.H0.abf = NA,
                                    PP.H1.abf = NA,
                                    PP.H2.abf = NA,
                                    PP.H3.abf = NA,
                                    PP.H4.abf = NA),
                  regional = NA, zscore = NA))
    }

    if (length(cres)) {
      cres[length(cres) + 1] <- f1
      names(cres)[length(cres)] <- "file.exposure"
      cres[length(cres) + 1] <- f2
      names(cres)[length(cres)] <- "file.outcome"
    }

    ret_cres <- as.data.frame(split(unname(cres), names(cres)))

    return(list(data.frame = ret_cres, regional = rp, zscore = zp))
  }, mc.cores = cores) %>%
    `c`(.)

  # This is not vectorised but could be!
  coloc.df <- tibble::tibble(file.exposure = character(),
                             file.outcome = character(),
                             nsnps = character(),
                             PP.H0.abf = character(),
                             PP.H1.abf = character(),
                             PP.H2.abf = character(),
                             PP.H3.abf = character(),
                             PP.H4.abf = character())
  for (i in 1:length(coloc_res)) {
    if (all(is.na(coloc_res[[i]]))) {
      next
    }
    coloc.df <- rbind(coloc.df, coloc_res[[i]]$data.frame)
  }

  # This vectorised version is iffy
  #coloc.df <- t(as.data.frame(sapply(coloc_res[which(!is.na(sapply(coloc_res, '[[', 1)))], '[[', 1)))
  #coloc.df <- as.data.frame(coloc.df)
  #rownames(coloc.df) <- NULL
  #coloc.df[, 3:8] <- sapply(coloc.df[, 3:8], as.numeric)

  regional.plots <- sapply(coloc_res, function(x) { x$regional })
  zscore.plots <- sapply(coloc_res, function(x) { x$zscore })

  return(list(res = coloc.df, regional = regional.plots, zscore = zscore.plots))
}

#' Sub-function for the colocalisation analyses
#'
#' @param dat1 SNPs, etc. from first dataset
#' @param dat2 SNPs, etc. from second dataset
#' @param min_snps Number of minimum SNPs to check for analysis to continue (Optional)
#' @param p1 p1 for coloc (Optional)
#' @param p2 p2 for coloc (Optional)
#' @param p12 p12 for coloc (Optional)
#' @param bfile Path to Plink bed/bim/fam files (Optional; required for SuSiE)
#' @param plink Path to Plink binary (Optional; required for SuSiE)
#' @param susie Run SuSiE? (Optional, boolean)
#' @param verbose Display verbose information (Optional, boolean)
#'
#' @return Results data.frame
#' @keywords Internal
#' @importFrom dplyr coalesce
#' @importFrom coloc coloc.abf
.coloc_sub <- function(dat1, dat2,
                       min_snps = 100,
                       p1 = 1e-4,
                       p2 = 1e-4,
                       p12 = 1e-5,
                       susie = FALSE,
                       bfile = NULL,
                       plink = NULL,
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
  #if (all(is.na(dat1$MAF)) || all(is.na(dat2$MAF)))
  #{
  #  .print_msg("Minor allele frequenices are required for coloc but were not found in these datasets.", verbose)
  #  return(NULL)
  #}
  dat1[is.na(dat1$MAF)] <- 0.5
  dat2[is.na(dat2$MAF)] <- 0.5

  if (length(dat1$snp) < min_snps || length(dat2$snp) < min_snps)
  {
    .print_msg(paste0("Datasets did not contain enough SNPs (< ", min_snps, ") for colocalisation. Skipping."), verbose)
    return(NULL)
  }

  if (susie)
  {
    if (is.null(plink) || is.null(bfile))
    {
      .print_msg("Coloc with SuSiE requires Plink and a reference panel; running only coloc.", verbose)
      susie <- FALSE
    } else
    {
      attempt <- try(cres <- .coloc_susie_sub(dat1, dat2,
                                              p1 = p1,
                                              p2 = p2,
                                              p12 = p12),
                     silent = T)

      if (inherits(attempt, "try-error") || is.na(cres)) {
        susie <- FALSE
      }
    }
  }

  # Not an "else if" so we can revert to just coloc if SuSiE cannot be run
  if (!susie)
  {
    attempt <- try(cres <- coloc::coloc.abf(dat1, dat2,
                                            p1 = p1,
                                            p2 = p2,
                                            p12 = p12),
                   silent = T)
  }

  if (inherits(attempt, "try-error")) {
    return(NA)
  }

  return(cres$summary)
}

#' Sub function to run SuSiE and coloc
#'
#' @param d1 Dataset 1
#' @param d2 Dataset 2
#' @param bfile Path to Plink bed/bim/fam files (Optional; required for SuSiE)
#' @param plink Path to Plink binary (Optional; required for SuSiE)
#' @param verbose Display verbose information (Optional, boolean)
#' @param ... Other arguments passed to coloc.susie and coloc.bf
#'
#' @return Results data.frame
#' @keywords Internal
#' @importFrom ieugwasr ld_matrix_local
#' @importFrom coloc runsusie coloc.susie
.coloc_susie_sub <- function(d1, d2,
                             bfile = NULL,
                             plink = NULL,
                             verbose = TRUE,
                             ...)
{
  # Datasets have already been harmonised so contain intersection of SNPs
  # Therefore, we can attempt to find the LD matrix straight away
  attempt <- tryCatch(
    expr = {
      ld <- ieugwasr::ld_matrix_local(d1$snp, bfile, plink, with_alleles = F)
    },
    error = function(e) {
      .print_msg(e, verbose)
    }
  )

  if (inherits(attempt, "error"))
  {
    .print_msg(".coloc_susie_sub: No LD matrix generated for coloc with SuSiE; running only coloc.", verbose)
    return(NA)
  }

  d1 <- d1[which(d1$snp %in% colnames(ld)), ]
  d1 <- d1[match(colnames(ld, d1$snp)), ]
  row.names(d1) <- 1:nrow(d1)
  d1 <- c(d1, ld = ld)

  d2 <- d2[which(d2$snp %in% colnames(ld)), ]
  d2 <- d2[match(colnames(ld, d2$snp)), ]
  row.names(d2) <- 1:nrow(d2)
  d2 <- c(d2, ld = ld)

  # Run SuSiE and store result
  s1 <- coloc::runsusie(d1)
  s2 <- coloc::runsusie(d2)

  return(coloc::coloc.susie(s1, s2, p1 = p1, p2 = p2, p12 = p12, ...))
}

#' Prepare gwasvcf files for coloc.
#' This method will extract SNPs from one file using one chrompos and then look
#' up those SNPs in the other file -- this is to ensure coloc can be conducted
#' upon two datasets of different genomic builds without the need of liftover.
#'
#' @param vcf1 VCF object or path to vcf file
#' @param vcf2 VCF object or path to vcf file
#' @param chrompos Character of the format chr:pos1-pos2
#'
#' @return list of coloc-ready data, or NA if failed
#' @importFrom tidyr drop_na
#' @importFrom gwasvcf query_chrompos_file query_chrompos_vcf query_gwas
#' @importFrom gwasglue gwasvcf_to_TwoSampleMR
#' @keywords Internal
.gwasvcf_to_coloc_rsid <- function(vcf1, vcf2, chrompos,
                                   type1 = NULL, type2 = NULL,
                                   build1 = "GRCh37", build2 = "GRCh37",
                                   verbose = TRUE)
{
  if (is.character(vcf1)) {
    r1 <- gwasvcf::query_chrompos_file(chrompos, vcf1, build = build1)
  }
  else if (class(vcf1) %in% c("CollapsedVCF", "ExpandedVCF")) {
    r1 <- gwasvcf::query_chrompos_vcf(chrompos, vcf1)
  }
  if (length(r1) < 1) {
    .pring_msg(paste0(".gwasvcf_to_coloc_rsid: Could not extract SNPs in region \"", chrompos, "\" for file, \"", vcf1, "\". Skipping."), verbose)
    return(NA)
  }
  r1 <- r1 %>% gwasglue::gwasvcf_to_TwoSampleMR() %>%
    tidyr::drop_na(pval.exposure, samplesize.exposure, beta.exposure, se.exposure, pval.exposure, SNP)

  r2 <- gwasvcf::query_gwas(vcf2, rsid = unique(r1$SNP), build = build2)
  if (length(r2) < 1) {
    .pring_msg(paste0(".gwasvcf_to_coloc_rsid: Could not extract SNPs in region \"", chrompos, "\" for file, \"", vcf2, "\". Skipping."), verbose)
    return(NA)
  }
  r2 <- r2 %>% gwasglue::gwasvcf_to_TwoSampleMR(type = "outcome") %>%
    tidyr::drop_na(pval.outcome, samplesize.outcome, beta.outcome, se.outcome, pval.outcome, SNP)

  # Get overlap
  # TODO Needs to be made better!
  #index <- as.character(tab1$REF) == as.character(tab2$ea) &
  #  as.character(tab1$ALT) == as.character(tab2$nea) &
  #  as.character(tab1$seqnames) == as.character(tab2$rsid) &
  #  tab1$start == tab2$position
  tab1 <- r1[r1$SNP %in% r2$SNP, ]
  tab2 <- r2[r2$SNP %in% tab1$SNP, ]
  tab1 <- tab1[tab1$SNP %in% tab2$SNP, ]

  if (nrow(tab1) < 1 || nrow(tab2) < 1) {
    .print_msg(paste0(".gwasvcf_to_coloc_rsid: No SNPs matched based on rsID between \"", vcf1, "\" and \"", vcf2, "\". Skipping coloc analysis for this pair."), verbose = verbose)
    return(NA)
  }

  # Need to convert AF to MAF
  tab1$maf.exposure <- ifelse(tab1$eaf.exposure > 0.5, 1 - tab1$eaf.exposure, tab1$eaf.exposure)
  tab1$maf.exposure <- ifelse(is.null(tab1$maf.exposure), 0.5, tab1$maf.exposure)
  tab2$maf.outcome <- ifelse(tab2$eaf.outcome > 0.5, 1 - tab2$eaf.outcome, tab2$eaf.outcome)
  tab2$maf.outcome <- ifelse(is.null(tab2$maf.outcome), 0.5, tab2$maf.outcome)

  # Attempt to determine from data types
  if (is.null(type1) && "ncase.exposure" %in% names(tab1) && all(!is.na(tab1$ncase.exposure))) {
    type1 <- "cc"
  } else {
    type1 <- "quant"
  }
  if (is.null(type2) && "ncase.outcome" %in% names(tab2) && all(!is.na(tab2$ncase.outcome))) {
    type2 <- "cc"
  } else {
    type2 <- "quant"
  }

  out1 <- tab1 %>%
    { list(pvalues = .$pval.exposure,
           N = .$samplesize.exposure,
           MAF = .$maf.exposure,
           beta = .$beta.exposure,
           varbeta = .$se.exposure^2,
           type = type1,
           snp = .$SNP,
           z = .$beta.exposure / .$se.exposure,
           #chr = .$chr.exposure,
           #pos = .$pos.exposure,
           id = vcf1)
    }

  if(type1 == "cc")
  {
    out1$s <- tab1$ncase.exposure / tab1$samplesize.exposure
  }

  out2 <- tab2 %>%
    { list(pvalues = .$pval.outcome,
           N = .$samplesize.outcome,
           MAF = .$maf.outcome,
           beta = .$beta.outcome,
           varbeta = .$se.outcome^2,
           type = type2,
           snp = .$SNP,
           z = .$beta.outcome / .$se.outcome,
           #chr = .$chr.outcome,
           #pos = .$pos.outcome,
           id = vcf2)
    }

  if(type2 == "cc")
  {
    out2$s <- tab2$ncase.outcome / tab2$samplesize.outcome
  }

  return(list(dataset1 = out1, dataset2 = out2))
}

#' Prepare gwasvcf files for PWCoCo
#'
#' Write files for PWCoCo where data are read from two VCF objects or files.
#'
#' @param vcf1 VCF object or path to VCF file
#' @param vcf2 VCF object or path to VCF file
#' @param chrompos Character of the format chr:pos1-pos2
#' @param type1 How to treat vcffile1 for coloc, either "quant" or "cc" (Optional)
#' @param type2 How to treat vcffile2 for coloc, either "quant" or "cc" (Optional)
#' @param outfile Path to output files, without file ending
#'
#' @return 0 if success, 1 if there was a problem
#' @keywords Internal
#' @importFrom gwasvcf vcflist_overlaps vcf_to_granges
#' @importFrom dplyr as_tibble select rename any_of
#' @importFrom VariantAnnotation header meta
.gwasvcf_to_pwcoco <- function(vcf1, vcf2, chrompos, type1=NULL, type2=NULL, outfile)
{
  overlap <- gwasvcf::vcflist_overlaps(list(vcf1, vcf2), chrompos)
  vcf1 <- overlap[[1]]
  vcf2 <- overlap[[2]]

  if (length(vcf1) == 0 || length(vcf2) == 0)
  {
    message(".gwasvcf_to_pwcoco: No overlaps for the given chrompos in ", ifelse(length(vcf1) == 0, "vcf1", "vcf2"), ".")
    return(1)
  }

  # vcf1
  tib1 <- vcf1 %>% gwasvcf::vcf_to_granges() %>% dplyr::as_tibble() %>%
    dplyr::select(dplyr::any_of(c("rsid", "ALT", "REF", "AF", "ES", "SE", "LP", "SS", "NC"))) %>%
    dplyr::rename(
      SNP = rsid,
      A1 = ALT,
      A2 = REF,
      freq = AF,
      b = ES,
      se = SE,
      p = LP,
      N = ss,
      N_case = NC,
      .cols = dplyr::any_of()
    )
  tib1$p <- 10^(-tib1$p)

  # Coloc type -- if study type is continuous then do not need the case column
  if (type1 == "quant" || VariantAnnotation::header(vcf1) %>% VariantAnnotation::meta() %>% {.[["SAMPLE"]][["StudyType"]]} == "Continuous")
  {
    tib1 <- tib1[c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")]
  }

  # vcf2
  tib2 <- vcf2 %>% gwasvcf::vcf_to_granges() %>% dplyr::as_tibble() %>%
    dplyr::select(dplyr::any_of(c("rsid", "ALT", "REF", "AF", "ES", "SE", "LP", "SS", "NC"))) %>%
    dplyr::rename(
      SNP = rsid,
      A1 = ALT,
      A2 = REF,
      freq = AF,
      b = ES,
      se = SE,
      p = LP,
      N = ss,
      N_case = NC,
      .cols = dplyr::any_of()
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

#' Prepare ieugwasr data for PWCoCo
#'
#' Write files for PWCoCo where data are read from the OpenGWAS DB.
#'
#' @param id1 ID for trait 1
#' @param id2 ID for trait 2
#' @param chrompos Character of the format chr:pos1-pos2
#' @param type1 How to treat vcffile1 for coloc, either "quant" or "cc" (Optional)
#' @param type2 How to treat vcffile2 for coloc, either "quant" or "cc" (Optional)
#' @param outfile Path to output files, without file ending
#'
#' @return 0 if success, 1 if there was a problem
#' @keywords Internal
#' @importFrom ieugwasr associations gwasinfo
#' @importFrom dplyr select rename
.ieugwasr_to_pwcoco <- function(id1, id2, chrompos, type1=NULL, type2=NULL, outfile)
{
  tib1 <- ieugwasr::associations(id=id1, variants=chrompos) %>% subset(., !duplicated(rsid))
  tib2 <- ieugwasr::associations(id=id2, variants=chrompos) %>% subset(., !duplicated(rsid))

  if (length(tib1) < 1 || length(tib2) < 1)
  {
    message(".ieugwasr_to_pwcoco: Data could not be read using the ieugwasr package for id1 = ", id1, " and id2 = ", id2, ".")
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

#' Sub-function to run PWCoCo
#'
#' @param bfile Path to Plink bed/bim/fam files
#' @param chrpos Character of the format chr:pos1-pos2
#' @param pwcoco Path to PWCoCo executible
#' @param maf MAF cut-off (Optional)
#' @param p1 p1 for coloc (Optional)
#' @param p2 p2 for coloc (Optional)
#' @param p12 p12 for coloc (Optional)
#' @param workdir Path to save temporary files (Optional)
#' @param verbose Display verbose information (Optional, boolean)
#'
#' @return Results data.frame
#' @keywords Internal
#' @importFrom glue glue
#' @importFrom data.table fread
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
