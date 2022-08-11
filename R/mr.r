#' Calculate F-statistic
#'
#' Calculates portion of variance explained and F-statistic.
#' If the data is lacking key information, i.e. allele frequencies, sample size
#' or consists of only one SNP, then the approximate F-statistic will be used
#' instead: \eqn{F = b ** 2 / SE ** 2}.
#'
#' @seealso [do_mr()]
#'
#' @param dat Data.frame from do_mr()
#' @param f_cutoff F-statistic cutoff (Optional)
#' @param force_approx Force to use the approximate F-statistic instead (Optional, boolean)
#' @param verbose Display verbose information (Optional, boolean)
#'
#' @return Modified `dat` data.frame (if f_cutoff > 0 supplied)
#' @export
#' @importFrom plyr ddply
calc_f_stat <- function(dat, f_cutoff = 10, force_approx = FALSE, verbose = TRUE)
{
  dat <- plyr::ddply(dat, c("id.exposure"), function(x)
  {
    full.f.stat <- !force_approx

    if (full.f.stat)
    {
      if (any(is.na(x$eaf.exposure))) {
        warning("Using approximate F-statistic as some allele frequencies are missing.")
        full.f.stat <- F
      }
      else if (any(is.na(x$samplesize.exposure))) {
        warning("Using approximate F-statistic as some sample sizes are missing.")
        full.f.stat <- F
      }
      else if (length(unique(x$SNP)) == 1) {
        full.f.stat <- F
      }
    }

    # PVE / F-statistic
    if (full.f.stat) {
      x$maf.exposure <- ifelse(x$eaf.exposure < 0.5, x$eaf.exposure, 1 - x$eaf.exposure)
      x$pve.exposure <- mapply(.calc_pve,
                               x$beta.exposure,
                               x$maf.exposure,
                               x$se.exposure,
                               x$samplesize.exposure)
      # From https://doi.org/10.1093/ije/dyr036
      x$f.stat.exposure <- ifelse(x$pve.exposure == -1, 0,
                                  ((x$samplesize.exposure - length(unique(x$SNP)) - 1) / length(unique(x$SNP))) * (x$pve.exposure / (1 - x$pve.exposure)))
    } else {
      x$maf.exposure <- NA
      x$pve.exposure <- NA
      x$f.stat.exposure <- x$beta.exposure ** 2 / x$se.exposure ** 2
    }
    x
  })

  rem.f.stat <- nrow(dat[dat$f.stat.exposure < f_cutoff, ])
  .print_msg(paste0("Amount of SNPs which did not meet threshold: ", rem.f.stat), verbose)
  #dat <- dat[dat$f.stat.exposure >= f_cutoff & !is.na(dat$f.stat.exposure), ]

  return(dat)
}

#' Calculate PVE
#'
#' Calculates proportion of variance explained.
#' From \href{https://doi.org/10.1371/journal.pone.0120758}{S1 Text}
#'
#' @param b Vector or number, beta
#' @param maf Vector or number, minor allele frequency
#' @param se Vector or number, standard error of beta
#' @param n Vector or number, sample size
#'
#' @return Vector or number, proportion of variance explained
#' @keywords Internal
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

#' WR Taylor Approx of SE
#'
#' Calculates the second term Taylor approximation for standard error of
#' the Wald ratio method.
#' From \href{https://doi.org/10.1101/2021.03.01.433439}{supplementary}
#'
#' @param object Harmonised data.frame
#'
#' @return Results data.frame
#' @keywords Internal
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

#' IVW weighted delta
#'
#' Calculates the inverse variance weighted delta method from the MendelianRandomization package
#'
#' @param object Harmonised data.frame
#'
#' @return Results data.frame
#' @keywords Internal
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

#' Run Mendelian randomisation analyses
#'
#' Runs Mendelian randomisation and related analyses:
#' \enumerate{
#' \item Wald ratio, see \link[mrpipeline]{.wr_taylor_approx()}
#' \item Inverse variance weighted, see \link[mrpipeline]{.ivw.delta()}
#' \item Steiger filtering, see TwoSampleMR::directionality_test()
#' }
#'
#' @seealso [.wr_taylor_approx()], [.ivw_delta()], [TwoSampleMR::directionality_test()]
#'
#' @param dat A data.frame of harmonised data
#' @param f_cutoff Define an F-statistic cutoff (Optional)
#' @param all_wr Should the Wald ratio be calculated for all SNPs, even if IVW can be used? (Optional)
#' @param verbose Display verbose information (Optional, boolean)
#'
#' @return A data.frame of MR results
#' @export
#' @importFrom plyr ddply
#' @importFrom TwoSampleMR generate_odds_ratios directionality_test
do_mr <- function(dat, f_cutoff = 10, all_wr = TRUE, verbose = TRUE)
{
  if (!is.null(f_cutoff)) {
    if ("f.stat.exposure" %in% names(dat)) {
      dat <- dat[dat$f.stat.exposure >= f_cutoff & !is.na(dat$f.stat.exposure), ]
    } else {
      .print_msg("F-statistic cut-off given but could not find column \"f.stat.exposure\". Calculating these now.", verbose = verbose)
      dat <- calc_f_stat(dat, f_cutoff = f_cutoff, verbose = verbose)
    }
  }

  if (!nrow(dat) || !length(dat)) {
    return(NA)
  }

  res <- plyr::ddply(dat, c("id.exposure", "id.outcome"), function(x1)
  {
    x <- subset(x1, mr_keep.exposure)

    nsnps <- nrow(x)

    # WR on all single SNPs
    if (nsnps == 1 || (nsnps > 1 && all_wr))
    {
      res <- lapply(1:nsnps, function(i)
      {
        with(x, .wr_taylor_approx(x[i, ]))
      })
      res <- do.call(rbind, res)
      wr_res <- TRUE
    }
    else {
      wr_res <- FALSE
    }

    # IVW if applicable
    if (nsnps > 1)
    {
      tryCatch(
        expr = {
          res_ <- .ivw_delta(x)
          if (wr_res) {
            res <- rbind(res, res_)
          } else {
            res <- res_
          }
        },
        error = function(e) {
          message("Error encounted in IVW, no result for this will be given: ")
          message(e)
        }
      )
    }

    if (nrow(res) > 0) {
      res <- subset(res, !(is.na(b) & is.na(se) & is.na(pval))) %>%
        TwoSampleMR::generate_odds_ratios()
    }
  })

  # MR-Egger interceptt
  egger_intercept_res <- TwoSampleMR::mr_pleiotropy_test(dat)
  if (nrow(egger_intercept_res) > 1) {
    res <- base::merge(res, egger_intercept_res, by = c("exposure", "outcome", "id.exposure", "id.outcome"), suffixes = c("", ".egger"))
  } else {
    res$egger_intercept <- NA
    res$se.egger <- NA
    res$pval.egger <- NA
  }

  # Steiger
  if (any(is.na(dat$samplesize.exposure)) || any(is.na(dat$samplesize.outcome))) {
    warning("Samplesizes are required for Steiger filtering.")
    res$snp_r2.exposure <- NA
    res$snp_r2.outcome <- NA
    res$correct_causal_direction <- NA
    res$steiger_pval <- NA
    res$steigerflag <- NA
  } else {
    steigerres <- TwoSampleMR::directionality_test(dat)

    if (is.null(steigerres)) {
      res$snp_r2.exposure <- NA
      res$snp_r2.outcome <- NA
      res$correct_causal_direction <- NA
      res$steiger_pval <- NA
      res$steigerflag <- NA
    }
    else {
      steigerres$steigerflag <- ifelse(steigerres$correct_causal_direction == T & steigerres$steiger_pval < 0.05,
                                       "True",
                                       ifelse(steigerres$correct_causal_direction == F & steigerres$steiger_pval < 0.05,
                                              "False",
                                              "Unknown"))

      res <- base::merge(res, steigerres, by = c("exposure", "outcome", "id.exposure", "id.outcome"), all.x = TRUE)
    }
  }

  return(res)
}
