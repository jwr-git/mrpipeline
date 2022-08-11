#' Harmonise exposure and outcomes.
#'
#' Wrapper function for \link[TwoSampleMR]{harmonise_data} function in the
#' `TwoSampleMR` package.
#' @seealso [TwoSampleMR::harmonise_data()]
#'
#' @param exposure Data.frame of exposure dataset(s)
#' @param outcome Data.frame of outcome dataset(s)
#' @param action How to harmonise alleles; see \link[TwoSampleMR]{harmonise_data}.
#' @param cores Number of cores for multi-threaded tasks (Optional)
#'              NB: Unavailable on Windows machines
#' @param verbose Display verbose information (Optional, boolean)
#'
#' @return Harmonised data.frame
#' @importFrom TwoSampleMR harmonise_data
#' @importFrom tidyr crossing
#' @importFrom parallel mclapply
#' @export
harmonise <- function(exposure,
                      outcome,
                      action = 1,
                      cores = 1,
                      verbose = TRUE)
{
  if (!action %in% 1:3) {
    stop("\"action\" not recognised -- must take a value between 1 and 3!")
  }

  if (action == 1) {
    warning("Using action 1 is very lenient and not recommended! See the vignettes for more details.")
  }

  pairs <- tidyr::crossing(exposure$id.exposure, outcome$id.outcome)
  names(pairs) <- c("id.exposure", "id.outcome")
  dat <- parallel::mclapply(1:nrow(pairs), function(i)
  {
    exp <- exposure[exposure$id.exposure == pairs[i, "id.exposure"][[1]], ]
    out <- outcome[outcome$id.outcome == pairs[i, "id.outcome"][[1]], ]

    harmonised <- TwoSampleMR::harmonise_data(exp, out, action = action)

    if (!nrow(harmonised) | all(is.na(harmonised))) {
      return(NULL)
    }

    return(harmonised)
  }, mc.cores = cores) %>%
    dplyr::bind_rows()
}

#' Check if SNP is palindromic
#'
#' @param A1 Allele 1
#' @param A2 Allele 2
#'
#' @return True/False if palindromic
is_palindromic <- function(A1, A2)
{
  (A1 == "A" & A2 == "T") |
    (A1 == "T" & A2 == "A") |
    (A1 == "C" & A2 == "G") |
    (A1 == "G" & A2 == "C")
}

#' Check if SNP frequencies are ambiguous
#'
#' @param freq Frequency
#' @param tol Tolerance around 0.5 (Optional)
#'
#' @return True/False if ambiguous
is_ambiguous <- function(freq, tol = 0.08)
{
  freq > 0.5 - tol & freq < 0.5 + tol
}

#' Performs pairwise harmonisation and analyses -- helpful when analysing many
#' exposure-outcome pairs as performing the standard "linear" approach will
#' be very slow.
#'
#' In this function, an exposure-outcome pair are harmonised, analyses are ran
#' on those data and those results are saved to a file. Analyses ran can be MR
#' or colocalisation, as desired.
#'
#' @param exposure Data.frame of exposure dataset(s)
#' @param outcome Data.frame of outcome dataset(s)
#' @param res_path Path to save result files
#' @param overwrite If the results file exists, should it be overwritten?
#' @param ... Other arguments for the following functions:
#'            \link[mrpipeline]{harmonise}
#'            \link[mrpipeline]{do_mr}
#'            \link[mrpipeline]{do_coloc}
#' @param do_coloc True/False run colocalisation analyses
#' @param cores Number of cores for multi-threaded tasks (Optional)
#'              NB: Unavailable on Windows machines
#' @param verbose Display verbose information (Optional, boolean)
#'
#' @importFrom tidyr crossing
#' @importFrom parallel mclapply
#'
#' @export
pairwise_analysis <- function(exposure,
                              outcome,
                              res_path,
                              overwrite = F,
                              action = 1,
                              f_cutoff = 10,
                              all_wr = TRUE,
                              do_coloc = FALSE,
                              cores = 1,
                              verbose = TRUE,
                              ...)
{
  if (!dir.exists(res_path)) {
    stop("Pathway does not exist: \"", res_path, "\".")
  }
  pairs <- tidyr::crossing(exposure$id.exposure, outcome$id.outcome)
  names(pairs) <- c("id.exposure", "id.outcome")
  parallel::mclapply(1:nrow(pairs), function(i)
  {
    exp <- exposure[exposure$id.exposure == pairs[i, "id.exposure"][[1]], ]
    out <- outcome[outcome$id.outcome == pairs[i, "id.outcome"][[1]], ]

    exp <- check_snps(exp, analyses = "mr")
    out <- check_snps(out, analyses = "mr")

    if (nrow(exp) == 0 || nrow(out) == 0 || length(exp) == 0 || length(out) == 0)
    {
      return(NULL)
    }

    exp_name <- basename(tools::file_path_sans_ext(pairs[i, "id.exposure"][[1]]))
    out_name <- basename(tools::file_path_sans_ext(pairs[i, "id.outcome"][[1]]))

    mr_file_name <- paste0(res_path,
                           "/",
                           exp_name,
                           "_",
                           out_name,
                           "_mr_results.txt")

    if (!file.exists(mr_file_name) || overwrite)
    {
      dat <- harmonise(exp, out, action = action, cores = cores, verbose = verbose)

      if (nrow(dat) == 0 || length(dat) == 0) {
        return(NULL)
      }

      mr_res <- do_mr(dat, f_cutoff = f_cutoff, all_wr = all_wr, verbose = verbose)
      if (is.na(mr_res)) {
        return(NULL)
      }
      write.table(mr_res, file = mr_file_name, ...)
    }

    coloc_file_name <- paste0(res_path,
                              "/",
                              exp_name,
                              "_",
                              out_name,
                              "_coloc_results.txt")

    if (do_coloc && (!file.exists(coloc_file_name) || overwrite))
    {
      coloc_res <- do_coloc(dat, ..., cores = cores, verbose = verbose) # TODO remove `...`
      write.table(coloc_res, file = coloc_file_name, ...)
    }

  }, mc.cores = cores)
}
