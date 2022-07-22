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
                              ...,
                              do_coloc = FALSE,
                              cores = 1,
                              verbose = TRUE)
{
  if (!dir.exists(res_path)) {
    stop("Pathway does not exist: \"", res_path, "\".")
  }
  pairs <- tidyr::crossing(exposure$id.exposure, outcomes$id.outcome)
  names(pairs) <- c("id.exposure", "id.outcome")
  parallel::mclapply(1:nrow(pairs), function(i)
  {
    exp <- exposure[exposure$id.exposure == pairs[i, "id.exposure"][[1]], ]
    out <- outcome[outcome$id.outcome == pairs[i, "id.outcome"][[1]], ]

    dat <- harmonise(exp, out, ..., cores, verbose)

    mr_res <- do_mr(dat, ..., verbose = verbose)
    write.table(mr_res,
                file = paste0(res_path,
                              "/",
                              pairs[i, "id.exposure"][[1]],
                              "_",
                              pairs[i, "id.outcome"][[1]],
                              "_mr_results.txt"),
                ...)

    if (do_coloc)
    {
      coloc_res <- do_coloc(dat, ..., cores = cores, verbose = verbose)
      write.table(coloc_res,
                  file = paste0(res_path,
                                "/",
                                pairs[i, "id.exposure"][[1]],
                                "_",
                                pairs[i, "id.outcome"][[1]],
                                "_coloc_results.txt"),
                  ...)
    }

  }, mc.cores = cores)
}
