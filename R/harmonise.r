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
    exp <- exposure[exposure$id.exposure == pairs[i, "id.exposure"], ]
    out <- outcome[outcome$id.outcome == pairs[i, "id.outcome"], ]

    TwoSampleMR::harmonise_data(exp, out, action = action)
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
