#' Combine MR and coloc results into one data.frame
#'
#' @param mr_res A data.frame of MR results from `do_mr()`
#' @param coloc_res A data.frame of coloc results from `do_coloc()`
#' @param mr_res.by MR columns to use for merging
#' @param coloc_res.by Coloc columns to use for merging
#'
#' @return Data.frame of merged results
#' @export
combine_results <- function(mr_res,
                            coloc_res,
                            mr_res.by = c("id.exposure", "id.outcome"),
                            coloc_res.by = c("file.exposure", "file.outcome")
)
{
  if (!length(mr_res) || !length(coloc_res)) {
    warning("Required at least results from MR and colocalisation to combine.")
    return(NA)
  }

  if ("plots" %in% names(coloc_res)) {
    coloc_res.df <- coloc_res$res
    coloc_res.plots <- coloc_res$plots
  }

  if (!all(mr_res.by %in% names(mr_res))) {
    warning("Could not find column names in MR results data.frame: ", paste(mr_res.by, collapse = ", "))
  }

  if (!all(coloc_res.by %in% names(coloc_res.df))) {
    warning("Could not find column names in coloc results data.frame: ", paste(coloc_res.by, collapse = ", "))
  }

  merged <- base::merge(mr_res, coloc_res.df, by.x = mr_res.by, by.y = coloc_res.by, all.x = TRUE)
  return(merged)
}
