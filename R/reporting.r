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
    warning("Could not find column names in MR results data.frame: ", paste(names, collapse = ", "))
  }

  if (!all(coloc_res.by %in% names(coloc_res.df))) {
    warning("Could not find column names in coloc results data.frame: ",)
  }

  merged <- base::merge(mr_res, coloc_res.df, by.x = mr_res.by, by.y = coloc_res.by, all.x = TRUE)
}
