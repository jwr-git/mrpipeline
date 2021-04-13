pipeline_f_plot <- function(dat, report,
                                  beta_col = "beta.exposure",
                                  se_col = "se.exposure",
                                  pval_col = "pval.exposure",
                                  trait_col = "exposure",
                                  id_col = "id.exposure")
{
  # TODO Could probably do this in a more "R-way" -- apply?
  for(trait in unique(get(trait_col, dat))) {
    plot_dat <- subset(dat, get(trait_col) == trait)
    id <- unique(plot_dat[[id_col]])

    f <- ggplot2::ggplot(data = plot_dat, aes(x = !!sym(beta_col) / !!sym(se_col),
                                                    y = -log10(!!sym(pval_col)))) +
          geom_point(size = 2, alpha = 0.8) +
          theme_bw() +
          #theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          #      legend.position = c(0.9, 0.9), legend.background = element_rect(size=0.5, linetype="solid",
          #                                                                      colour ="darkgrey"),
          #      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          #      legend.text=element_text(size=20), legend.title = element_text(size=20)) +
          geom_vline(xintercept = 0.0, col = "black", linetype = "dotted") +
          geom_hline(yintercept = -log10(5e-8), col = "darkgrey", linetype = "twodash", size = 1) +
          xlab("Beta / SE") +
          ylab("-log10(P value)") +
          ggtitle(paste0("Plot of beta^2/SE^2 against P value for ",trait)) +
          stat_smooth(method = "lm", formula = y ~ poly(x, 2), alpha = 0.5, colour = "grey")

    report$add_plot(list(id1 = id, id2 = NA, type = "f"), f)
  }
}

volcano_plot <- function(res, report)
{
  if (!nrow(res)) {
    return()
  }

  id <- unique(res[["id.exposure"]])
  res <- res[res$method %in% c("Wald ratio", "Inverse variance weighted"), ]

  if (is.na(id) || !nrow(res)) {
    return()
  }

  p <- ggplot2::ggplot(data = res, aes(x = b, y = -log10(pval), label = outcome)) +
    ggplot2::geom_point(aes(shape = method), size=2) +
    ggrepel::geom_label_repel(box.padding=0.35,
                     point.padding=0.5,
                     size=6,
                     segment.colour="grey50") +
    ggplot2::theme_bw() +
    ggplot2::xlab("Odds ratio") +
    ggplot2::ylab("-log10(P value)")

  report$add_plot(list(id1 = id, id2 = NA, type = "volcano"), p)
}

phewas_plot <- function(res, report)
{
  if (!nrow(res)) {
    return()
  }

  id <- unique(res[["id.exposure"]])
  res <- res[res$method %in% c("Wald ratio", "Inverse variance weighted"), ]

  if (is.na(id) || !nrow(res)) {
    return()
  }

  p <- res %>%
    dplyr::mutate(or = scale(or)) %>%
    ggplot2::ggplot(data = res, mapping = aes(x = outcome, y = -log10(pval), size = or)) +
    ggplot2::geom_point() +
    ggplot2::scale_size_binned() +
    ggplot2::theme_bw() +
    ggplot2::xlab("Outcome") +
    ggplot2::ylab("-log10(P value)") +
    ggplot2::labs(size = "Standardised OR")

  report$add_plot(list(id1 = id, id2 = NA, type = "phewas"), p)
}

regional_plot <- function(dat, report)
{
  if (require("gassocplot"))
  {
    tryCatch(
      expr = {
        dat <- gwasglue::coloc_to_gassocplot(dat)
        p <- gassocplot::stack_assoc_plot(dat$markers, dat$z, dat$corr, traits = dat$traits)

        report$add_plot(list(id1 = dat$traits[1], id2 = dat$traits[2], type = "regional"), p)
      },
      error = function(e) {
        report$add_plot(list(id1 = dat$traits[1], id2 = dat$traits[2], type = "regional"), NA)
      }
    )
  }
}

scatter_plot <- function(res, keys, dat, report)
{
  # Mapping went wrong or no results
  if (!ncol(keys) || ncol(keys) > 2 || !nrow(res))
  {
    return()
  }

  # Subset using the mapping IDs
  dat <- dat[dat$id.exposure == keys[[1]] & dat$id.outcome == keys[[2]],]

  # Not enough SNPs to plot?
  if (nrow(dat) < 2)
  {
    report$add_plot(list(id1 = keys[[1]], id2 = keys[[2]], type = "scatter"), NA)
  }

  # Plot!
  p <- TwoSampleMR::mr_scatter_plot(res, dat)
  report$add_plot(list(id1 = keys[[1]], id2 = keys[[2]], type = "scatter"), p)
}
