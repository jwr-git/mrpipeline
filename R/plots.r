f_plot <- function(dat, report,
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
          ggtitle(paste0("Plot of F-statistic against P value for ",trait)) +
          stat_smooth(method = "lm", formula = y ~ poly(x, 2), alpha = 0.5, colour = "grey")

    report$add_plot(list(id1 = id,
                         id2 = NA,
                         name1 = unique(plot_dat[[trait_col]]),
                         name2 = "",
                         type = "f"),
                    f)
  }
}

interactive_z_plot <- function(dat)
{
  f <- ggplot2::ggplot(data = dat, aes(x = beta.exposure / se.exposure,
                                            y = -log10(pval.exposure))) +
    geom_point(size = 2, alpha = 0.8, aes(text = sprintf("SNP: %s", SNP))) +
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
    ggtitle(paste0("Plot of Z scores against P value")) +
    stat_smooth(method = "lm", formula = y ~ poly(x, 2), alpha = 0.5, colour = "grey")

  return(f)
}

volcano_plot <- function(res, report)
{
  if (!nrow(res)) {
    return()
  }

  id <- unique(res[["id.exposure"]])
  name <- unique(res[["exposure"]])
  res <- res[res$method %in% c("Wald ratio", "Inverse variance weighted"), ]

  if (is.na(id) || !nrow(res)) {
    return()
  }

  p <- ggplot2::ggplot(data = res, aes(x = b, y = -log10(pval), label = outcome)) +
    ggplot2::geom_point(aes(shape = method), size=2) +
    #ggrepel::geom_label_repel(box.padding=0.35,
    #                 point.padding=0.5,
    #                 size=6,
    #                 segment.colour="grey50") +
    ggplot2::theme_bw() +
    ggplot2::xlab("Odds ratio") +
    ggplot2::ylab("-log10(P value)")

  report$add_plot(list(id1 = id,
                       id2 = NA,
                       name1 = name,
                       name2 = "",
                       type = "volcano"),
                  p)
}

phewas_plot <- function(res, report)
{
  if (!nrow(res)) {
    return()
  }

  id <- unique(res[["id.exposure"]])
  name <- unique(res[["exposure"]])
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

  report$add_plot(list(id1 = id,
                       id2 = NA,
                       name1 = name,
                       name2 = "",
                       type = "phewas"),
                  p)
}

regional_plot <- function(dat, pairs, i, report)
{
  if (require("gassocplot"))
  {
    tryCatch(
      expr = {
        dat <- gwasglue::coloc_to_gassocplot(dat)
        p <- gassocplot::stack_assoc_plot(dat$markers, dat$z, dat$corr, traits = dat$traits)

        # Annotate titles
        # TODO Better way of referencing these than hard-coded numbers?
        tryCatch({
          p$grobs[[1]]$grobs[[22]]$label <- pairs[i, "trait.x"] # Top
          p$grobs[[1]]$grobs[[65]]$label <- pairs[i, "trait.y"] # Bottom
        })

        # a <- grid.grabExpr(gassocplot::stack_assoc_plot(dat$markers, dat$z, dat$corr, traits = dat$traits)) %>%
        # editGrob()

        report$add_plot(list(id1 = dat$traits[1],
                             id2 = dat$traits[2],
                             name1 = pairs[i, "trait.x"],
                             name2 = pairs[i, "trait.y"],
                             type = "regional"),
                        p)
      },
      error = function(e) {
        report$add_plot(list(id1 = dat$traits[1],
                             id2 = dat$traits[2],
                             name1 = pairs[i, "trait.x"],
                             name2 = pairs[i, "trait.y"],
                             type = "regional"),
                        NA)
      }
    )
  }
}

mr_scatter_plot <- function(res, keys, dat, report)
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
    report$add_plot(list(id1 = keys[[1]],
                         id2 = keys[[2]],
                         name1 = "",
                         name2 = "",
                         type = "scatter"), NA)
  }

  # Plot!
  p <- TwoSampleMR::mr_scatter_plot(res, dat)
  report$add_plot(list(id1 = keys[[1]],
                       id2 = keys[[2]],
                       name1 = "",
                       name2 = "",
                       type = "scatter"),
                  p)
}

interactive_scatter_plot <- function(dat, id.exposure, id.outcome)
{
  # Each unique pair of traits need to be colocalised
  subdat <- dat[dat$id.exposure == id.exposure & dat$id.outcome == id.outcome, ]

  if (!length(subdat)) {
    return(NULL)
  }

  f <- ggplot2::ggplot(data = subdat, aes(x = beta.exposure,
                                          y = beta.outcome)) +
    geom_point(size = 2, alpha = 0.8, aes(text = sprintf("SNP: %s", SNP))) +
    geom_errorbar(alpha = 0.3,
                  width = 0.001,
                  aes(ymin = beta.outcome - 1.96 * se.outcome,
                      ymax = beta.outcome + 1.96 * se.outcome)) +
    geom_errorbarh(alpha = 0.3,
                   height = 0.001,
                   aes(xmin = beta.exposure - 1.96 * se.exposure,
                       xmax = beta.exposure + 1.96 * se.exposure)) +
    theme_bw() +
    geom_vline(xintercept = 0.0, col = "black", linetype = "dashed") +
    geom_hline(yintercept = 0.0, col = "darkgrey", linetype = "dashed") +
    xlab("Exposure Beta") +
    ylab("Outcome Beta") +
    ggtitle(paste0("Plot of SNP effects on Exposure vs Outcome"))

  return(f)
}

forest_plot <- function(res, dat, report)
{

}
