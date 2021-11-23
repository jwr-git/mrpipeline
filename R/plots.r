#' Plots F-statistic against P value
#' UNIMPLEMENTED
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

#' Creates Z-score plot for SNPs
#'
#' @param dat A data.frame of harmonised data
#'
#' @return Plot
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

#' Plots Volcano plot
#' UNIMPLEMENTED
volcano_plot <- function(res, report)
{
  if (!nrow(res) || !length(res)) {
    return()
  }

  id.exposure <- unique(res[["id.exposure"]])
  id.outcome <- unique(res[["id.outcome"]])
  exposure <- unique(res[["exposure"]])
  outcome <- unique(res[["outcome"]])

  if (is.na(id.exposure) || !nrow(res)) {
    return()
  }

  p <- ggplot2::ggplot(data = res, aes(x = b, y = -log10(pval), label = outcome)) +
    ggplot2::geom_point(aes(text = sprintf("SNP: %s", snp)), size=2) +
    #ggrepel::geom_label_repel(box.padding=0.35,
    #                 point.padding=0.5,
    #                 size=6,
    #                 segment.colour="grey50") +
    ggplot2::theme_bw() +
    ggplot2::xlab("ln(Odds ratio)") +
    ggplot2::ylab("-log10(P value)")

  report$add_plot(list(id1 = id.exposure,
                       id2 = id.outcome,
                       name1 = exposure,
                       name2 = outcome,
                       type = "volcano"),
                  p)
}

#' Plots and returns an interactive volcano plot
#'
#' @param dat A data.frame of harmonised data
#'
#' @return Plot
interactive_volcano_plot <- function(dat)
{
  if (!nrow(dat) || !length(dat)) {
    return()
  }

  r <- dat %>%
    TwoSampleMR::mr_singlesnp() %>%
    TwoSampleMR::generate_odds_ratios() %>%
    dplyr::filter(!startsWith(getElement(., "SNP"), "All"))

  p <- r %>%
    plotly::highlight_key(~outcome) %>%
    plotly::plot_ly(
      x = ~b, y = ~-log10(p),
      type = 'scatter',
      mode = "markers+text",
      textposition = "top",
      hovertemplate = paste(
        "<b>", r$SNP, " : ", r$outcome, "</b><br><br>",
        "MR Beta = %{x:.2f}<br>",
        "OR (95% CI) = ", signif(r$or, 3), "(", signif(r$or_lci95, 3), ",", signif(r$or_uci95, 3), ")"
        )
    ) %>%
    plotly::highlight(on = "plotly_click", selectize = TRUE, dynamic = TRUE)

  return(p)
}

#' Plots a PheWAS-like plot of all exposures against one outcome
#'
#' @param res A data.frame of results from MR
#' @param report QCReport class of results, etc. for reporting
#'
#' @return NULL
phewas_plot <- function(res, report)
{
  if (!nrow(res)) {
    return()
  }

  id <- unique(res[["id.outcome"]])
  name <- unique(res[["outcome"]])
  res <- res[res$method %in% c("Wald ratio", "Inverse variance weighted"), ]

  if (is.na(id) || nrow(res) < 2) {
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

  report$add_plot(list(id1 = NA,
                       id2 = id,
                       name1 = name,
                       name2 = "",
                       type = "phewas"),
                  p)
}

#' Plots a regional plot of the area being tested for colocalisation
#'
#' @param dat A data.frame of harmonised data
#' @param exposure Character, name of exposure
#' @param outcome Character, name of outcome
#' @param bfile Path to Plink bed/bim/fam files
#' @param plink Path to Plink binary
#' @param verbose Print messages or not
#'
#' @return NULL
regional_plot <- function(dat, exposure, outcome, bfile = NULL, plink = NULL, verbose = TRUE)
{
  p <- NA
  if (require("gassocplot"))
  {
    tryCatch(
      expr = {
        dat <- gwasglue::coloc_to_gassocplot(dat, bfile = bfile, plink_bin = plink)
        p <- gassocplot::stack_assoc_plot(dat$markers, dat$z, dat$corr, traits = dat$traits)

        # Annotate titles
        # TODO Better way of referencing these than hard-coded numbers?
        tryCatch({
          p$grobs[[1]]$grobs[[22]]$label <- exposure # Top
          p$grobs[[1]]$grobs[[65]]$label <- outcome # Bottom
        })

        # a <- grid.grabExpr(gassocplot::stack_assoc_plot(dat$markers, dat$z, dat$corr, traits = dat$traits)) %>%
        # editGrob()
      },
      error = function(e) {
        .print_msg(paste0("Could not generate regional plot for \"", exposure, "\" and \"", outcome, "\".)"), verbose = verbose)
      }
    )
  }

  return(p)
}

#' Plots a scatter plot of MR results
#' UNIMPLEMENTED
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

#' Plots an interactive scatter plot of SNP-exposure and SNP-outcome effects
#'
#' @param dat A data.frame of harmonised data
#' @param id.exposure Exposure ID
#' @param id.outcome Outcome ID
#'
#' @return Plot
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

#' Forest plot of MR results
#' UNIMPLEMENTED
forest_plot <- function(res, dat, report)
{
  if (!nrow(res) || !length(res)) {
    return()
  }

  #res <- res[res$method %in% c("Wald ratio", "Inverse variance weighted"), ]

  #if (!nrow(res)) {
  #  return()
  #}

  # Reverse order of res so that WR are first followed by IVW
  res <- res[seq(dim(res)[1], 1), ]

  id <- unique(res[["id.exposure"]])
  name <- unique(res[["exposure"]])

  fontsize <- 0.7
  cpositions <- c(0.02, 0.22, 0.4)
  xlab <- ""
  ylab <- ""
  nodigits <- 2 # Number of digits to show
  main <- "Forest Plot of MR Results"

  # Build clean df
  toShow <- res
  toShow <- toShow[c("exposure", "outcome", "method", "snp", "b", "se", "pval", "lo_ci", "up_ci")]
  toShow$estimate <- round(exp(toShow$b), nodigits)
  toShow$conf.low <- round(exp(toShow$lo_ci), nodigits)
  toShow$conf.high <- round(exp(toShow$up_ci), nodigits)
  toShow$ci <- paste0("(",toShow$conf.low," - ",toShow$conf.high,")")
  toShow$p.value <- signif(toShow$pval, nodigits+1)
  toShow$stars <- paste0(round(toShow$pval, 3), " ",
                         ifelse(toShow$pval < 0.05, "*",""),
                         ifelse(toShow$pval < 0.01, "*",""),
                         ifelse(toShow$pval < 0.001,"*",""))
  toShow$stars[which(toShow$pval < 0.001)] = "<0.001 ***"
  toShow$ci[is.na(toShow$b)] = ""
  toShow$estimate[is.na(toShow$b)] = 0

  terms <- ifelse(toShow$method == "Wald ratio",
                  paste0(toShow$exposure, " | ", toShow$outcome, " | ", toShow$method, " | ", toShow$snp),
                  paste0(toShow$exposure, " | ", toShow$outcome, " | ", toShow$method))
  coeffs <- res$b
  lower <- res$lo_ci
  upper <- res$up_ci

  # Graphing stuff
  rangeb <- range(lower, upper, na.rm = TRUE)
  breaks <- axisTicks(rangeb/2, log = TRUE, nint = 7)
  rangeplot <- rangeb
  rangeplot[1] <- rangeplot[1] - diff(rangeb)
  rangeplot[2] <- rangeplot[2] + .15 * diff(rangeb)

  width <- diff(rangeplot)
  # y co-ords for labels:
  y_variable <- rangeplot[1] + cpositions[1] * width
  y_nlevel <- rangeplot[1] + cpositions[2] * width
  y_cistring <- rangeplot[1] + cpositions[3] * width
  y_stars <- rangeb[2]
  x_annotate <- seq_along(terms)
  annot_size_mm <- fontsize *
    as.numeric(convertX(unit(theme_get()$text$size, "pt"), "mm"))

  p <- ggplot(data=res, aes(seq_along(terms), exp(coeffs))) +
    geom_rect(aes(xmin = seq_along(terms) - .5, xmax = seq_along(terms) + .5,
                  ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]),
                  fill = ordered(seq_along(terms) %% 2 + 1))) +
    scale_fill_manual(values = c("#FFFFFF33", "#00000033"), guide = "none") +
    geom_point(pch = 15, size = 4) +
    xlab(ylab) +
    ylab(xlab) +
    geom_errorbar(aes(ymin = exp(lower), ymax = exp(upper)), width = 0.15) +
    geom_hline(yintercept = 1, linetype = 3) +
    coord_flip(ylim = exp(rangeplot)) +
    ggtitle(main) +
    scale_y_log10(
      labels = sprintf("%g", breaks),
      expand = c(0.02, 0.02),
      breaks = breaks) +
    theme_light() +
    theme(panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.position = "none",
          panel.border=element_blank(),
          #axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.x=element_text(hjust=0.7),
          plot.title = element_text(hjust = 0.5)) +
    annotate(geom = "text", x = x_annotate, y = exp(y_variable),
             label = terms, fontface = "bold", hjust = 0,
             size = annot_size_mm) +
    annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
             label = toShow$estimate, size = annot_size_mm,
             vjust = ifelse(toShow$estimate == "reference", .5, -0.1)) +
    annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
             label = toShow$ci, size = annot_size_mm,
             vjust = 1.1,  fontface = "italic") +
    annotate(geom = "text", x = x_annotate, y = exp(y_stars),
             label = toShow$stars, size = annot_size_mm,
             hjust = -0.2,  fontface = "italic")

  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"

  report$add_plot(list(id1 = id,
                       id2 = NA,
                       name1 = name,
                       name2 = "",
                       type = "forest"),
                  gt)
}
