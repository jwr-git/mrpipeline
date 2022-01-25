#' Z plot
#'
#' Creates a Z-score plot for SNPs, where \eqn{Z = b / SE}.
#' This should follow a parabollic shape and so can be used to find certain
#' SNPs which may not follow this shape.
#' If `plotly` is installed, the plot will be returned interactive.
#'
#' @param dat A data.frame of data
#' @param snp_col Column of SNP names (Optional)
#' @param beta_col Column of MR beta estimates (Optional)
#' @param se_col Column of standard errors for the beta estimates (Optional)
#' @param pval_col Column of P values (Optional)
#' @param force_static True for forcing the plot to be returned as a static plot (Optional)
#'
#' @return Plot
#' @export
z_plot <- function(dat,
                   snp_col = "SNP",
                   beta_col = "beta.exposure",
                   se_col = "se.exposure",
                   pval_col = "pval.exposure",
                   force_static = FALSE)
{
  if (nrow(dat) < 1 || length(dat) < 1) {
    warning("\"dat\" is empty or lacks data. Please check before continuing.")
    return(NA)
  }

  if (!require("plotly") && !force_static) {
    warning("plotly not found! Plot will be static. Please install plotly for interactive plots.")
  }

  if (!(snp_col %in% names(dat))
      || !(beta_col %in% names(dat))
      || !(se_col %in% names(dat))
      || !(pval_col %in% names(dat))) {
    warning("Could not find SNP, beta, SE or P value column in dat. Please check this before continuing.")
    return(NA)
  }

  p <- ggplot(data = dat, aes(x = get(beta_col) / get(se_col),
                              y = -log10(get(pval_col)))) +
    geom_point(size = 2, alpha = 0.8, aes(text = sprintf("SNP: %s", get(snp_col)))) +
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

  if (require("plotly") && !force_static) {
    p <- plotly::ggplotly(p)
  }

  return(p)
}

#' Volcano plot
#'
#' Creates a volcano plot of Wald ratios from [do_mr()].
#' This function will take all of the Wald ratios in the given data.frame and plot these.
#' If the plot is too crowded, subsetting the results before passing them to this plotter will help.
#' If `plotly` is installed, the plot will be returned interactive.
#'
#' @param res A data.frame of MR results
#' @param label Column whose values will be used to group results by (Optional)
#' @param snp_col Column name for SNPs (Optional)
#' @param beta_col Column name for beta (Optional)
#' @param pval_col Column name for P value (Optional)
#' @param or_col Column name for odds ratio (Optional)
#' @param or_lci_col Column name for lower CI of OR (Optional)
#' @param or_uci_col Column name for upper CI of OR (Optional)
#' @param method_col Column name which contains the MR method (Optional)
#' @param force_static True for forcing the plot to be returned as a static plot (Optional)
#'
#' @return Plot
#' @export
volcano_plot <- function(res,
                         label = "outcome",
                         snp_col = "snp",
                         beta_col = "b",
                         pval_col = "pval",
                         or_col = "or",
                         or_lci_col = "or_lci95",
                         or_uci_col = "or_uci95",
                         method_col = "method",
                         force_static = FALSE)
{
  if (nrow(res) < 1 || length(res) < 1) {
    warning("\"res\" is empty or lacks data. Please check before continuing.")
    return(NA)
  }

  res <- res[res[[method_col]] == "Wald ratio", ]
  if (nrow(res) < 1 || length(res) < 1) {
    warning("\"res\" is empty or lacks data after attempting to find Wald ratios. Please check before continuing.")
    return(NA)
  }

  if (require(plotly) && !force_static) {
    p <- res %>%
      plotly::highlight_key(~res[[label]]) %>%
      plotly::plot_ly(
        x = ~res[[beta_col]], y = ~-log10(res[[pval_col]]),
        type = 'scatter',
        mode = "markers+text",
        textposition = "top",
        hovertemplate = paste(
          "<b>", res[[snp_col]], " : ", res[[label]], "</b><br><br>",
          "MR Beta = %{x:.2f}<br>",
          "OR (95% CI) = ", signif(res[[or_col]], 3), "(", signif(res[[or_lci_col]], 3), ",", signif(res[[or_uci_col]], 3), ")"
        )
      ) %>%
      plotly::highlight(on = "plotly_click", selectize = TRUE, dynamic = TRUE) %>%
      plotly::layout(title = "Volcano plot of Wald ratios",
                     xaxis = list(title = "Beta estimate"),
                     yaxis = list(title = "-log10(P)"))
  } else {
    p <- ggplot2::ggplot(data = res, aes(x = get(beta_col), y = -log10(get(pval_col)))) +
      ggplot2::geom_point(size = 2, alpha = 0.8, aes(color = get(label))) +
      ggplot2::scale_colour_viridis_d("Outcomes") +
      ggplot2::scale_fill_viridis_d("Outcomes") +
      ggplot2::theme_bw() +
      ggplot2::xlab("Beta") +
      ggplot2::ylab("-log10(P value)") +
      ggplot2::ggtitle("Volcano plot of Wald ratios")
  }

  return(p)
}

#' Forest plot
#'
#' Creates a forest plot of MR estimates from do_mr().
#' Will plot both the Wald ratios for all SNPs which form the instrument and
#' inverse variance weighted method. However, if you wish for only the "discovery"
#' results to be plotted (i.e. WR for single-SNP instruments and only IVW for
#' multi-SNP instruments), then setting `plot_all_res` to \code{FALSE} will achieve this.
#' If the plot is too crowded, subsetting the results before passing them to this plotter will help.
#'
#' @seealso [do_mr()]
#'
#' @param res A data.frame of MR results
#' @param snp_col Column name for SNPs (Optional)
#' @param beta_col Column name for beta (Optional)
#' @param se_col Column name for standard error (Optional)
#' @param pval_col Column name for P value (Optional, unused for now)
#' @param or_col Column name for odds ratio (Optional)
#' @param or_lci_col Column name for lower CI of OR (Optional)
#' @param or_uci_col Column name for upper CI of OR (Optional)
#' @param method_col Column name which contains the MR method (Optional)
#' @param exposure_col Column name for exposure names (Optional)
#' @param outcome_col Column name for outcome names (Optional)
#' @param plot_all_res For multi-SNP instruments, also plot the Wald ratios for
#'                     all SNPs (\code{TRUE}) or just the inverse variance weighted
#'                     result (\code{FALSE}).
#'
#' @return Plot
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbarh geom_vline facet_grid xlab ylab theme_bw
#' @export
forest_plot <- function(res,
                        snp_col = "snp",
                        beta_col = "b",
                        se_col = "se",
                        pval_col = NULL,
                        or_col = "or",
                        or_lci_col = "or_lci95",
                        or_uci_col = "or_uci95",
                        method_col = "method",
                        exposure_col = "exposure",
                        outcome_col = "outcome",
                        plot_all_res = TRUE)
{
  if (nrow(res) < 1 || length(res) < 1) {
    warning("\"res\" is empty or lacks data. Please check before continuing.")
    return(NA)
  }

  # If !plot_all_res, then we plot only IVW or WR if only one result
  # Otherwise, all results are plotted in the "summary" style
  if (!plot_all_res) {
    to_plot <- res %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c(exposure_col, outcome_col)))) %>%
      dplyr::mutate(.n = n()) %>%
      dplyr::mutate(.plot = ifelse(.n ==1 | (!!sym(method_col)) == "Inverse variance weighted", TRUE, FALSE))

    if (all(to_plot$.plot== FALSE)) {
      warning("No data retained from using plot_all_res, consider changing this or checking your result.")
      return(NA)
    }

    # Some aesthetics
    to_plot$label <- to_plot[[snp_col]]
    to_plot <- to_plot[to_plot$.plot == TRUE, ]
    to_plot$shape <- 1
  } else {
    to_plot <- res

    # Some aesthetics
    to_plot$label <- ifelse(to_plot[[method_col]] == "Inverse variance weighted",
                            "Inverse variance weighted",
                            to_plot[[snp_col]])
    to_plot$shape <- ifelse(to_plot[[method_col]] == "Inverse variance weighted",
                            5,
                            1)
  }

  p <- to_plot %>%
    dplyr::arrange((!!exposure_col), (!!outcome_col), desc(label)) %>%
    ggplot(data = .,
           aes(x = get(or_col),
               y = label,
               group = interaction(get(exposure_col), get(outcome_col)))) +
    geom_point(aes(shape = as.factor(shape))) +
    geom_errorbarh(aes(xmin = get(or_lci_col), xmax = get(or_uci_col))) +
    geom_vline(xintercept = 1, lty = 2) +
    facet_grid(as.formula(paste(outcome_col, "+", exposure_col, "~.")), scales = "free", space = "free") +
    xlab("Odds Ratio (95% CI)") +
    ylab(" ") +
    theme_bw()

  return(p)
}

# Testing
forest_plot_single <- function()
{
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
  g <- grid::grid.draw(gt)
  return(g)
}
