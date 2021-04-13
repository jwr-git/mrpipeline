QCReport <- setRefClass("QCReport",
                        fields = list(exposures = "data.frame",
                                      outcomes = "data.frame",
                                      results = "data.frame",
                                      sresults = "data.frame",
                                      cresults = "data.frame",
                                      bonferroni = "numeric",
                                      plots = "data.frame",
                                      raw.plots = "list"),
                        methods = list(

                          initialise = function() {
                            exposures <<- data.frame(Trait.Name = character(),
                                                    Trait.ID = character(),
                                                    SNPs = integer(),
                                                    Preclumped = logical(),
                                                    Rem.F.Stat = integer(),
                                                    Rem.Clumped = integer(),
                                                    SNPs.Formatted = integer(),
                                                    stringsAsFactors = F
                            )
                            outcomes <<- data.frame(Trait.Name = character(),
                                                    Trait.ID = character(),
                                                    SNPs = integer(),
                                                    Num.Proxies = integer(),
                                                    SNPs.Formatted = integer(),
                                                    stringsAsFactors = F
                            )
                            results <<- data.frame(id.exposure = character(),
                                                   id.outcome = character(),
                                                   outcome = character(),
                                                   exposure = character(),
                                                   method = character(),
                                                   nsnp = integer(),
                                                   pval = integer(),
                                                   or = integer(),
                                                   or_lci95 = integer(),
                                                   or_uci95 = integer()
                            )
                            sresults <<- data.frame(results)
                            cresults <<- data.frame(id1 = character(),
                                                    id2 = character(),
                                                    nsnps = integer(),
                                                    H0 = integer(),
                                                    H1 = integer(),
                                                    H2 = integer(),
                                                    H3 = integer(),
                                                    H4 = integer(),
                                                    chrpos = character())
                            plots <<- data.frame(id1 = character(),
                                                 id2 = character(),
                                                 type = character())
                            raw.plots <<- list()
                          },

                          get_exposures_name = function() {
                            return("exposures")
                          },

                          get_outcomes_name = function() {
                            return("outcomes")
                          },

                          add_dataset = function(name, x) {
                            temp <- get(name)
                            temp[nrow(temp) + 1, ] <- c(x, rep(0, ncol(temp) - length(x)))
                            assign(name, temp, inherits = T)
                          },

                          update_cell = function(name, trait, col, val) {
                            temp <- get(name)
                            temp[temp$Trait.Name == trait, col] <- val
                            assign(name, temp, inherits = T)
                          },

                          add_plot = function(p, raw) {
                            raw.plots[[nrow(plots) + 1]] <<- raw
                            plots[nrow(plots) + 1, ] <<- p
                          },

                          add_results = function(r) {
                            results <<- rbind(results, r)
                            #results[nrow(results) + 1, ] <<- r
                          },

                          add_sresults = function(r) {
                            sresults <<- rbind(results, r)
                          },

                          add_cresults = function(r) {
                            cresults <<- rbind(cresults, r)
                          }
                        )
                        )

qc_cleaning <- function(dat, report,
                        beta_col = "beta",
                        se_col = "se",
                        pval_col = "pval",
                        trait_col = "trait")
{
  # First we plot b/se against P value in a "volcano" plot
  pipeline_volcano_plot(dat, report,
                        beta_col = beta_col,
                        se_col = se_col,
                        pval_col = pval_col,
                        trait_col = trait_col)

  return(dat)
}

make_report <- function(filepath, filename, dat, report)
{
  rmarkdown::render(input = "report/report.rmd",
                    output_dir = paste0(filepath, "\\", filename),
                    output_file = "report.html",
                    params = list(report = report,
                                  dat = dat,
                                  name = filename)
                    )
}

make_results <- function(filepath, filename, report, id.exposure)
{
  rmarkdown::render(input = "report/results.rmd",
                    output_dir = paste0(filepath, "\\", filename, "\\traits"),
                    output_file = paste0(id.exposure, ".html"),
                    params = list(report = report,
                                  id.exposure = id.exposure,
                                  new_title = paste0(id.exposure, " results"))
                    )
}
