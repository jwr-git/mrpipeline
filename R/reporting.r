#' A Reference Class for the reporting of QC, results, etc.
#'
#' @field exposures A data.frame of QC results for the exposure datasets
#' @field outcomes A data.frame of QC results for the outcome datasets
#' @field results A data.frame of MR results
#' @field sresults A data.frame of MR sensitivity results
#' @field hetresults A data.frame of heterogeneity results
#' @field pleioresults A data.frame of pleiotropy results
#' @field steigerresults A data.frame of Steiger filtering results
#' @field cresults A data.frame of colocalisation results
#' @field bonferroni Numeric, Bonferroni-corrected P value threshold
#' @field plots A data.frame of plot names, IDs, etc.
#' @field raw.plots A list of all plots
QCReport <- setRefClass("QCReport",
                        fields = list(exposures = "data.frame",
                                      outcomes = "data.frame",
                                      results = "data.frame",
                                      sresults = "data.frame",
                                      hetresults = "data.frame",
                                      pleioresults = "data.frame",
                                      steigerresults = "data.frame",
                                      cresults = "data.frame",
                                      otresults = "data.frame",
                                      dgidbresults = "data.frame",
                                      dgidbdrugs = "data.frame",
                                      bonferroni = "numeric",
                                      plots = "data.frame",
                                      raw.plots = "list"),
                        methods = list(

                          initialise = function() {
                            "Initialise all fields"
                            exposures <<- data.frame(Trait.Name = character(),
                                                    Trait.ID = character(),
                                                    Trait.Annotated = character(),
                                                    SNPs = integer(),
                                                    Preclumped = logical(),
                                                    Rem.F.Stat = integer(),
                                                    Rem.Clumped = integer(),
                                                    SNPs.Formatted = integer(),
                                                    From.IEUGWASDB = integer(),
                                                    stringsAsFactors = F
                            )
                            outcomes <<- data.frame(Trait.Name = character(),
                                                    Trait.ID = character(),
                                                    Trait.Annotated = character(),
                                                    SNPs = integer(),
                                                    Num.Proxies = integer(),
                                                    SNPs.Formatted = integer(),
                                                    From.IEUGWASDB = integer(),
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
                            hetresults <<- data.frame(id.exposure = character(),
                                                      id.outcome = character(),
                                                      exposure = character(),
                                                      outcome = character(),
                                                      method = character(),
                                                      Q = integer(),
                                                      Q_df = integer(),
                                                      Q_pval = integer()
                            )
                            pleioresults <<- data.frame(id.exposure = character(),
                                                        id.outcome = character(),
                                                        outcome = character(),
                                                        exposure = character(),
                                                        method = character(),
                                                        egger_intercept = character(),
                                                        se = character(),
                                                        pval = character()
                            )
                            steigerresults <<- data.frame(id.exposure = character(),
                                                          id.outcome = character(),
                                                          exposure = character(),
                                                          outcome = character(),
                                                          snp_r2.exposure = integer(),
                                                          snp_r2.outcome = integer(),
                                                          correct_causal_direction = logical(),
                                                          steiger_pval = integer(),
                                                          flag = character()
                            )
                            cresults <<- data.frame(id1 = character(),
                                                    id2 = character(),
                                                    name1 = character(),
                                                    name2 = character(),
                                                    nsnps = integer(),
                                                    H0 = integer(),
                                                    H1 = integer(),
                                                    H2 = integer(),
                                                    H3 = integer(),
                                                    H4 = integer(),
                                                    chrpos = character(),
                                                    plot = character()
                            )
                            otresults <<- data.frame(id1 = character(),
                                                     id2 = character(),
                                                     name1 = character(),
                                                     name2 = character(),
                                                     overall.ot = integer(),
                                                     literature.ot = integer(),
                                                     rna_expression.ot = integer(),
                                                     genetic_assoc.ot = integer(),
                                                     somatic_mute.ot = integer(),
                                                     known_drug.ot = integer(),
                                                     animal_model.ot = integer(),
                                                     affected_pathway.ot = integer()
                            )
                            dgidbresults <<- data.frame(id = character(),
                                                        trait = character(),
                                                        Druggable.Genome = logical(),
                                                        Clinically.Actionable = logical(),
                                                        Drug.Resistant = logical()
                            )
                            dgidbdrugs <<- data.frame(Interaction.Types = character(),
                                                      Drug.Name = character(),
                                                      CHEMBL.ID = character(),
                                                      Sources = character(),
                                                      PMIDs = character(),
                                                      Score = numeric(),
                                                      trait = character(),
                                                      id = character()
                            )
                            plots <<- data.frame(id1 = character(),
                                                 id2 = character(),
                                                 name1 = character(),
                                                 name2 = character(),
                                                 type = character())
                            raw.plots <<- list()
                          },

                          get_exposures_name = function() {
                            "Internal function to get name of exposure data.frame"
                            return("exposures")
                          },

                          get_outcomes_name = function() {
                            "Internal function to get name of outcome data.frame"
                            return("outcomes")
                          },

                          add_dataset = function(name, x) {
                            "Add dataset to either exposure/outcome data.frame"
                            temp <- base::get(name)
                            temp[nrow(temp) + 1, ] <- c(x, rep(0, ncol(temp) - length(x)))
                            assign(name, temp, inherits = T)
                          },

                          update_cell = function(name, trait, col, val) {
                            "Update cell of exposure/outcome data.frame"
                            temp <- base::get(name)
                            temp[temp$Trait.Name == trait, col] <- val
                            assign(name, temp, inherits = T)
                          },

                          add_plot = function(p, raw) {
                            "Adds information and raw plot to report"
                            raw.plots[[nrow(plots) + 1]] <<- raw
                            plots[nrow(plots) + 1, ] <<- p
                          },

                          add_results = function(r) {
                            "Adds MR results to data.frame"
                            results <<- rbind(results, r)
                            #results[nrow(results) + 1, ] <<- r
                          },

                          add_sresults = function(name, x) {
                            "Adds sensitivity MR results to data.frame"
                            temp <- base::get(name)
                            temp <- rbind(temp, x)
                            assign(name, temp, inherits = T)
                          },

                          add_cresults = function(r) {
                            "Adds colocalisation results to data.frame"
                            cresults <<- rbind(cresults, r)
                          },
                          add_hresults = function(r) {
                            "Adds heterogeneity results to data.frame"
                            hetresults <<- rbind(hetresults, r)
                          },
                          add_otresults = function(r) {
                            "Adds Open Targets results to data.frame"
                            otresults <<- rbind(otresults, r)
                          },
                          add_dgidbresults = function(r) {
                            "Adds DGIdb results to data.frame"
                            dgidbresults <<- rbind(dgidbresults, r)
                          },
                          add_dgidbdrugs = function(r) {
                            "Adds drugs from DGIdb to data.frame"
                            dgidbdrugs <<- rbind(dgidbdrugs, r)
                          }
                        )
                        )

#' Make "main" markdown report
#'
#' @param filepath Character, path to save
#' @param filename Character, name of file
#' @param dat A data.frame of harmonised data
#' @param report QCReport class of results, etc. for reporting
#'
#' @return NULL
make_report <- function(filepath, filename, dat, report, conf)
{
  rmarkdown::render(input = "report/report.rmd",
                    output_dir = paste0(filepath, "\\", filename),
                    output_file = "report.html",
                    params = list(report = report,
                                  conf = conf,
                                  dat = dat,
                                  name = filename)
                    )
}

#' Make results markdown report for each exposure
#'
#' @param filepath Character, path to save
#' @param filename Character, name of file
#' @param report QCReport class of results, etc. for reporting
#' @param dat A data.frame of harmonised data
#' @param id.exposure ID of exposure
#' @param trait.name Annotated name of exposure
#'
#' @return NULL
make_results <- function(filepath, filename, report, conf, dat, id.exposure, trait.name)
{
  rmarkdown::render(input = "report/results.rmd",
                    output_dir = paste0(filepath, "\\", filename, "\\traits"),
                    output_file = paste0(trait.name, ".html"),
                    params = list(report = report,
                                  conf = conf,
                                  dat = dat,
                                  id.exposure = id.exposure,
                                  trait.name = trait.name,
                                  new_title = paste0(trait.name, " results"))
                    )
}

#' Make markdown report for each outcome
#'
#' @param filepath Character, path to save
#' @param filename Character, name of file
#' @param report QCReport class of results, etc. for reporting
#' @param dat A data.frame of harmonised data
#' @param id.outcome ID of outcome
#' @param trait.name Annotated name of exposure
#'
#' @return NULL
make_outcomes <- function(filepath, filename, report, dat, id.outcome, trait.name)
{
  rmarkdown::render(input = "report/outcomes.rmd",
                    output_dir = paste0(filepath, "\\", filename, "\\traits"),
                    output_file = paste0(trait.name, ".html"),
                    params = list(report = report,
                                  dat = dat,
                                  id.outcome = id.outcome,
                                  trait.name = trait.name,
                                  new_title = paste0(trait.name, " results"))
                    )
}
