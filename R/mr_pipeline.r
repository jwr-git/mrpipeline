#' Entry point for the pipeline.
#'
#' @param ids1 IDs or filenames for summary statistics
#' @param ids2 IDs or filenames for summary statistics
#' @param out_name Name of the analysis given to the markdown report,
#'                 defaults to time and date
#' @param config_file Path to config.yml file. Defaults to file that comes with the
#'                    package. Please see that file for more details.
#' @param ... Other arguments for plotting, NOT YET IMPLEMENTED
#'
#' @return list of results for debugging
#' @export
mr_pipeline <- function(ids1, ids2,
                        out_name = "",
                        config_file = "",
                        ...)
{
  # Config file loading
  if (config_file != "") {
    # Attempt to load config from given path
    conf <- try(config::get(file = config_file), silent = T)
  }

  if (exists("conf") == F) {
    message("Using default config file.")
    conf <- try(config::get(system.file("config", "config.yml", package = "mrpipeline")), silent = T)
  }

  if (class(conf) == "try-error" | exists("conf") == F) {
    stop("No config file found -- has the default been moved or deleted?")
  }

  conf <- validate_config(conf)

  #if (is.null(p_cutoff) && is.null(chrompos_query) && is.null(gene_query))
  #{
  #  stop("One of the following must be provided: a P value cutoff, a chromosome/position lookup, or a gene lookup. ",
  #       "Chromosome/positions and genes can be given as a data.frame for multiple lookups.")
  #}

  # Set up report class containing meta-data
  report <- QCReport$new()
  report$initFields()
  report$initialise()

  # Parsing ids/filename given
  id1 <- DatasetIDs$new()
  id1$initFields()
  id1$initialise(ids1)

  id2 <- DatasetIDs$new()
  id2$initFields()
  id2$initialise(ids2)

  # Read datasets, format and harmonise in TwoSampleMR format
  dat <- read_datasets(id1, id2,
                       report,
                       conf)

  # Annotate the IDs
  dat <- annotate_data(dat, id1, id2)

  # MR analyses
  res <- do_mr(dat, report, conf)

  # Colocalisation
  cres <- do_coloc(dat, id1, id2, report, conf)

  # Heterogeneity between instruments
  #hret <- do_heterogeneity(dat, report)

  # Drug target-related
  #ot_res <- ot_linkage(id1, id2, report)
  dgidb_linkage(id1, id2, report)
  #clinvar_linkage(id1, id2, report)
  epigdb_linkage(id1, id2, report)
  intermine(id1, id2, report)

  # Where are we saving the reports?
  out_path <- conf$out_path

  if (out_name == "" || dir.exists(paste0(out_path, "\\", out_name)))
  {
    out_name <- conf$out_name
    message("No output name given or already exists; defaulting to: ", out_name)
  }
  dir.create(file.path(out_path, out_name))
  dir.create(file.path(out_path, out_name, "traits"))

  # Create each per-trait report first
  for (id.exposure in id1$info$id)
  {
    trait.name <- id1$info[id1$info$id == id.exposure, ]$trait
    #make_results(out_path, out_name, report, conf, dat[dat$id.exposure == id.exposure, ], id.exposure, trait.name)
  }

  # Per-outcome report
  for (id.outcome in id2$info$id)
  {
    trait.name <- id2$info[id2$info$id == id.outcome, ]$trait
    #make_outcomes(out_path, out_name, report, dat[dat$id.outcome == id.outcome, ], id.outcome, trait.name)
  }

  # Main report file
  #make_report(out_path, out_name, dat, report, conf)

  return(list(dat, res, cres, report, id1, id2, conf))
}

#' Annotates the data using given IDs
#'
#' @param dat Data.frame of data from vcf files or OpenGWAS DB
#' @param id1 DatasetsID class of exposure IDs
#' @param id2 DatasetsID class of outcome IDs
#'
#' @return Data.frame of annotated dat
annotate_data <- function(dat, id1, id2)
{
  dat$ogdb_eid <- dat$id.exposure
  dat$ogdb_oid <- dat$id.outcome

  dat <- left_join(dat, id1$info, by=c("id.exposure"="id")) %>%
    mutate(exposure = ifelse(is.na(trait), exposure, trait)) %>%
    select(-filename, -trait)

  dat <- left_join(dat, id2$info, by=c("id.outcome"="id")) %>%
    mutate(outcome = ifelse(is.na(trait), outcome, trait)) %>%
    select(-filename, -trait)

  return(dat)
}

#' Validates parameters in config file
#'
#' @param conf config::config file of parameters
#'
#' @return Validated config class
validate_config <- function(conf)
{
  if (file.access(conf$out_path)) {
    stop(paste0("Cannot access output path given or invalid permissions directory: ", conf$out_path, "."))
  }

  if (class(conf$cores) != "numeric") {
    warning("Cores parameter is not numeric or not set in environment; defaulting to 1.")
    conf$cores <- 1
  }

  if (conf$pop != "EUR" & conf$pop != "SAS" & conf$pop != "EAS" & conf$pop != "AFR" & conf$pop != "AMR") {
    warning("Invalid population given, valid options are: EUR, SAS, EAS, AFR, AMR. Defaulting to EUR.")
    conf$pop <- "EUR"
  }

  # LD-related
  if (Sys.info()['sysname'] == "Windows"
      & (!is.null(conf$dbfile_path)
         | !is.null(conf$bfile_path)
         | !is.null(conf$plink_path)
         | !is.null(conf$bcf_path))
      )
    {
    message("Many local LD-related functions require OSX or Unix to use. These have been disabled as you seem to be using Windows.")
    conf$bfile_path <- NULL
    conf$dbfile_path <- NULL
    conf$plink_path <- NULL
    conf$bcf_path <- NULL
  }

  if (!is.null(conf$dbfile_path)) {
    if (!file.exists(conf$dbfile_path)) {
      warning("Cannot access LD tagging database at: ", conf$dbfile_path, ".")
      conf$dbfile_path <- NULL
    }
  }

  if (!is.null(conf$bfile_path)) {
    if (!file.exists(paste0(conf$bfile_path), ".bim")) {
      warning("Cannot access reference panel at: ", conf$bfile_path, ". Did you include the file endings? Local LD panel will not be used for any of these analyses.")
      conf$bfile_path <- NULL
      conf$plink_path <- NULL
    }
  }

  if (!is.null(conf$plink_path)) {
    if (!file.exists(paste0(conf$plink_path))) {
      warning("Cannot access Plink at: ", conf$plink_path, ". Plink will not be used for any of these analyses.")
      conf$bfile_path <- NULL
      conf$plink_path <- NULL
    }
  }

  if (!is.null(conf$bcf_path)) {
    if (conf$bcf_path == "") {
      a <- requireNamespace("genetics.binaRies")
      if (a) {
        message("Path not provided, using binaries in the explodecomputer/genetics.binaRies package")
        conf$bcftools <- genetics.binaRies::get_bcftools_binary()
      }
      else {
        message("Provide a path to bcftools binary or run devtools::install_github('explodecomputer/genetics.binaRies')")
      }
    }
  }

  if (class(conf$proxies) != "NULL" & conf$proxies != 0) {
    conf$proxies <- TRUE
  }
  else {
    conf$proxies <- FALSE
  }

  return(conf)
}

# Examples
#retval <- mr_pipeline(c("C:\\Users\\jr18055\\Downloads\\IEU-a-2.vcf.gz", "eqtl-a-ENSG00000000457", "eqtl-a-ENSG00000114013", "eqtl-a-ENSG00000116815", "eqtl-a-ENSG00000135541"), c("ieu-b-18", "ieu-b-2"))
#retval <- mr_pipeline_(c("C:\\Users\\jr18055\\Downloads\\IEU-a-2.vcf.gz", "ieu-a-1"), c("ieu-b-85", "ieu-a-1082"), nthreads=4)

# Biogen test
#retval <- mr_pipeline(c("eqtl-a-ENSG00000196296", "eqtl-a-ENSG00000140623", "eqtl-a-ENSG00000067836","eqtl-a-ENSG00000267795", "eqtl-a-ENSG00000061656", "eqtl-a-ENSG00000256762", "eqtl-a-ENSG00000185294", "eqtl-a-ENSG00000120071", "eqtl-a-ENSG00000185829", "eqtl-a-ENSG00000238083", "eqtl-a-ENSG00000228696", "eqtl-a-ENSG00000225190", "eqtl-a-ENSG00000153574", "eqtl-a-ENSG00000203760", "eqtl-a-ENSG00000166037", "eqtl-a-ENSG00000187323", "eqtl-a-ENSG00000101306", "eqtl-a-ENSG00000005059", "eqtl-a-ENSG00000168214"), c("ieu-a-1044", "ieu-a-1045", "ieu-a-1041", "ieu-a-1046", "ieu-a-1047", "ieu-a-1048"))

# Denis stuff
#retval <- mr_pipeline(c("eqtl-a-ENSG00000254858", "eqtl-a-ENSG00000215912", "eqtl-a-ENSG00000173226", "eqtl-a-ENSG00000162695", "eqtl-a-ENSG00000159640", "eqtl-a-ENSG00000203710"), c("ieu-b-18", "ieu-b-2")))
