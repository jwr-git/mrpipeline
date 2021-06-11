#' Entry point for the pipeline.
#'
#' @param ids1 IDs or filenames for summary statistics
#' @param ids2 IDs or filenames for summary statistics
#' @param p_cutoff P value cutoff for "tophits" from summary statistics,
#'                 default: 5e-8
#' @param f_cutoff F-statistic cutoff for exposure data, default: 10
#' @param chrompos_query Select SNPs within this region only,
#'                       format: "chrom:start-end"
#' @param gene_query Select SNPs nearby this gene only, NOT YET IMPLEMENTED
#' @param markdown Generate markdown report, default: TRUE
#' @param out_path Path to save markdown report, defaults to user home directory
#' @param out_name Name of the analysis given to the markdown report,
#'                 defaults to time and date
#' @param nthreads Number of threads for use in multithreaded functions,
#'                 default: 1
#' @param bfile Plink 1.x format reference files for use in local clumping,
#'              etc, default: NULL
#' @param plink_bin Plink 1.x location, default: NULL
#' @param ... Other arguments for plotting, NOT YET IMPLEMENTED
#'
#' @return list of results for debugging
#' @export
mr_pipeline <- function(ids1, ids2,
                        p_cutoff = 5e-8,
                        f_cutoff = 10,
                        chrompos_query = NULL,
                        gene_query = NULL,
                        markdown = T,
                        out_path = "",
                        out_name = "",
                        nthreads = 1,
                        bfile = NULL,
                        plink_bin = NULL,
                        ...)
{
  if (is.null(p_cutoff) && is.null(chrompos_query) && is.null(gene_query))
  {
    stop("One of the following must be provided: a P value cutoff, a chromosome/position lookup, or a gene lookup. ",
         "Chromosome/positions and genes can be given as a data.frame for multiple lookups.")
  }

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
                       p_cutoff, chrompos_query, gene_query,
                       f_cutoff = f_cutoff,
                       bfile = bfile,
                       plink_bin = plink_bin,
                       nthreads = nthreads)

  # Annotate the IDs
  dat <- annotate_data(dat, id1, id2)

  # MR analyses
  res <- do_mr(dat, report)

  # Colocalisation
  cres <- do_coloc(dat, id1, id2, report)

  # Heterogeneity between instruments
  #hret <- do_heterogeneity(dat, report)

  # Where are we saving the reports?
  if (out_path == "")
  {
    out_path <- Sys.getenv("HOME")
    if (file.access(out_path))
    {
      stop(paste0("No output path given and invalid permissions for User's Home directory: ", out_path, "."),
           "Please supply an output path using the `out_path_ argument.")
    }
    message(paste0("No output path given, therefore defaulting to User's Home directory: ", out_path, "."))
  }

  if (out_name == "" || dir.exists(paste0(out_path, "\\", out_name)))
  {
    out_name <- format(Sys.time(), "%Y-%m-%d-%s")
    message("No output name given or already exists; defaulting to: ", out_name)
  }
  dir.create(file.path(out_path, out_name))
  dir.create(file.path(out_path, out_name, "traits"))

  # Create each per-trait report first
  for (id.exposure in id1$info$id)
  {
    trait.name <- id1$info[id1$info$id == id.exposure, ]$trait
    make_results(out_path, out_name, report, dat[dat$id.exposure == id.exposure, ], id.exposure, trait.name)
  }

  # Per-outcome report
  for (id.outcome in id2$info$id)
  {
    trait.name <- id2$info[id2$info$id == id.outcome, ]$trait
    make_outcomes(out_path, out_name, report, dat[dat$id.outcome == id.outcome, ], id.outcome, trait.name)
  }

  # Main report file
  make_report(out_path, out_name, dat, report)

  return(list(dat, res, cres, report, id1, id2))
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

# Examples
#retval <- mr_pipeline(c("C:\\Users\\jr18055\\Downloads\\IEU-a-2.vcf.gz", "eqtl-a-ENSG00000000457", "eqtl-a-ENSG00000114013", "eqtl-a-ENSG00000116815", "eqtl-a-ENSG00000135541"), c("ieu-b-18", "ieu-b-2"))
#retval <- mr_pipeline_(c("C:\\Users\\jr18055\\Downloads\\IEU-a-2.vcf.gz", "ieu-a-1"), c("ieu-b-85", "ieu-a-1082"), nthreads=4)

# Biogen test
#retval <- mr_pipeline(c("eqtl-a-ENSG00000196296", "eqtl-a-ENSG00000140623", "eqtl-a-ENSG00000067836","eqtl-a-ENSG00000267795", "eqtl-a-ENSG00000061656", "eqtl-a-ENSG00000256762", "eqtl-a-ENSG00000185294", "eqtl-a-ENSG00000120071", "eqtl-a-ENSG00000185829", "eqtl-a-ENSG00000238083", "eqtl-a-ENSG00000228696", "eqtl-a-ENSG00000225190", "eqtl-a-ENSG00000153574", "eqtl-a-ENSG00000203760", "eqtl-a-ENSG00000166037", "eqtl-a-ENSG00000187323", "eqtl-a-ENSG00000101306", "eqtl-a-ENSG00000005059", "eqtl-a-ENSG00000168214"), c("ieu-a-1044", "ieu-a-1045", "ieu-a-1041", "ieu-a-1046", "ieu-a-1047", "ieu-a-1048"))
