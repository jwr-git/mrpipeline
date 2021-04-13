# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

mr_pipeline <- function(id1, id2,
                         p_cutoff = 5e-8,
                         f_cutoff = 10,
                         chrompos_query = NULL,
                         gene_query = NULL,
                         markdown = T,
                         out_path = "",
                         out_name = "",
                         nthreads = 1,
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
  id1 <- load_ids(id1)
  id2 <- load_ids(id2)

  # Read datasets, format and harmonise in TwoSampleMR format
  dat <- read_datasets(id1, id2, report, p_cutoff, chrompos_query, gene_query, f_cutoff = f_cutoff, nthreads = nthreads)

  # MR analyses
  res <- do_mr(dat, report)

  # Colocalisation
  cres <- do_coloc(id1, id2, dat, report)

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
  for (id.exposure in id1$id)
  {
    make_results(out_path, out_name, report, id.exposure)
  }

  make_report(out_path, out_name, dat, report)

  return(list(dat, res, cres, report))
}

load_ids <- function(ids)
{
  # IDs may be local files/filepaths or simply IDs for ieugwasr
  dat <- dplyr::tibble(id = ids, trait = ids, filename = ids)

  loc <- file.exists(ids)
  if (any(loc))
  {
    dat[loc,]$filename <- file.path(dat[loc,]$filename)
    dat[loc,]$trait <- tools::file_path_sans_ext(basename(dat[loc,]$filename))
    dat[loc,]$id <- tools::file_path_sans_ext(basename(dat[loc,]$filename))
  }
  return(dat)
}

# C:\\Users\\jr18055\\Downloads\\IEU-a-2.vcf.gz
#retval <- mr_pipeline(c("C:\\Users\\jr18055\\Downloads\\IEU-a-2.vcf.gz", "eqtl-a-ENSG00000000457", "eqtl-a-ENSG00000000460", "eqtl-a-ENSG00000134460", "eqtl-a-ENSG00000114013", "eqtl-a-ENSG00000116815", "eqtl-a-ENSG00000135541", "eqtl-a-ENSG00000164691"), c("ieu-b-18", "ieu-b-2"))
#retval <- mr_pipeline_(c("C:\\Users\\jr18055\\Downloads\\IEU-a-2.vcf.gz", "ieu-a-1"), c("ieu-b-85", "ieu-a-1082"), nthreads=4)
