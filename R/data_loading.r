#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' Convert file(s) to gwasvcf format.
#'
#' Function to convert a file (or files) to gwasvcf format.
#' @seealso [dat_to_gwasvcf()] for converting data.frames
#'
#' @param file Path to file
#' @param chr_col Column name for chromosome
#' @param pos_col Column name for position
#' @param nea_col Column name for non-effect allele
#' @param ea_col Column name for effect allele
#' @param snp_col Column name for SNP (Optional)
#' @param eaf_col Column name for effect allele frequency (Optional)
#' @param beta_col Column name for beta (Optional)
#' @param se_col Column name for standard error (Optional)
#' @param pval_col Column name for P value (Optional)
#'                 NB: P values will be saved as 10^-P
#' @param n Sample size (Optional), can be int or column name
#' @param n_case Number of cases (Optional), can be int or column name
#' @param name Trait name (Optional), can be string or column name
#' @param header Whether file has a header or not (Optional, boolean)
#' @param sep File separater (Optional)
#' @param cores Number of cores for multi-threaded tasks (Optional)
#'              NB: Unavailable on Windows machines
#' @param bcf_tools Path to bcf_tools (Optional)
#' @param verbose Display verbose information (Optional, boolean)
#'
#' @return gwasvcf object(s)
#' @importFrom parallel mclapply
#' @export
file_to_gwasvcf <- function(file,
                            chr_col,
                            pos_col,
                            nea_col,
                            ea_col,
                            snp_col = NULL,
                            eaf_col = NULL,
                            beta_col = NULL,
                            se_col = NULL,
                            pval_col = NULL,
                            n = NULL,
                            n_case = NULL,
                            name = NULL,
                            header = TRUE,
                            sep = "\t",
                            cores = 1,
                            bcf_tools = NULL,
                            verbose = TRUE
                            )
{
  outputs <- parallel::mclapply(1:nrow(file), function(i)
  {
    f <- file[i]

    if (!file.exists(f)) {
      .print_msg(paste0("Could not find file \"", f, "\". Skipping."), verbose)
      return(NA)
    }

    dat <- read.table(f, header = header, sep = sep, stringsAsFactors = F)

    if (nrow(dat) < 1) {
      .print_msg(paste0("Failed to load data from file \"", f, "\". Skipping."), verbose)
      return(NA)
    }

    out <- paste0(file_path_sans_ext(file), ".vcf")

    dat_to_gwasvcf(dat, out, chr_col, pos_col, nea_col, ea_col,
                  snp_col = snp_col, eaf_col = eaf_col, beta_col = beta_col,
                  se_col = se_col, pval_col = pval_col, n = n, n_case = n_case,
                  name = name,
                  bcf_tools = bcf_tools,
                  verbose = verbose)
  }, mc.cores = cores)

  return(outputs)
}

#' Convert data.frame to gwasvcf format.
#'
#' Function to convert a data.frame to gwasvcf format.
#' @seealso [file_to_gwasvcf()] for converting files.
#'
#' @param dat Data.frame
#' @param out Path to save output
#' @param chr_col Column name for chromosome
#' @param pos_col Column name for position
#' @param nea_col Column name for non-effect allele
#' @param ea_col Column name for effect allele
#' @param snp_col Column name for SNP (Optional)
#' @param eaf_col Column name for effect allele frequency (Optional)
#' @param beta_col Column name for beta (Optional)
#' @param se_col Column name for standard error (Optional)
#' @param pval_col Column name for P value (Optional)
#'                 NB: P values will be saved as 10^-P
#' @param n Sample size (Optional), can be int or column name
#' @param n_case Number of cases (Optional), can be int or column name
#' @param name Trait name (Optional), can be string or column name
#' @param bcf_tools Path to bcf_tools (Optional)
#' @param verbose Display verbose information (Optional, boolean)
#'
#' @return gwasvcf object
#' @importFrom gwasvcf create_vcf create_rsidx_index_from_vcf set_bcftools
#' @importFrom VariantAnnotation writeVcf indexVcf
#' @importFrom tools file_path_sans_ext
#' @importFrom Rsamtools bgzip
#' @export
dat_to_gwasvcf <- function(dat,
                           out,
                           chr_col,
                           pos_col,
                           nea_col,
                           ea_col,
                           snp_col = NULL,
                           eaf_col = NULL,
                           beta_col = NULL,
                           se_col = NULL,
                           pval_col = NULL,
                           n = NULL,
                           n_case = NULL,
                           name = NULL,
                           bcf_tools = NULL,
                           verbose = TRUE)
{
  if (nrow(dat) < 1) {
    return(NA)
  }

  # Check if these are columns or single values
  # Unsure if this is necessary or whether the gwasvcf::create_vcf can handle
  # this by itself without any pre-cleaning...
  if (!is.null(n) && n %in% names(dat)) {
    n_col <- n
  } else {
    n_col <- NA
  }

  if (!is.null(n_case) && n_case %in% names(dat)) {
    ncase_col <- n_case
  } else {
    ncase_col <- NA
  }

  if (!is.null(name) && name %in% names(dat)) {
    name_col <- name
  } else {
    name_col <- NA
  }

  # Attempt to set bcf_tools -- could be improved upon no doubt
  if (!is.null(bcf_tools)) {
    gwasvcf::set_bcftools(bcf_tools)
  }

  # Cannot have duplicated rsIDs
  # TODO maybe this could be improved somehow?
  dat <- dat[!duplicated(dat$rsid), ]

#  attempt <- tryCatch(
#    expr = {
      vcf <- gwasvcf::create_vcf(chrom = dat[[chr_col]],
                            pos = dat[[pos_col]],
                            nea = dat[[nea_col]],
                            ea = dat[[ea_col]],
                            snp = dat[[snp_col]],
                            ea_af = dat[[eaf_col]],
                            effect = dat[[beta_col]],
                            se = dat[[se_col]],
                            pval = dat[[pval_col]],
                            n = ifelse(is.na(n_col), n, dat[[n_col]]),
                            ncase = ifelse(is.na(ncase_col), n_case, dat[[ncase_col]]),
                            name = ifelse(is.na(name_col), name, dat[[name_col]]))
#    },
#    error = function(e) {
#      .print_msg(paste0("Error creating vcf:"), verbose)
#      .print_msg(e, verbose)
#      return(NA)
#    })

  VariantAnnotation::writeVcf(vcf, file = out)
  if (gwasvcf::check_bcftools())
  {
    header_string <- "##SAMPLE=<"
    header_string <- paste0(header_string, "TotalVariants=", nrow(dat))
    header_string <- paste0(header_string, ",StudyType=", ifelse(is.na(ncase_col), "Continuous", "CaseControl"))
    header_string <- paste0(header_string, ">")

    # Save to header and write using bcf_tools
    header_path <- tempfile()
    fcon <- file(header_path)
    writeLines(header_string, fcon)
    close(fcon)

    # Temp file
    temp <- paste0(tools::file_path_sans_ext(out), ".temp.vcf")
    cmd <- glue::glue("bcftools annotate -h {header_path} -o {temp} {out}")
    system(cmd)

    # Swap over as bcf_tools will not do this in place
    if (file.exists(temp)) {
      file.remove(out, showWarnings = F)
      file.rename(temp, out)
    }
  }

  # Some file fudgery
  # .vcf -> .vcf.bgz -> index -> rsidx index
  out.bgz <- paste0(out, ".bgz")
  Rsamtools::bgzip(out, overwrite = T)
  file.remove(out, showWarnings = F)
  VariantAnnotation::indexVcf(out.bgz)
  gwasvcf::create_rsidx_index_from_vcf(out.bgz, paste0(tools::file_path_sans_ext(out), ".rsidx"))

  return(out.bgz)
}

#' Helper function for message printing.
#'
#' @param msg Message
#' @param verbose Display message or suppress
#' @keywords Internal
.print_msg <- function(msg, verbose)
{
  if (verbose == TRUE) {
    message(msg)
  }
}

#' Read exposures
#'
#' Read a dataset (or datasets) as exposure.
#' Can accept both OpenGWAS IDs or file paths. Accepts only gwasvcf file formats.
#' This function can clump data locally, if supplied with the `plink` and `bfile`
#' arguments. If these are not supplied, clumping will take place on the
#' OpenGWAS servers \emph{only} for OpenGWAS IDs.
#' If you would like to clump local files, please provide paths to Plink and bfiles.
#'
#' @param ids List of OpenGWAS IDs or file paths (to gwasvcf files)
#' @param pval Threshold to extract SNPs (Optional)
#' @param plink Path to Plink binary (Optional)
#' @param bfile Path to Plink .bed/.bim/.fam files (Optional)
#' @param clump_r2 r2 threshold for clumping SNPs (Optional)
#' @param clump_kb Distance outside of which SNPs are considered in linkage equilibrium (Optional)
#' @param pop Population (Optional, used only for clumping on OpenGWAS)
#' @param cores Number of cores for multi-threaded tasks (Optional)
#'              NB: Unavailable on Windows machines
#' @param verbose Display verbose information (Optional, boolean)
#'
#' @return Data.frame of exposure datasets
#' @export
read_exposure <- function(ids,
                          pval = 5e-8,
                          plink = NULL,
                          bfile = NULL,
                          clump_r2 = 0.01,
                          clump_kb = 10000,
                          pop = "EUR",
                          cores = 1,
                          verbose = TRUE)
{
  exp <- .read_dataset(ids = ids,
                       pval = pval,
                       proxies = FALSE,
                       plink = plink,
                       bfile = bfile,
                       clump_r2 = clump_r2,
                       clump_kb = clump_kb,
                       pop = pop,
                       type = "exposure",
                       cores = cores,
                       verbose = verbose)

  return(exp)
}

#' Read outcomes
#'
#' Read a dataset (or datasets) as outcome.
#' Can accept both OpenGWAS IDs or file paths. Accepts only gwasvcf file formats.
#' This function can search for proxy SNPs locally, if supplied with the `plink`
#' and `bfile` arguments. If these are not supplied, proxy searching will take
#' place on the OpenGWAS servers \emph{only} for OpenGWAS IDs.
#' If you would like to search for proxies locally, please provide paths to Plink
#' and bfiles.
#'
#' @param ids List of OpenGWAS IDs or file paths (to gwasvcf files)
#' @param rsids List of SNP rsIDs to extract
#' @param proxies Whether to search for proxies (Optional, boolean)
#' @param proxy_rsq R2 threshold to use when searching for proxies (Optional)
#' @param proxy_kb kb threshold to use when searching for proxies (Optional)
#' @param proxy_nsnp Number of SNPs when searching for proxies (Optional)
#' @param plink Path to Plink binary (Optional)
#' @param bfile Path to Plink .bed/.bim/.fam files (Optional)
#' @param cores Number of cores for multi-threaded tasks (Optional)
#'              NB: Unavailable on Windows machines
#' @param cores_proxy Number of cores for multi-threaded proxy searching (Optional)
#'                    NB: Unavailable on Windows machines
#'                    NB: Should not be more than `cores` argument!
#' @param verbose Display verbose information (Optional, boolean)
#'
#' @return Data.frame of outcome datasets
#' @export
read_outcome <- function(ids,
                         rsids,
                         proxies = TRUE,
                         proxy_rsq = 0.8,
                         proxy_kb = 5000,
                         proxy_nsnp = 5000,
                         plink = NULL,
                         bfile = NULL,
                         cores = 1,
                         cores_proxy = 1,
                         verbose = TRUE)
{
  out <- .read_dataset(ids = ids,
                       rsids = rsids,
                       proxies = proxies,
                       plink = plink,
                       bfile = bfile,
                       type = "outcome",
                       cores = cores,
                       verbose = verbose)

  return(out)
}

#' Read datasets
#'
#' Helper function that is called from \link[mrpipeline]{read_exposure} and \link[mrpipeline]{read_outcome}.
#' Extracts exposure and outcome data according to arguments.
#' Should not be called directly.
#'
#' @seealso [read_exposure()], [read_outcome()]
#'
#' @param ids List of OpenGWAS IDs or file paths (to gwasvcf files)
#' @param rsids List of SNP rsIDs to extract
#' @param pval Threshold to extract SNPs (Optional)
#' @param proxies Whether to search for proxies (Optional, boolean)
#' @param proxy_rsq R2 threshold to use when searching for proxies (Optional)
#' @param proxy_kb kb threshold to use when searching for proxies (Optional)
#' @param proxy_nsnp Number of SNPs when searching for proxies (Optional)
#' @param plink Path to Plink binary (Optional)
#' @param bfile Path to Plink .bed/.bim/.fam files (Optional)
#' @param clump_r2 r2 threshold for clumping SNPs (Optional)
#' @param clump_kb Distance outside of which SNPs are considered in linkage equilibrium (Optional)
#' @param pop Population (Optional, used only for clumping on OpenGWAS)
#' @param type Type of data (Optional, "exposure" or "outcome")
#' @param cores Number of cores for multi-threaded tasks (Optional)
#'              NB: Unavailable on Windows machines
#' @param cores_proxy Number of cores for multi-threaded proxy searching (Optional)
#'                    NB: Unavailable on Windows machines
#'                    NB: Should not be more than `cores` argument!
#' @param verbose Display verbose information (Optional, boolean)
#'
#' @return Data.frame of datasets
#' @importFrom parallel mclapply
#' @importFrom gwasvcf query_gwas
#' @importFrom gwasglue gwasvcf_to_TwoSampleMR ieugwasr_to_TwoSampleMR
#' @importFrom dplyr mutate bind_rows
#' @importFrom ieugwasr tophits associations gwasinfo ld_clump
#' @importFrom VariantAnnotation readVcf
.read_dataset <- function(ids,
                          rsids = NULL,
                          pval = 5e-8,
                          proxies = TRUE,
                          proxy_rsq = 0.8,
                          proxy_kb = 5000,
                          proxy_nsnp = 5000,
                          plink = NULL,
                          bfile = NULL,
                          clump_r2 = 0.01,
                          clump_kb = 10000,
                          pop = "EUR",
                          type = "exposure",
                          cores = 1,
                          cores_proxy = 1,
                          verbose = TRUE)
{
  if (Sys.info()['sysname'] == "Windows") {
    .print_msg("Running this pipeline on Windows disables the following features:", verbose)
    .print_msg("Local clumping, proxy searching and other Plink operations; parallelisation.", verbose)
    plink <- NULL
    cores <- 1
    cores_proxy <- 1
  }

  if (is.null(rsids) && is.null(pval)) {
    stop("Please provide either rsIDs or a P value cut-off for extraction.")
  }

  if (is.null(plink) && any(file.exists(ids)) && any(tools::file_ext(ids) == "vcf")) {
    .print_msg("Path to Plink not given; LD-related functions will be run through the OpenGWAS API. It may be preferential to run local clumping instead.", verbose)
  }

  if (!is.null(plink) && is.null(bfile)) {
    .print_msg("Path to Plink given but no bfile; please provide this for local LD-related functions.", verbose)
  }

  if (type != "exposure" && type != "outcome") {
    type <- "exposure"
    .print_msg(paste0("Unknown type \"", type, "\" given; defaulting to exposure."), verbose)
  }

  if (type == "outcome" && is.null(rsids)) {
    stop("Type \"outcome\" requires a list of rsIDs to extract.")
  }

  if (!is.null(rsids)) {
    rsids <- unique(rsids)
  }

  if (cores_proxy > cores) {
    cores <- cores / 2 # Halve cores
    cores_proxy <- cores / 4 # And give two cores per to proxying
    .print_msg(paste0("Number of cores to use for proxy searching (", cores_proxy, ") is more than total cores (", cores, "). Reducing."), verbose)
  } else if (cores_proxy > 1) {
    cores <- cores / cores_proxy
  }

  local_clump <- proxies <- !(is.null(plink) && is.null(bfile))

  dat <- parallel::mclapply(1:length(ids), function(i)
  {
    id <- ids[i]

    # Extract data from local GWAS VCF format file
    if (file.exists(id) && (tools::file_ext(id) == "vcf" || tools::file_ext(id) == "bgz" || tools::file_ext(id) == "gz"))
    {
      .print_msg(paste0("Reading \"", id, "\" as GWAS VCF file."), verbose)

      if (type == "exposure") {
        dat <- gwasvcf::query_gwas(vcf = id,
                                   pval = pval)
      } else {
        dat <- gwasvcf::query_gwas(vcf = id,
                                   rsid = rsids,
                                   proxies = ifelse(proxies == TRUE, "yes", "no"),
                                   bfile = bfile,
                                   tag_kb = proxy_kb,
                                   tag_nsnp = proxy_nsnp,
                                   tag_r2 = proxy_rsq,
                                   threads = cores_proxy)
      }

      if (!exists("dat") || !length(dat)) {
        warning("Data extraction failed. Please re-check parameters and file paths.")
        return(NULL)
      }

      if (nrow(dat) < 1) {
        warning("No SNPs extracted with the given parameters. Please re-check these and try again.")
        return(NULL)
      }

      dat <- dat %>%
        gwasglue::gwasvcf_to_TwoSampleMR(., type = type)

      dat[[paste0("id.", type)]] <- id

    } else {
      # Extract data from OpenGWAS DB
      if (type == "exposure") {
        dat <- ieugwasr::tophits(id,
                                 pval = pval,
                                 clump = !local_clump,
                                 r2 = clump_r2,
                                 kb = clump_kb,
                                 pop = pop)
      } else {
        dat <- ieugwasr::associations(rsids,
                                      id,
                                      proxies = proxies,
                                      r2 = proxy_rsq)
      }

      if (!exists("dat") || !length(dat) || inherits(dat, "response")) {
        warning("Data extraction failed for ", id, ". Please re-check parameters and file paths.")
        return(NULL)
      }

      if (nrow(dat) < 1) {
        warning("No SNPs extracted with the given parameters for \"", id, "\".")
        return(NULL)
      }

      dat <- dat %>%
        gwasglue::ieugwasr_to_TwoSampleMR(., type = type)

      if (type == "exposure") {
        dat$exposure <- ieugwasr::gwasinfo(id)$trait
      }
    }

    # Clumping for exposure data
    if (type == "exposure")
    {
      attempt <- try(dat <- dplyr::mutate(dat,
                                          rsid = SNP,
                                          pval = pval.exposure) %>%
                       ieugwasr::ld_clump(clump_kb = clump_kb,
                                          clump_r2 = clump_r2,
                                          pop = pop,
                                          bfile = bfile,
                                          plink_bin = plink),
                     silent = T)

      if (class(attempt) == "try-error" || nrow(dat) < 1) {
        warning("Error with clumping, or no SNPs retained thereafter for exposure: ", id)
        return(NULL)
      }
    }

    dat[[paste0("file.", type)]] <- id

    # Sometimes, ieugwasr returns multiple copies of the same SNP -- who knows why?
    dat <- dat[!duplicated(dat), ]
    return(dat)
  }, mc.cores = cores) %>%
    dplyr::bind_rows()

  return(dat)
}

#' Convert ENSG IDs -> gene names
#'
#' Attempts to convert ENSG IDs to gene names (hgnc_symbol).
#' This is attempting using biomaRt's service and thus requires the optional
#' biomaRt package to be installed.
#'
#' @param dat Data.frame of data
#' @param ensg_col Column name containing ENSG IDs (Optional)
#' @param new_col Column to append to `dat` with converted names (Optional)
#' @param build Genomic build (Optional)
#'
#' @return Data.frame with appended column for names
#' @importFrom biomaRt useMart getBM
#' @keywords Internal
.ensg_to_name <- function(dat,
                          ensg_col = "trait",
                          new_col = "hgnc_symbol",
                          build = "grch37")
{
  if (!require("biomaRt")) {
    warning("biomaRt not found! To annotate ENSG IDs, the biomaRt package must be installed.")
    return(dat)
  }

  if (build == "grch37") {
    biomart <- "ENSEMBL_MART_ENSEMBL"
    host <- "grch37.ensembl.org"
    dataset <- "hsapiens_gene_ensembl"
  } else {
    biomart <- "ensembl"
    host <- "https://www.ensembl.org"
    dataset <- "hsapiens_gene_ensembl"
  }

  attempt <- tryCatch({
    mart.gene <- biomaRt::useMart(biomart = biomart,
                                  host = host,
                                  dataset = dataset)
  }, error = function(e) {
    warning("biomaRt is currently unavailable and so ENSG IDs will not be annotated.")
  })

  if (inherits(attempt, "error")) {
    return(dat)
  }

  results <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position"),
                            filters = "ensembl_gene_id",
                            values = unique(dat[[ensg_col]]),
                            mart = mart.gene)

  if (length(results) && nrow(results)) {
    results <- subset(results, select = c(ensembl_gene_id, hgnc_symbol))
    names(results)[names(results) == "hgnc_symbol"] <- new_col

    dat <- base::merge(dat, results, by.x = ensg_col, by.y = "ensembl_gene_id", all.x = TRUE)
  }

  dat[is.na(dat[[new_col]]), new_col] <- dat[is.na(dat[[new_col]]), ][[ensg_col]]

  return(dat)
}

#' Annotate ENSG IDs with gene names
#'
#' Attempts to convert ENSG IDs to gene names (hgnc_symbol).
#' This is attempting using biomaRt's service and thus requires the optional
#' biomaRt package to be installed.
#'
#' @param dat Data.frame of data
#' @param col Column name containing ENSG IDs (Optional)
#' @param gene_name_col Column to append to `dat` with converted names (Optional)
#' @param build Genomic build (Optional)
#'
#' @return Data.frame with appended column for names
#' @export
annotate_ensg <- function(dat,
                          column = "exposure",
                          gene_name_col = "hgnc_symbol",
                          build = "grch37")
{
  # Attempt to annotate ENSG IDs
  lookup <- dat[grepl("ENSG[0-9]+", dat[[column]]), ][[column]]

  if (length(lookup)) {
    dat <- .ensg_to_name(dat,
                         ensg_col = column,
                         new_col = ifelse(is.null(gene_name_col), column, gene_name_col),
                         build = build)
  }

  return(dat)
}

#' Annotate diseases with EFO IDs
#'
#' Attempts to annotate disease names with EFO IDs using the
#' `epigraphdb` package.
#' Note that the matching is fuzzy and some disease names will have multiple
#' associated EFOs which may differ in definition slightly.
#'
#' @param dat Data.frame of data
#' @param column Column name containing disease names
#'
#' @return Data.frame with appended column for EFO IDs
#' @importFrom epigraphdb ontology_gwas_efo
#' @export
annotate_efo <- function(dat,
                         column = "outcome")
{
  dat$efo_id <- NULL

  lapply(1:nrow(dat), function(i) {
    disease <- iconv(dat[i, column], from = 'UTF-8', to = 'ASCII//TRANSLIT')
    attempt <- try(r <- epigraphdb::ontology_gwas_efo(disease, fuzzy = T, mode = "table"), silent = T)

    if (!inherits(attempt, "try-error") && length(r)) {
      efo <- r[r$r.score > 0.95,][1, ]$efo.id
      efo <- tail(strsplit(efo, "/", fixed = T)[[1]], n = 1)

      dat[i, "efo_id"] <<- efo
    }
  })

  return(dat)
}

#' Annotate cis/trans SNPs
#'
#' Attempts to annotate SNPs as cis or trans depending on their
#' location to the gene coding region.
#' This is achieved using the `biomaRt` R package.
#'
#' @param dat Data.frame of data
#' @param cis_region Cis region definition (Optional, in kb)
#' @param chr_col Column name for chromosome (Optional)
#' @param pos_col Column name for position (Optional)
#' @param snp_col Column name for SNP rsIDs (Optional)
#' @param values_col Column name for gene names or ENSG IDs (Optional)
#'                   NB: Choice must match the `filter` value
#' @param filter How to search for genes in biomaRt, either:
#'               \enumerate{
#'               \item "ensembl_gene_id" for ENSG IDs, or
#'               \item "hgnc_symbol" for gene names
#'               }
#'               (Optional)
#' @param missing "Include" or "drop" SNPs which could not be matched (Optional)
#' @param build Genomic build (Optional)
#'
#' @return Data.frame with appended column for cis/trans status
#' @importFrom biomaRt useMart getBM useEnsembl
#' @export
cis_trans <- function(dat,
                      cis_region = 500000,
                      chr_col = "chr.exposure",
                      pos_col = "position.exposure",
                      snp_col = NULL,
                      values_col = "exposure",
                      filter = "ensembl_gene_id",
                      missing = "include",
                      build = "grch37")
{
  if (!(missing %in% c("include", "drop"))) {
    warning("Must select one of the following options for cis/trans annotation: include, drop.")
    return(dat)
  }

  if (!(filter %in% c("ensembl_gene_id", "hgnc_symbol"))) {
    warning("Must select one of the following options for search: ensembl_gene_id, hgnc_symbol")
    return(dat)
  }

  if (!require("biomaRt")) {
    warning("biomaRt not found! To annotate ENSG IDs, the biomaRt package must be installed.")
    return(dat)
  }

  if ((is.null(snp_col) && (is.null(chr_col) || is.null(pos_col)))) {
    stop("Must give either snp_col with rsIDs or chr_col and pos_col with chromosome positions.")
    return(dat)
  }

  if (build == "grch37") {
    biomart <- "ENSEMBL_MART_ENSEMBL"
    host <- "grch37.ensembl.org"
    dataset <- "hsapiens_gene_ensembl"

    grch <- 37
  } else {
    biomart <- "ensembl"
    host <- "https://www.ensembl.org"
    dataset <- "hsapiens_gene_ensembl"

    grch <- NULL
  }

  attempt <- tryCatch({
    mart.gene <- biomaRt::useMart(biomart = biomart,
                                  host = host,
                                  dataset = dataset)
  }, error = function(e) {
    warning("biomaRt is currently unavailable and so cis/trans status will not be annotated.")
  })

  if (inherits(attempt, "error")) {
    return(dat)
  }

  results <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position"),
                            filters = filter,
                            values = unique(dat[[values_col]]),
                            mart = mart.gene)

  if (length(results) && nrow(results)) {
    dat <- base::merge(dat, results, by.x = values_col, by.y = filter, all.x = TRUE, suffixes = c("", ".y"))
    dat$cistrans <- "U"

    if (!is.null(chr_col) && !is.null(pos_col)) {
      dat$cistrans <- ifelse(
        is.na(dat$chromosome_name) | is.na(dat$start_position) | is.na(dat$end_position),
        "U",
          ifelse(dat[[chr_col]] == dat$chromosome_name & as.numeric(dat[[pos_col]]) >= as.numeric(dat$start_position) - cis_region & as.numeric(dat[[pos_col]]) <= as.numeric(dat$end_position) + cis_region,
            "C",
            "T"
          )
      )
    } else {
      snp_res <- biomaRt::getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end"),
                                filters = "snp_filter",
                                values = unique(dat[[snp_col]]),
                                mart = biomaRt::useEnsembl("snp", dataset = "hsapiens_snp", GRCh = grch))

      if (length(snp_res) && nrow(snp_res)) {
        dat <- base::merge(dat, snp_res, by.x = snp_col, by.y = refsnp_id, all.x = TRUE, suffixes = c("", ".z"))

        dat$cistrans <- ifelse(
          is.na(dat$chr_name) | is.na(dat$chrom_start) | is.na(dat$chrom_end),
          "U",
          ifelse(dat$chr_name == dat$chromosome_name & as.numeric(dat$chrom_start) >= as.numeric(dat$start_position) - cis_region & as.numeric(dat$chrom_start) <= as.numeric(dat$end_position) + cis_region,
                 "C",
                 "T"
          )
        )

        # Clean up SNP merged cols
      }
    }

    # Clean up rest of merged cols
    if ("hgnc_symbol.y" %in% names(dat)) {
      dat <- subset(dat, select = -hgnc_symbol.y)
    }
  }

  return(dat)
}

#' Helper function to extract colocalisation regions for when one dataset comes
#' from a local file and another from OpenGWAS.
#' @seealso [gwasglue::gwasvcf_to_coloc()], [gwasglue::ieugwasr_to_coloc()]
#'
#' @param f1 File path or OpenGWAS ID for trait 1
#' @param f2 File path or OpenGWAS ID for trait 2
#' @param chrpos Character of the format chr:pos1-pos2
#' @param verbose Display verbose information (Optional, boolean)
#'
#' @return list of coloc-ready data
#' @keywords Internal
#' @importFrom VariantAnnotation readVcf header meta
#' @importFrom gwasvcf query_gwas vcf_to_tibble
#' @importFrom ieugwasr associations gwasinfo
.cdat_from_mixed <- function(f1, f2,
                             chrpos,
                             verbose = TRUE)
{
  if (file.exists(f1)) {
    vcf <- f1
    opengwas <- f2
  } else {
    vcf <- f2
    opengwas <- f1
  }

  # First, get data from vcf file
  vcffile <- VariantAnnotation::readVcf(vcf)
  vcf_range <- gwasvcf::query_gwas(vcffile, chrompos = c(chrpos))

  if (length(vcf_range) < 1) {
    .print_msg(paste0("No data extracted from vcf file: \"", vcf, "\" for range: \"", chrpos, "\". Cannot run coloc for this dataset."), verbose = verbose)
    return(NA)
  }

  # Code from gwasglue::gwasvcf_to_coloc()
  tab1 <- vcf_range %>%
    gwasvcf::vcf_to_tibble()

  type1 <- ifelse(VariantAnnotation::header(vcffile) %>%
                    VariantAnnotation::meta() %>%
                    {.[["SAMPLE"]][["StudyType"]]} == "Continuous", "quant", "cc")

  # No header? Try to determine from present data
  if (length(type1) == 0) {
    if ("NC" %in% names(tab1) && !all(is.na(tab1$NC))) {
      type1 <- "cc"
    } else {
      type1 <- "quant"
    }
  }

  tab1$AF[is.na(tab1$AF)] <- 0.5

  # Next, from OpenGWAS
  opengwas_range <- ieugwasr::associations(id = f2, variants = chrpos) %>%
    subset(., !duplicated(rsid))

  if (length(opengwas_range) < 1 || nrow(opengwas_range) < 1) {
    .print_msg(paste0("No data extracted from OpenGWAS ID: \"", opengwas, "\" for range: \"", chrpos, "\". Cannot run coloc for this dataset."), verbose = verbose)
    return(NA)
  }

  # Code from gwasglue::ieugwasr_to_coloc()
  tab2 <- opengwas_range

  tab2$eaf <- as.numeric(tab2$eaf)
  tab2$eaf[which(tab2$eaf > 0.5)] <- 1 - tab2$eaf[which(tab2$eaf > 0.5)]
  s <- sum(is.na(tab2$eaf))
  if(s > 0)
  {
    .print_msg(paste0(s, " out of ", nrow(tab2), " variants have missing allele frequencies in ", opengwas, ". Setting to 0.5"), verbose = verbose)
    tab2$eaf[is.na(tab2$eaf)] <- 0.5
  }

  info2 <- ieugwasr::gwasinfo(opengwas)
  if (is.na(info2$unit)) {
    if (!("ncase" %in% names(info2))) {
      info2$ncase <- NA
    }

    if (is.na(info2$ncase)) {
      type2 <- "quant"
    } else {
      type2 <- "cc"
    }
  } else {
    type2 <- ifelse(info2$unit %in% c("logOR", "log odds"), "cc", "quant")
  }

  tab2$n[is.na(tab2$n)] <- info2$sample_size

  # Get overlap
  # TODO Needs to be made better!
  #index <- as.character(tab1$REF) == as.character(tab2$ea) &
  #  as.character(tab1$ALT) == as.character(tab2$nea) &
  #  as.character(tab1$seqnames) == as.character(tab2$rsid) &
  #  tab1$start == tab2$position
  tab1 <- tab1[tab1$rsid %in% tab2$rsid, ]
  tab2 <- tab2[tab2$rsid %in% tab1$rsid, ]
  tab1 <- tab1[tab1$rsid %in% tab2$rsid, ]

  if (nrow(tab1) < 1 || nrow(tab2) < 1) {
    .print_msg(paste0("No SNPs matched based on rsID between \"", vcf, "\" and \"", opengwas, "\". Skipping coloc analysis for this pair."), verbose = verbose)
    return(NA)
  }

  out1 <- tab1 %>%
    { list(pvalues = 10^-.$LP,
           N = .$SS,
           MAF = .$AF,
           beta = .$ES,
           varbeta = .$SE^2,
           type = type1,
           snp = .$rsid,
           z = .$ES / .$SE,
           chr = .$seqnames,
           pos = .$start,
           id = vcf)
    }

  if(type1 == "cc")
  {
    out1$s <- mean(tab1$NC / tab1$SS, na.rm=TRUE)
  }

  out2 <- tab2 %>%
    { list(pvalues = .$p,
           N = .$n,
           MAF = .$eaf,
           beta = .$beta,
           varbeta = .$se^2,
           type = type2,
           snp = .$rsid,
           z = .$beta / .$se,
           chr = .$chr,
           pos = .$position,
           id = opengwas)
    }

  if(type2 == "cc")
  {
    out2$s <- info2$ncase / info2$sample_size
  }

  return(list(dataset1 = out1, dataset2 = out2))
}

#' Get column names from agnostic but formatted data.frame
#'
#' @param df Data.frame of formatted data (exposure or outcome)
#' @param data Name of column to find
#' @param type Type of data (exposure or outcome); NA to check automatically
#'
#' @keywords Internal
#' @return Name of column formatted for given data.frame
get_col_name <- function(df, data, type = NA)
{
  if (is.na(type)) {
    name <- ifelse(
      "exposure" %in% names(df),
      paste0(data, ".exposure"),
      paste0(data, ".outcome")
    )
  } else {
    name <- paste0(data, ".", type)
  }

  if (!name %in% names(df))
  {
    stop("Could not find column for \"", data, "\" in given data.")
  }

  return(name)
}

#' Check if SNPs are good for use in analyses and mark them as such.
#'
#' @param dat A data.frame of formatted data (exposure or outcome)
#' @param analyses Which analyses should be checked?
#' @param drop Whether to drop SNPs if they failed the check
#'
#' @details List of analyses and what data are checked for:
#' \itemize{
#'  \item{"MR"} {beta, SE, P value}
#'  \item{"coloc"} {chromosome, position, P value}
#' }
#'
#' @return Data.frame
#' @export
check_snps <- function(dat,
                       analyses = c("mr", "coloc"),
                       drop = T)
{
  if (length(dat) == 0 || nrow(dat) < 1) {
    stop("No data were given.")
  }

  if (length(analyses) == 0) {
    warning("No analyses given to check!")
    return(dat)
  }

  dat$include <- T
  # MR needs: beta, se, P
  # MR nice to have: eaf, alleles
  # Coloc needs: chr, pos, P
  if ("mr" %in% analyses || "coloc" %in% analyses)
  {
    dat[is.na(dat[[get_col_name(dat, "pval")]]), "include"] <- F
  }

  if ("mr" %in% analyses || "fstat" %in% analyses)
  {
    dat[is.na(dat[[get_col_name(dat, "beta")]]), "include"] <- F
    dat[is.na(dat[[get_col_name(dat, "se")]]), "include"] <- F
  }

  if ("coloc" %in% analyses)
  {
    dat[is.na(dat[[get_col_name(dat, "chr")]]), "include"] <- F
    dat[is.na(dat[[get_col_name(dat, "position")]]), "include"] <- F
  }

  if (drop)
  {
    dat <- dat[dat$include, ]
  }

  return(dat)
}
