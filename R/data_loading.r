# Structure
# 1. Prepare and clean data
# - Function to load and prepare OpenGWAS IDs
# - Function to load and prepare gwasvcf files
# - Annotate data:
# -- ENSGs -> gene names
# -- Gene names -> ENSGs
# -- Diseases -> EFO
# -- EFo -> disease
# --
# 2. Harmonise data
# - F-statistic cutoff into harmonisation
# - Plotting data at this stage?
# 3. MR
# - WR Taylor expansion
# - IVW delta
# - Plots
# 4. Coloc
# - Coloc with SuSiE
# - PWCoCo
# 5. Report

file_to_gwasvcf <- function(file,
                            chr_col,
                            pos_col,
                            nea_col,
                            ea_col,
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
                            verbose = TRUE
                            )
{
  outputs <- parallel::mclapply(1:nrow(file), function(i)
  {
    f <- file[i]

    if (!file.exists(f)) {
      .print_msg(paste0("Could not find file \"", f, "\". Skipping."), verbose)
      next
    }

    dat <- read.table(f, header = header, sep = sep, stringsAsFactors = F)

    if (nrow(dat) < 1) {
      .print_msg(paste0("Failed to load data from file \"", f, "\". Skipping."), verbose)
      next
    }

    out <- paste0(file_path_sans_ext(file), ".vcf")

    dat_to_gwasvcf(dat, out, chr_col, pos_col, nea_col, ea_col,
                  eaf_col = eaf_col, beta_col = beta_col, se_col = se_col,
                  pval_col = pval_col, n = n, n_case = n_case, name = name,
                  verbose = verbse)
  }, mc.cores = cores)

  return(outputs)
}

dat_to_gwasvcf <- function(dat,
                           out,
                           chr_col,
                           pos_col,
                           nea_col,
                           ea_col,
                           eaf_col = NULL,
                           beta_col = NULL,
                           se_col = NULL,
                           pval_col = NULL,
                           n = NULL,
                           n_case = NULL,
                           name = NULL,
                           verbose = TRUE)
{
  if (nrow(dat) < 1) {
    return(NA)
  }

  # Check if these are columns or single values
  # Unsure if this is necessary or whether the gwasvcf::create_vcf can handle
  # this by itself without any pre-cleaning...
  if (n %in% names(dat)) {
    n_col <- n
  } else {
    n_col <- NA
  }

  if (n_case %in% names(dat)) {
    ncase_col <- n_case
  } else {
    ncase_col <- NA
  }

  if (name %in% names(dat)) {
    name_col <- name
  } else {
    name_col <- NA
  }

  attempt <- tryCatch(
    expr = {
      vcf <- dat %>%
        gwasvcf::create_vcf(chrom = chr_col,
                            pos = pos_col,
                            nea = nea_col,
                            ea = ea_col,
                            snp = snp_col,
                            ea_af = eaf_col,
                            effect = beta_col,
                            se = se_col,
                            pval = pval_col,
                            n = ifelse(is.na(n_col), n, n_col),
                            ncase = ifelse(is.na(ncase_col), n_case, ncase_col),
                            name = ifelse(is.na(name_col), name, name_col))
    },
    error = function(e) {
      .print_msg(paste0("Error creating vcf:"), verbose)
      .print_msg(e, verbose)
      return(NA)
    })

  VariantAnnotation::writeVcf(vcf, file = out)
  return(out)
}

.print_msg <- function(msg, verbose)
{
  if (verbose == TRUE) {
    message(msg)
  }
}

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

read_outcome <- function(ids,
                         rsids,
                         proxies = TRUE,
                         proxy_rsq = 0.8,
                         plink = NULL,
                         bfile = NULL,
                         cores = 1,
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

.read_dataset <- function(ids,
                          rsids = NULL,
                          pval = 5e-8,
                          proxies = TRUE,
                          proxy_rsq = 0.8,
                          plink = NULL,
                          bfile = NULL,
                          clump_r2 = 0.01,
                          clump_kb = 10000,
                          pop = "EUR",
                          type = "exposure",
                          cores = 1,
                          verbose = TRUE)
{
  if (Sys.info()['sysname'] == "Windows") {
    .print_msg("Running this pipeline on Windows disables the following features:", verbose)
    .print_msg("Local clumping, proxy searching and other Plink operations; parallelisation.", verbose)
    plink <- NULL
    cores <- 1
  }

  if (is.null(rsids) && is.null(pval)) {
    stop("Please provide either rsIDs or a P value cut-off for extraction.")
  }

  if (is.null(plink) && any(file.exists(ids)) && any(tools::file_ext(ids) == "vcf")) {
    .print_msg("Path to Plink not given; no LD-related functions will be run for local files.", verbose)
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

  local_clump <- !(is.null(plink) && is.null(bfile))

  dat <- parallel::mclapply(1:length(ids), function(i)
  {
    id <- ids[i]

    # Extract data from local GWAS VCF format file
    if (file.exists(id) && tools::file_ext(id) == "vcf")
    {
      .print_msg(paste0("Reading \"", id, "\" as GWAS VCF file."), verbose)

      if (type == "exposure") {
        dat <- gwasvcf::query_gwas(vcf = VariantAnnotation::readVcf(id),
                                   pval = pval)
      } else {
        dat <- gwasvcf::query_gwas(vcf = VariantAnnotation::readVcf(id),
                                   rsid = rsids,
                                   proxies = ifelse(proxies == TRUE, "yes", "no"),
                                   bfile = bfile,
                                   tag_kb = tag_kb,
                                   tag_nsnp = tag_nsnp,
                                   tag_r2 = tag_r2)
      }

      if (!exists("dat") || !length(dat)) {
        warning("Data extraction failed. Please re-check parameters and file paths.")
        return(NULL)
      }

      dat <- dat %>%
        gwasglue::gwasvcf_to_TwoSampleMR(., type = type)

      if (nrow(dat) < 1) {
        warning("No SNPs extracted with the given parameters. Please re-check these and try again.")
        return(NULL)
      }

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
    if (local_clump && type == "exposure")
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
    return(dat)
  }, mc.cores = cores) %>%
    dplyr::bind_rows()

  return(dat)
}

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
                                mart = useEnsembl("snp", dataset = "hsapiens_snp", GRCh = grch))

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

harmonise <- function(exposure,
                      outcome,
                      action = 2)
{
  dat <- TwoSampleMR::harmonise_data(exposure, outcome, action = action)
}
