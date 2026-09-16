# prepare_data.R
# Refresh package reference data from the current MatriSpace beta sources.
# Run this from the matrispace package directory.

stapp_data_dir <- "../STapp_beta/data"

required_files <- c(
  "matrisome_v2.rds",
  "ecm_ucell_signatures.rds",
  "ultimate_ecm_interactions_DEDUPLICATED.rds",
  "ECM_domains_transformed4ScType.xlsx"
)

missing_files <- required_files[!file.exists(file.path(stapp_data_dir, required_files))]
if (length(missing_files) > 0) {
  stop(
    "Missing beta reference files: ",
    paste(missing_files, collapse = ", ")
  )
}

standardize_genes <- function(genes) {
  genes <- as.character(genes)
  genes <- trimws(genes)
  genes <- genes[genes != "" & !is.na(genes)]

  if (length(genes) == 0) {
    return(character(0))
  }

  standardized <- tryCatch({
    if (!requireNamespace("scCustomize", quietly = TRUE)) {
      stop("scCustomize not installed")
    }
    result <- scCustomize::Updated_HGNC_Symbols(
      unique(genes),
      verbose = FALSE,
      case_check_as_warn = TRUE
    )
    mapping <- setNames(result$Output_Features, result$input_features)
    unname(mapping[genes])
  }, error = function(e) {
    checked <- HGNChelper::checkGeneSymbols(unique(genes), unmapped.as.na = FALSE)
    mapping <- setNames(checked$Suggested.Symbol, checked$x)
    out <- unname(mapping[genes])
    out[is.na(out) | out == ""] <- genes[is.na(out) | out == ""]
    out
  })

  standardized[is.na(standardized) | standardized == ""] <- genes[is.na(standardized) | standardized == ""]
  standardized
}

subcategory_map <- c(
  "basement_membrane" = "Basement Membrane",
  "hemostasis" = "Hemostasis",
  "elastic_fibers" = "Elastic fibers",
  "growth_factor-binding" = "Growth Factor-binding",
  "laminin_-_basement_membrane" = "Laminins",
  "matricellular" = "Matricellular proteins",
  "syndecan" = "Syndecan",
  "glypican" = "Glypican",
  "annexin" = "Annexin",
  "cathepsin" = "Cathepsin",
  "ccn_family" = "CCN Family",
  "cystatin" = "Cystatin",
  "facit" = "FACIT",
  "fibulin" = "Fibulin",
  "galectin" = "Galectin",
  "mucin" = "Mucin",
  "plexin" = "Plexin",
  "semaphorin" = "Semaphorin",
  "perivascular" = "Peri-vascular ECM"
)

parse_marker_genes <- function(value) {
  genes <- unlist(strsplit(ifelse(is.na(value), "", value), ","))
  genes <- trimws(genes)
  genes <- genes[genes != ""]
  unique(standardize_genes(genes))
}

# 1. Matrisome reference database
matrisome_db <- readRDS(file.path(stapp_data_dir, "matrisome_v2.rds"))
matrisome_db$gene <- standardize_genes(matrisome_db$gene)

# 2. Matrisome signatures derived from the reference table
matrisome_signatures <- list(
  ecm_glycoproteins = matrisome_db$gene[matrisome_db$notes == "ECM Glycoproteins"],
  collagens = matrisome_db$gene[matrisome_db$notes == "Collagens"],
  proteoglycans = matrisome_db$gene[matrisome_db$notes == "Proteoglycans"],
  ecm_regulators = matrisome_db$gene[matrisome_db$notes == "ECM Regulators"],
  secreted_factors = matrisome_db$gene[matrisome_db$notes == "Secreted Factors"],
  `ecm-affiliated_proteins` = matrisome_db$gene[matrisome_db$notes == "ECM-affiliated Proteins"]
)

for (sig_name in names(subcategory_map)) {
  matrisome_signatures[[sig_name]] <- matrisome_db$gene[
    grepl(subcategory_map[[sig_name]], matrisome_db$ecm_subcategory, fixed = TRUE)
  ]
}
matrisome_signatures <- lapply(matrisome_signatures, function(genes) unique(stats::na.omit(genes)))

# 3. ECM niche signatures
ecm_signatures <- readRDS(file.path(stapp_data_dir, "ecm_ucell_signatures.rds"))
ecm_signatures <- lapply(ecm_signatures, function(genes) unique(standardize_genes(genes)))

# 4. MatriComDB pairs: retain one row per unique heterotypic pair
matrisome_pairs <- readRDS(
  file.path(stapp_data_dir, "ultimate_ecm_interactions_DEDUPLICATED.rds")
)
names(matrisome_pairs)[seq_len(2)] <- c("Gene1", "Gene2")
matrisome_pairs$Gene1 <- standardize_genes(matrisome_pairs$Gene1)
matrisome_pairs$Gene2 <- standardize_genes(matrisome_pairs$Gene2)
matrisome_pairs <- matrisome_pairs[
  !is.na(matrisome_pairs$Gene1) & !is.na(matrisome_pairs$Gene2) &
    matrisome_pairs$Gene1 != "" & matrisome_pairs$Gene2 != "" &
    matrisome_pairs$Gene1 != matrisome_pairs$Gene2,
  ,
  drop = FALSE
]
pair_gene_1 <- pmin(matrisome_pairs$Gene1, matrisome_pairs$Gene2)
pair_gene_2 <- pmax(matrisome_pairs$Gene1, matrisome_pairs$Gene2)
matrisome_pairs$Gene1 <- pair_gene_1
matrisome_pairs$Gene2 <- pair_gene_2
matrisome_pairs <- matrisome_pairs[
  !duplicated(matrisome_pairs[, c("Gene1", "Gene2")]),
  ,
  drop = FALSE
]
rownames(matrisome_pairs) <- NULL

# 5. ECM markers and bundled workbook
xlsx_path <- file.path(stapp_data_dir, "ECM_domains_transformed4ScType.xlsx")
cell_markers <- openxlsx::read.xlsx(xlsx_path)
cell_markers <- cell_markers[cell_markers$tissueType == "ECM", ]

gs_positive <- lapply(cell_markers$geneSymbolmore1, parse_marker_genes)
names(gs_positive) <- cell_markers$cellName
gs_negative <- lapply(cell_markers$geneSymbolmore2, parse_marker_genes)
names(gs_negative) <- cell_markers$cellName
ecm_markers <- list(gs_positive = gs_positive, gs_negative = gs_negative)

# Keep the public package signature names stable even if the beta source
# changes how the niche signatures are stored.
normalized_ecm_signatures <- list(
  Interstitial = ecm_signatures[["Interstitial"]],
  Basement = ecm_signatures[["Basement"]],
  Vascular = gs_positive[["Vascular_ECM"]]
)
normalized_ecm_signatures <- normalized_ecm_signatures[!vapply(normalized_ecm_signatures, is.null, logical(1))]
ecm_signatures <- lapply(normalized_ecm_signatures, unique)

dir.create("inst/extdata", showWarnings = FALSE, recursive = TRUE)
file.copy(
  xlsx_path,
  file.path("inst", "extdata", "ECM_domains_transformed4ScType.xlsx"),
  overwrite = TRUE
)

usethis::use_data(matrisome_db, overwrite = TRUE)
usethis::use_data(matrisome_signatures, overwrite = TRUE)
usethis::use_data(ecm_signatures, overwrite = TRUE)
usethis::use_data(matrisome_pairs, overwrite = TRUE)
usethis::use_data(ecm_markers, overwrite = TRUE)

message("Data preparation complete.")
