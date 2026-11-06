#' @title ECM Differential Expression Analysis
#' @name differential-expression
#' @description Functions for differential expression analysis of ECM/matrisome genes
#' in spatial transcriptomics data.
NULL

# Declare global variables to avoid R CMD check NOTEs
utils::globalVariables(c("gene", "division", "notes", "p_val_adj", "cluster"))

#' Get ECM gene list
#'
#' @param categories Character vector of matrisome categories to filter.
#'   Can include division names ("Core matrisome", "Matrisome-associated") or
#'   subcategory names ("Collagens", "ECM Glycoproteins", etc.). If NULL, returns all.
#' @return Character vector of gene symbols
#' @keywords internal
#' @noRd
.get_ecm_genes <- function(categories = NULL) {
  db <- matrisome_db[!matrisome_db$gene %in% c("F8", "USH2A"), , drop = FALSE]

  if (is.null(categories)) {
    return(unique(db$gene))
  }

  # Match against both division and notes (subcategory) columns
  matched <- db$division %in% categories | db$notes %in% categories

  if (sum(matched) == 0) {
    warning("No genes found matching specified categories: ",
            paste(categories, collapse = ", "))
    return(character(0))
  }

  unique(db$gene[matched])
}

#' Annotate DE results with matrisome categories
#'
#' @param de_results Data frame with 'gene' column (or rownames as genes)
#' @return Data frame with added 'category' and 'subcategory' columns
#' @keywords internal
#' @noRd
.annotate_ecm_results <- function(de_results) {
  if (is.null(de_results) || nrow(de_results) == 0) {
    return(de_results)
  }

  db <- matrisome_db

  # Handle gene column or rownames
  if (!"gene" %in% colnames(de_results)) {
    de_results$gene <- rownames(de_results)
  }

  # Create lookup
  lookup <- db[, c("gene", "division", "notes")]
  colnames(lookup) <- c("gene", "category", "subcategory")
  lookup <- lookup[!duplicated(lookup$gene), ]

  # Merge
  result <- merge(de_results, lookup, by = "gene", all.x = TRUE, sort = FALSE)

  # Restore original order
  result <- result[match(de_results$gene, result$gene), ]
  rownames(result) <- result$gene

  result
}

#' Compute mean expression on the linear scale
#'
#' @param values Numeric vector
#' @param log_transformed Whether values are stored on a log1p scale
#' @return Numeric scalar
#' @keywords internal
#' @noRd
.mean_on_linear_scale <- function(values, log_transformed = TRUE) {
  values <- as.numeric(values)
  if (log_transformed) {
    return(mean(expm1(values)))
  }
  mean(values)
}

#' Find differentially expressed ECM genes between groups
#'
#' Identifies differentially expressed matrisome genes between two groups
#' (clusters, conditions, etc.) using Seurat's FindMarkers framework.
#'
#' @param seurat_obj A Seurat object
#' @param ident.1 Identity of the first group (cluster name, cell names, or condition)
#' @param ident.2 Identity of the second group. If NULL, compares ident.1 vs all other cells.
#' @param group.by Metadata column to use for grouping. If NULL, uses active identities.
#' @param assay Assay to use. If NULL, auto-detects (SCT > Spatial > first available).
#' @param test.use Statistical test: "wilcox" (default), "MAST", "DESeq2", "LR", "t"
#' @param min.pct Minimum fraction of cells expressing the gene in either group (default: 0.05)
#' @param logfc.threshold Minimum log2 fold change threshold (default: 0.25)
#' @param categories Matrisome categories to include. Can be divisions ("Core matrisome")
#'   or subcategories ("Collagens"). If NULL, uses all ECM genes.
#' @param only.pos Only return genes with positive log2FC (upregulated in ident.1)
#' @param p.adjust.method Method for p-value adjustment (default: "BH")
#' @param verbose Print progress messages
#'
#' @return Data frame with columns:
#' \describe{
#'   \item{gene}{Gene symbol}
#'   \item{p_val}{Raw p-value}
#'   \item{avg_log2FC}{Average log2 fold change}
#'   \item{pct.1}{Fraction of cells in ident.1 expressing the gene}
#'   \item{pct.2}{Fraction of cells in ident.2 expressing the gene}
#'   \item{p_val_adj}{Adjusted p-value}
#'   \item{category}{Matrisome division (Core/Associated)}
#'   \item{subcategory}{Matrisome subcategory}
#' }
#'
#' @examples
#' \dontrun{
#' # Compare fibroblasts vs all other cells
#' markers <- find_ecm_markers(seurat_obj, ident.1 = "Fibroblasts",
#'                             group.by = "cell_type")
#'
#' # Compare two specific clusters
#' markers <- find_ecm_markers(seurat_obj, ident.1 = "1", ident.2 = "2",
#'                             group.by = "seurat_clusters")
#'
#' # Filter to collagens only
#' collagen_de <- find_ecm_markers(seurat_obj, ident.1 = "Stroma",
#'                                  categories = "Collagens")
#' }
#'
#' @export
find_ecm_markers <- function(seurat_obj,
                             ident.1,
                             ident.2 = NULL,
                             group.by = NULL,
                             assay = NULL,
                             test.use = "wilcox",
                             min.pct = 0.05,
                             logfc.threshold = 0.25,
                             categories = NULL,
                             only.pos = FALSE,
                             p.adjust.method = "BH",
                             verbose = TRUE) {
  .assert_seurat_object(seurat_obj, "seurat_obj")

  if (!is.null(group.by) && !group.by %in% colnames(seurat_obj@meta.data)) {
    stop("group.by '", group.by, "' not found in metadata")
  }

  # Determine assay
  assay <- .resolve_assay(seurat_obj, assay)

  .message_if(verbose, "Using assay: ", assay)

  # Get ECM genes
  ecm_genes <- .get_ecm_genes(categories)
  if (length(ecm_genes) == 0) {
    stop("No ECM genes found for specified categories")
  }

  # Filter to genes present in the object
  # Try data layer first, fallback to counts for Seurat v5
  genes_in_data <- .safe_get_rownames(seurat_obj, assay = assay)
  ecm_genes_present <- intersect(ecm_genes, genes_in_data)

  if (length(ecm_genes_present) == 0) {
    stop("No ECM genes found in the expression data")
  }

  .message_if(
    verbose,
    sprintf(
      "Testing %d ECM genes (of %d total matrisome genes)",
      length(ecm_genes_present),
      length(ecm_genes)
    )
  )

  # Set identities if group.by specified
  if (!is.null(group.by)) {
    Seurat::Idents(seurat_obj) <- seurat_obj[[group.by, drop = TRUE]]
  }

  # Validate ident.1
  all_idents <- levels(Seurat::Idents(seurat_obj))
  all_cells <- colnames(seurat_obj)
  is_valid_selection <- function(ident) {
    if (is.null(ident) || length(ident) == 0) {
      return(FALSE)
    }
    all(ident %in% all_idents) || all(ident %in% all_cells)
  }

  if (!is_valid_selection(ident.1)) {
    stop(
      "ident.1 was not found in identities or cell names. Available identities: ",
      paste(all_idents, collapse = ", ")
    )
  }
  if (!is.null(ident.2) && !is_valid_selection(ident.2)) {
    stop(
      "ident.2 was not found in identities or cell names. Available identities: ",
      paste(all_idents, collapse = ", ")
    )
  }

  # Prepare SCT if needed
  if (assay == "SCT") {
    tryCatch({
      seurat_obj <- Seurat::PrepSCTFindMarkers(seurat_obj, verbose = FALSE)
    }, error = function(e) {
      .message_if(verbose, "PrepSCTFindMarkers not needed or failed: ", e$message)
    })
  }

  # Run FindMarkers
  .message_if(verbose, "Running differential expression test: ", test.use)

  de_results <- tryCatch({
    Seurat::FindMarkers(
      seurat_obj,
      ident.1 = ident.1,
      ident.2 = ident.2,
      features = ecm_genes_present,
      test.use = test.use,
      min.pct = min.pct,
      logfc.threshold = logfc.threshold,
      only.pos = only.pos,
      verbose = verbose
    )
  }, error = function(e) {
    stop("FindMarkers failed: ", e$message)
  })

  if (is.null(de_results) || nrow(de_results) == 0) {
    .message_if(verbose, "No significant markers found")
    return(data.frame(
      gene = character(),
      p_val = numeric(),
      avg_log2FC = numeric(),
      pct.1 = numeric(),
      pct.2 = numeric(),
      p_val_adj = numeric(),
      category = character(),
      subcategory = character()
    ))
  }

  # Add gene column from rownames
  de_results$gene <- rownames(de_results)

  # Recalculate adjusted p-values for ECM gene subset
  de_results$p_val_adj <- stats::p.adjust(de_results$p_val, method = p.adjust.method)

  # Annotate with matrisome categories
  de_results <- .annotate_ecm_results(de_results)

  # Reorder columns
  col_order <- c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2",
                 "p_val_adj", "category", "subcategory")
  de_results <- de_results[, col_order]

  # Sort by adjusted p-value
  de_results <- de_results[order(de_results$p_val_adj), ]
  rownames(de_results) <- de_results$gene

  .message_if(verbose, "Found ", nrow(de_results), " differentially expressed ECM genes")

  de_results
}

#' Find spatially corrected ECM markers between two groups
#'
#' Runs a two-stage workflow designed for within-sample comparisons where
#' spatial autocorrelation matters. ECM genes are first screened with a
#' Wilcoxon rank-sum test, then significant hits are re-tested with a
#' \pkg{spaMM} spatial mixed model using Matern covariance.
#'
#' @param seurat_obj A Seurat object
#' @param group.by Metadata column defining the comparison groups
#' @param ident.1 Group of interest
#' @param ident.2 Reference group
#' @param assay Assay to use. If NULL, auto-detects.
#' @param layer Expression layer to test. Defaults to `"data"`.
#' @param min.pct Minimum fraction of spots expressing the gene in either group.
#'   Defaults to `0.05`.
#' @param logfc.threshold Minimum absolute log2 fold-change required to send a
#'   gene to the spatial model. Defaults to `0.25`.
#' @param wilcox.padj.threshold Adjusted p-value threshold for Wilcoxon
#'   pre-screening. Defaults to `0.05`.
#' @param categories Optional matrisome categories to restrict testing to.
#' @param sp_topfrac Optional fraction of Wilcoxon-significant genes to pass to
#'   \pkg{spaMM}, ranked by adjusted p-value. Defaults to `1`.
#' @param verbose Print progress messages
#'
#' @return Data frame containing Wilcoxon statistics, spaMM statistics, and
#'   matrisome annotations.
#'
#' @details
#' This function mirrors the differential expression strategy used in the
#' MatriSpace case study: Wilcoxon pre-screening on ECM genes followed by
#' spatial correction with `spaMM::fitme(gene_expr ~ group + Matern(1 | x + y))`.
#'
#' @examples
#' \dontrun{
#' spatial_de <- find_spatially_corrected_ecm_markers(
#'   seurat_obj,
#'   group.by = "CompositionCluster_CC",
#'   ident.1 = "CC1",
#'   ident.2 = "CC3"
#' )
#' }
#'
#' @export
find_spatially_corrected_ecm_markers <- function(seurat_obj,
                                                 group.by,
                                                 ident.1,
                                                 ident.2,
                                                 assay = NULL,
                                                 layer = "data",
                                                 min.pct = 0.05,
                                                 logfc.threshold = 0.25,
                                                 wilcox.padj.threshold = 0.05,
                                                 categories = NULL,
                                                 sp_topfrac = 1,
                                                 verbose = TRUE) {
  .assert_seurat_object(seurat_obj, "seurat_obj")

  if (!requireNamespace("spaMM", quietly = TRUE)) {
    stop("Package 'spaMM' is required for spatially corrected differential expression")
  }
  if (!group.by %in% colnames(seurat_obj@meta.data)) {
    stop("group.by '", group.by, "' not found in metadata")
  }

  assay <- .resolve_assay(seurat_obj, assay)
  expr <- .safe_get_assay_data(seurat_obj, assay = assay, slot = layer)
  if (is.null(expr)) {
    stop("Could not retrieve assay data from layer '", layer, "'")
  }

  groups <- as.character(seurat_obj[[group.by, drop = TRUE]])
  idx_1 <- which(groups == ident.1)
  idx_2 <- which(groups == ident.2)
  idx_all <- c(idx_1, idx_2)

  if (length(idx_1) < 3 || length(idx_2) < 3) {
    stop("Need at least 3 spots per group for spatially corrected DE")
  }

  ecm_genes <- intersect(.get_ecm_genes(categories), rownames(expr))
  if (length(ecm_genes) == 0) {
    stop("No ECM genes found in the selected assay layer")
  }

  pct_1 <- Matrix::rowMeans(expr[ecm_genes, idx_1, drop = FALSE] > 0)
  pct_2 <- Matrix::rowMeans(expr[ecm_genes, idx_2, drop = FALSE] > 0)
  keep_genes <- ecm_genes[pmax(pct_1, pct_2) >= min.pct]

  if (length(keep_genes) == 0) {
    stop("No ECM genes passed the min.pct filter")
  }

  .message_if(verbose, "Running Wilcoxon pre-screen on ", length(keep_genes), " ECM genes")

  log_transformed <- layer != "counts"
  wilcox_results <- lapply(keep_genes, function(gene_name) {
    vals_1 <- as.numeric(expr[gene_name, idx_1])
    vals_2 <- as.numeric(expr[gene_name, idx_2])
    p_value <- tryCatch(
      stats::wilcox.test(vals_1, vals_2)$p.value,
      error = function(e) 1
    )

    mean_1 <- .mean_on_linear_scale(vals_1, log_transformed = log_transformed)
    mean_2 <- .mean_on_linear_scale(vals_2, log_transformed = log_transformed)

    data.frame(
      gene = gene_name,
      avg_expr_1 = mean(vals_1),
      avg_expr_2 = mean(vals_2),
      mean_linear_1 = mean_1,
      mean_linear_2 = mean_2,
      avg_log2FC = log2(mean_1 + 1) - log2(mean_2 + 1),
      pct.1 = pct_1[gene_name],
      pct.2 = pct_2[gene_name],
      wilcox_pval = p_value,
      stringsAsFactors = FALSE
    )
  })

  result <- dplyr::bind_rows(wilcox_results)
  result$wilcox_padj <- stats::p.adjust(result$wilcox_pval, method = "BH")
  result$wilcox_sig <- result$wilcox_padj < wilcox.padj.threshold &
    abs(result$avg_log2FC) > logfc.threshold

  sig_genes <- result$gene[result$wilcox_sig]
  if (sp_topfrac < 1 && length(sig_genes) > 0) {
    n_keep <- max(1, ceiling(length(sig_genes) * sp_topfrac))
    sig_genes <- result[result$wilcox_sig, ]
    sig_genes <- sig_genes[order(sig_genes$wilcox_padj), "gene"][seq_len(n_keep)]
  }

  coords <- .safe_get_tissue_coordinates(seurat_obj)[idx_all, , drop = FALSE]
  comparison_group <- factor(
    ifelse(groups[idx_all] == ident.1, ident.1, ident.2),
    levels = c(ident.2, ident.1)
  )
  spot_df <- data.frame(
    comparison_group = comparison_group,
    xpos = coords$imagecol,
    ypos = coords$imagerow
  )

  .message_if(verbose, "Running spaMM spatial correction on ", length(sig_genes), " genes")

  spamm_results <- lapply(seq_along(sig_genes), function(i) {
    gene_name <- sig_genes[i]
    spot_df$gene_expr <- as.numeric(expr[gene_name, idx_all])

    fit <- tryCatch(
      suppressMessages(
        spaMM::fitme(
          gene_expr ~ comparison_group + Matern(1 | xpos + ypos),
          data = spot_df,
          fixed = list(nu = 0.5),
          method = "REML",
          control.HLfit = list(algebra = "decorr"),
          verbose = FALSE
        )
      ),
      error = function(e) NULL
    )

    if (isTRUE(verbose) && (i %% 20 == 0 || i == length(sig_genes))) {
      message("  Processed ", i, " / ", length(sig_genes), " spaMM fits")
    }

    if (is.null(fit)) {
      return(data.frame(
        gene = gene_name,
        spamm_coef = NA_real_,
        spamm_se = NA_real_,
        spamm_pval = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    beta_tab <- summary.HLfit(fit, details = c(p_value = "Wald"), verbose = FALSE)$beta_table
    coef_row <- grep("^comparison_group", rownames(beta_tab))[1]
    if (is.na(coef_row)) {
      return(data.frame(
        gene = gene_name,
        spamm_coef = NA_real_,
        spamm_se = NA_real_,
        spamm_pval = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    data.frame(
      gene = gene_name,
      spamm_coef = beta_tab[coef_row, 1],
      spamm_se = beta_tab[coef_row, 2],
      spamm_pval = beta_tab[coef_row, 4],
      stringsAsFactors = FALSE
    )
  })

  if (length(spamm_results) > 0) {
    spamm_df <- dplyr::bind_rows(spamm_results)
    spamm_df$spamm_padj <- stats::p.adjust(spamm_df$spamm_pval, method = "BH")
  } else {
    spamm_df <- data.frame(
      gene = character(),
      spamm_coef = numeric(),
      spamm_se = numeric(),
      spamm_pval = numeric(),
      spamm_padj = numeric()
    )
  }

  result <- result |>
    dplyr::left_join(spamm_df, by = "gene") |>
    dplyr::mutate(
      spamm_sig = !is.na(.data$spamm_padj) & .data$spamm_padj < 0.05,
      survived_spatial = .data$wilcox_sig & .data$spamm_sig
    )

  result <- .annotate_ecm_results(result)
  result <- result[order(result$wilcox_padj, result$spamm_padj), ]
  rownames(result) <- result$gene
  result
}

#' Find ECM markers for all clusters
#'
#' Identifies differentially expressed ECM genes for each cluster compared to
#' all other cells. Wrapper around \code{find_ecm_markers} for all clusters.
#'
#' @inheritParams find_ecm_markers
#' @param min.diff.pct Minimum difference in expression percentage between groups
#'
#' @return Data frame with columns from \code{find_ecm_markers} plus:
#' \describe{
#'   \item{cluster}{Cluster/group identifier}
#' }
#'
#' @examples
#' \dontrun{
#' # Find ECM markers for all clusters
#' all_markers <- find_all_ecm_markers(seurat_obj, group.by = "seurat_clusters")
#'
#' # Filter to core matrisome only
#' core_markers <- find_all_ecm_markers(seurat_obj,
#'                                       categories = "Core matrisome")
#' }
#'
#' @export
find_all_ecm_markers <- function(seurat_obj,
                                  group.by = NULL,
                                  assay = NULL,
                                  test.use = "wilcox",
                                  min.pct = 0.05,
                                  logfc.threshold = 0.25,
                                  min.diff.pct = -Inf,
                                  categories = NULL,
                                  only.pos = TRUE,
                                  p.adjust.method = "BH",
                                  verbose = TRUE) {
  .assert_seurat_object(seurat_obj, "seurat_obj")

  if (!is.null(group.by) && !group.by %in% colnames(seurat_obj@meta.data)) {
    stop("group.by '", group.by, "' not found in metadata")
  }

  # Set identities
  if (!is.null(group.by)) {
    Seurat::Idents(seurat_obj) <- seurat_obj[[group.by, drop = TRUE]]
  }

  all_clusters <- levels(Seurat::Idents(seurat_obj))

  if (length(all_clusters) < 2) {
    stop("Need at least 2 clusters for comparison. Found: ", length(all_clusters))
  }

  .message_if(verbose, "Finding ECM markers for ", length(all_clusters), " clusters")

  # Run for each cluster
  all_results <- lapply(seq_along(all_clusters), function(i) {
    cluster_name <- all_clusters[i]

    .message_if(
      verbose,
      sprintf("  [%d/%d] Processing cluster: %s", i, length(all_clusters), cluster_name)
    )

    result <- tryCatch({
      find_ecm_markers(
        seurat_obj,
        ident.1 = cluster_name,
        ident.2 = NULL,
        group.by = group.by,
        assay = assay,
        test.use = test.use,
        min.pct = min.pct,
        logfc.threshold = logfc.threshold,
        categories = categories,
        only.pos = only.pos,
        p.adjust.method = "none",  # Will adjust globally
        verbose = FALSE
      )
    }, error = function(e) {
      if (isTRUE(verbose)) warning("Failed for cluster ", cluster_name, ": ", e$message)
      return(NULL)
    })

    if (!is.null(result) && nrow(result) > 0) {
      result$cluster <- cluster_name

      # Apply min.diff.pct filter
      if (is.finite(min.diff.pct)) {
        result <- result[result$pct.1 - result$pct.2 >= min.diff.pct, ]
      }
    }

    result
  })

  # Combine results
  combined <- dplyr::bind_rows(all_results)

  if (nrow(combined) == 0) {
    .message_if(verbose, "No significant markers found in any cluster")
    return(data.frame(
      gene = character(),
      p_val = numeric(),
      avg_log2FC = numeric(),
      pct.1 = numeric(),
      pct.2 = numeric(),
      p_val_adj = numeric(),
      category = character(),
      subcategory = character(),
      cluster = character()
    ))
  }

  # Global FDR correction
  combined$p_val_adj <- stats::p.adjust(combined$p_val, method = p.adjust.method)

  # Reorder columns
  col_order <- c("cluster", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2",
                 "p_val_adj", "category", "subcategory")
  combined <- combined[, col_order]

  # Sort by cluster then p-value
  combined <- combined[order(combined$cluster, combined$p_val_adj), ]
  rownames(combined) <- NULL

  .message_if(verbose, "Found ", nrow(combined), " total ECM markers across all clusters")

  combined
}

#' Differential expression of ECM genes between samples or conditions
#'
#' Compares ECM gene expression between different samples or experimental conditions.
#' Supports both spot-level and pseudo-bulk analysis approaches.
#'
#' @param seurat_obj A Seurat object
#' @param condition.col Metadata column containing condition/sample labels
#' @param condition.1 First condition to compare
#' @param condition.2 Second condition to compare
#' @param cluster.col Optional metadata column for cluster-aware analysis.
#'   In spot-level mode, DE is run separately within each cluster. In
#'   pseudobulk mode, counts are aggregated at the sample-by-cluster level and
#'   cluster is added to the DESeq2 design.
#' @param assay Assay to use. If NULL, auto-detects.
#' @param test.use Statistical test: "wilcox" (default), "MAST", "DESeq2", "LR"
#' @param min.pct Minimum fraction expressing (default: 0.05)
#' @param logfc.threshold Minimum log2 fold change (default: 0.25)
#' @param categories Matrisome categories to include. If NULL, uses all ECM genes.
#' @param pseudobulk If TRUE, aggregates to pseudo-bulk per sample before DE.
#'   Recommended for multi-sample experiments.
#' @param sample.col Required if pseudobulk=TRUE. Column identifying individual samples.
#' @param blocking.col Optional metadata column added to the DESeq2 design in
#'   pseudobulk mode to control for paired subjects, batches, or other
#'   blocking effects.
#' @param min.cells.per.sample Minimum cells per sample for pseudo-bulk (default: 10)
#' @param verbose Print progress messages
#'
#' @return Data frame with DE results. If cluster.col is provided, includes
#'   a 'cluster' column for cluster-stratified results.
#'
#' @details
#' ## Analysis Modes
#'
#' **Spot-level (pseudobulk=FALSE):**
#' Treats each spot as independent observation. Fast but may inflate
#' significance due to pseudo-replication within samples.
#'
#' **Pseudo-bulk (pseudobulk=TRUE):**
#' Aggregates counts by sample, then compares between conditions with DESeq2.
#' When `cluster.col` is supplied, profiles are aggregated at the
#' sample-by-cluster level and cluster is included in the design. Use
#' `blocking.col` to adjust for repeated subjects or batch effects.
#'
#' ## Cluster-Stratified Analysis
#'
#' When cluster.col is provided, DE is run separately within each cluster.
#' This reveals cluster-specific expression changes between conditions.
#'
#' @examples
#' \dontrun{
#' # Basic comparison between conditions
#' de_results <- find_ecm_de_samples(seurat_obj,
#'                                    condition.col = "condition",
#'                                    condition.1 = "tumor",
#'                                    condition.2 = "normal")
#'
#' # Cluster-stratified analysis
#' stratified_de <- find_ecm_de_samples(seurat_obj,
#'                                       condition.col = "condition",
#'                                       condition.1 = "tumor",
#'                                       condition.2 = "normal",
#'                                       cluster.col = "cell_type")
#'
#' # Pseudo-bulk analysis (recommended for multiple samples)
#' pseudobulk_de <- find_ecm_de_samples(seurat_obj,
#'                                       condition.col = "condition",
#'                                       condition.1 = "treated",
#'                                       condition.2 = "control",
#'                                       pseudobulk = TRUE,
#'                                       sample.col = "sample_id",
#'                                       blocking.col = "patient_id")
#' }
#'
#' @export
find_ecm_de_samples <- function(seurat_obj,
                                 condition.col,
                                 condition.1,
                                 condition.2,
                                 cluster.col = NULL,
                                 assay = NULL,
                                 test.use = "wilcox",
                                 min.pct = 0.05,
                                 logfc.threshold = 0.25,
                                 categories = NULL,
                                 pseudobulk = FALSE,
                                 sample.col = NULL,
                                 blocking.col = NULL,
                                 min.cells.per.sample = 10,
                                 verbose = TRUE) {
  .assert_seurat_object(seurat_obj, "seurat_obj")

  if (!condition.col %in% colnames(seurat_obj@meta.data)) {
    stop("condition.col '", condition.col, "' not found in metadata")
  }

  conditions <- seurat_obj[[condition.col, drop = TRUE]]
  unique_conditions <- unique(conditions)

  if (!condition.1 %in% unique_conditions) {
    stop("condition.1 '", condition.1, "' not found. Available: ",
         paste(unique_conditions, collapse = ", "))
  }
  if (!condition.2 %in% unique_conditions) {
    stop("condition.2 '", condition.2, "' not found. Available: ",
         paste(unique_conditions, collapse = ", "))
  }

  if (pseudobulk && is.null(sample.col)) {
    stop("sample.col is required when pseudobulk = TRUE")
  }

  if (!is.null(sample.col) && !sample.col %in% colnames(seurat_obj@meta.data)) {
    stop("sample.col '", sample.col, "' not found in metadata")
  }
  if (!is.null(blocking.col) && !blocking.col %in% colnames(seurat_obj@meta.data)) {
    stop("blocking.col '", blocking.col, "' not found in metadata")
  }

  # Determine assay
  assay <- .resolve_assay(seurat_obj, assay)

  # Get ECM genes
  ecm_genes <- .get_ecm_genes(categories)
  genes_in_data <- .safe_get_rownames(seurat_obj, assay = assay)
  ecm_genes_present <- intersect(ecm_genes, genes_in_data)

  if (length(ecm_genes_present) == 0) {
    stop("No ECM genes found in the expression data")
  }

  .message_if(verbose, sprintf("Comparing %s vs %s", condition.1, condition.2))
  .message_if(verbose, sprintf("Testing %d ECM genes", length(ecm_genes_present)))

  # Handle pseudo-bulk analysis
  if (pseudobulk) {
    return(.run_pseudobulk_de(
      seurat_obj = seurat_obj,
      condition.col = condition.col,
      condition.1 = condition.1,
      condition.2 = condition.2,
      sample.col = sample.col,
      cluster.col = cluster.col,
      blocking.col = blocking.col,
      assay = assay,
      ecm_genes = ecm_genes_present,
      min.cells.per.sample = min.cells.per.sample,
      verbose = verbose
    ))
  }

  # Spot-level analysis
  if (is.null(cluster.col)) {
    # Simple two-condition comparison
    Seurat::Idents(seurat_obj) <- seurat_obj[[condition.col, drop = TRUE]]

    de_results <- find_ecm_markers(
      seurat_obj,
      ident.1 = condition.1,
      ident.2 = condition.2,
      group.by = NULL,
      assay = assay,
      test.use = test.use,
      min.pct = min.pct,
      logfc.threshold = logfc.threshold,
      categories = categories,
      only.pos = FALSE,
      verbose = verbose
    )

    return(de_results)

  } else {
    # Cluster-stratified analysis
    if (!cluster.col %in% colnames(seurat_obj@meta.data)) {
      stop("cluster.col '", cluster.col, "' not found in metadata")
    }

    clusters <- unique(seurat_obj[[cluster.col, drop = TRUE]])

    .message_if(verbose, "Running cluster-stratified analysis for ", length(clusters), " clusters")

    all_results <- lapply(clusters, function(cl) {
      .message_if(verbose, "  Processing cluster: ", cl)

      # Subset to cluster
      cells_in_cluster <- colnames(seurat_obj)[seurat_obj[[cluster.col, drop = TRUE]] == cl]

      if (length(cells_in_cluster) < 10) {
        .message_if(verbose, "    Skipping - too few cells (", length(cells_in_cluster), ")")
        return(NULL)
      }

      cluster_obj <- subset(seurat_obj, cells = cells_in_cluster)

      # Check both conditions have cells
      cond_table <- table(cluster_obj[[condition.col, drop = TRUE]])
      if (!condition.1 %in% names(cond_table) || !condition.2 %in% names(cond_table)) {
        .message_if(verbose, "    Skipping - missing condition in cluster")
        return(NULL)
      }
      if (cond_table[condition.1] < 3 || cond_table[condition.2] < 3) {
        .message_if(verbose, "    Skipping - too few cells in one condition")
        return(NULL)
      }

      # Run DE
      Seurat::Idents(cluster_obj) <- cluster_obj[[condition.col, drop = TRUE]]

      result <- tryCatch({
        find_ecm_markers(
          cluster_obj,
          ident.1 = condition.1,
          ident.2 = condition.2,
          assay = assay,
          test.use = test.use,
          min.pct = min.pct,
          logfc.threshold = logfc.threshold,
          categories = categories,
          only.pos = FALSE,
          p.adjust.method = "none",
          verbose = FALSE
        )
      }, error = function(e) {
        if (isTRUE(verbose)) warning("    Failed: ", e$message)
        return(NULL)
      })

      if (!is.null(result) && nrow(result) > 0) {
        result$cluster <- cl
      }

      result
    })

    # Combine
    combined <- dplyr::bind_rows(all_results)

    if (nrow(combined) == 0) {
      .message_if(verbose, "No significant DE genes found")
      return(data.frame(
        gene = character(),
        cluster = character(),
        p_val = numeric(),
        avg_log2FC = numeric(),
        pct.1 = numeric(),
        pct.2 = numeric(),
        p_val_adj = numeric(),
        category = character(),
        subcategory = character()
      ))
    }

    # Global FDR correction
    combined$p_val_adj <- stats::p.adjust(combined$p_val, method = "BH")

    # Sort
    combined <- combined[order(combined$cluster, combined$p_val_adj), ]
    rownames(combined) <- NULL

    .message_if(verbose, "Found ", nrow(combined), " DE genes across clusters")

    return(combined)
  }
}

#' Run pseudo-bulk differential expression
#'
#' @keywords internal
#' @noRd
.run_pseudobulk_de <- function(seurat_obj,
                               condition.col,
                               condition.1,
                               condition.2,
                               sample.col,
                               cluster.col,
                               blocking.col,
                               assay,
                               ecm_genes,
                               min.cells.per.sample,
                               verbose) {
  .message_if(verbose, "Running pseudo-bulk analysis...")

  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required for pseudobulk differential expression")
  }

  # Get expression matrix
  expr_matrix <- .safe_get_assay_data(seurat_obj, assay = assay, slot = "counts")
  if (is.null(expr_matrix)) {
    stop("Could not retrieve counts data for pseudobulk analysis")
  }

  # Filter to ECM genes
  ecm_genes_present <- intersect(ecm_genes, rownames(expr_matrix))
  expr_matrix <- expr_matrix[ecm_genes_present, , drop = FALSE]

  # Get metadata
  meta <- seurat_obj@meta.data
  meta$cell_id <- rownames(meta)

  # Define grouping
  if (is.null(cluster.col)) {
    meta$group <- meta[[sample.col]]
  } else {
    meta$group <- paste(meta[[sample.col]], meta[[cluster.col]], sep = "_")
  }

  # Filter to conditions of interest
  cells_keep <- meta[[condition.col]] %in% c(condition.1, condition.2)
  meta <- meta[cells_keep, ]
  expr_matrix <- expr_matrix[, meta$cell_id, drop = FALSE]

  # Aggregate by group
  groups <- unique(meta$group)

  .message_if(verbose, "  Aggregating ", length(groups), " sample-groups...")

  pseudobulk_list <- lapply(groups, function(g) {
    cells <- meta$cell_id[meta$group == g]
    if (length(cells) < min.cells.per.sample) return(NULL)

    counts <- Matrix::rowSums(expr_matrix[, cells, drop = FALSE])
    condition <- unique(meta[[condition.col]][meta$group == g])
    sample_id <- unique(meta[[sample.col]][meta$group == g])

    if (!is.null(cluster.col)) {
      cluster_id <- unique(meta[[cluster.col]][meta$group == g])
    } else {
      cluster_id <- "all"
    }

    list(
      counts = counts,
      condition = condition[1],
      sample = sample_id[1],
      cluster = cluster_id[1],
      blocking = if (!is.null(blocking.col)) unique(meta[[blocking.col]][meta$group == g])[1] else "all",
      n_cells = length(cells)
    )
  })

  # Remove NULL entries
  pseudobulk_list <- pseudobulk_list[!sapply(pseudobulk_list, is.null)]

  if (length(pseudobulk_list) < 4) {
    warning("Too few samples for pseudo-bulk analysis (need >= 2 per condition)")
    return(data.frame(
      gene = character(),
      p_val = numeric(),
      avg_log2FC = numeric(),
      p_val_adj = numeric(),
      category = character(),
      subcategory = character()
    ))
  }

  # Build count matrix
  count_matrix <- do.call(cbind, lapply(pseudobulk_list, function(x) x$counts))
  colnames(count_matrix) <- sapply(pseudobulk_list, function(x) {
    paste(x$sample, x$cluster, sep = "_")
  })

  sample_info <- data.frame(
    sample = sapply(pseudobulk_list, function(x) x$sample),
    condition = factor(sapply(pseudobulk_list, function(x) x$condition)),
    cluster = sapply(pseudobulk_list, function(x) x$cluster),
    blocking = sapply(pseudobulk_list, function(x) x$blocking),
    n_cells = sapply(pseudobulk_list, function(x) x$n_cells)
  )

  sample_info$condition <- stats::relevel(sample_info$condition, ref = condition.2)
  sample_info$cluster <- factor(sample_info$cluster)
  sample_info$blocking <- factor(sample_info$blocking)

  condition_counts <- table(as.character(sample_info$condition))
  counts_to_check <- vapply(
    c(condition.1, condition.2),
    function(label) if (label %in% names(condition_counts)) condition_counts[[label]] else 0L,
    integer(1)
  )
  if (any(counts_to_check < 2)) {
    warning("Need at least 2 pseudobulk profiles per condition for DESeq2 analysis")
    return(data.frame(
      gene = character(),
      p_val = numeric(),
      avg_log2FC = numeric(),
      p_val_adj = numeric(),
      category = character(),
      subcategory = character()
    ))
  }

  design_terms <- c()
  if (!is.null(cluster.col)) {
    design_terms <- c(design_terms, "cluster")
  }
  if (!is.null(blocking.col)) {
    design_terms <- c(design_terms, "blocking")
  }
  design_terms <- c(design_terms, "condition")
  design_formula <- stats::as.formula(paste("~", paste(design_terms, collapse = " + ")))

  .message_if(verbose, "  DESeq2 design: ", deparse(design_formula))

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(count_matrix),
    colData = sample_info,
    design = design_formula
  )
  dds <- DESeq2::DESeq(dds, quiet = !isTRUE(verbose))

  res <- DESeq2::results(dds, contrast = c("condition", condition.1, condition.2))
  de_results <- as.data.frame(res)
  de_results$gene <- rownames(de_results)
  rownames(de_results) <- NULL
  de_results$p_val <- de_results$pvalue
  de_results$p_val_adj <- de_results$padj
  de_results$avg_log2FC <- de_results$log2FoldChange

  de_results <- .annotate_ecm_results(de_results)
  de_results <- de_results[order(de_results$p_val_adj, de_results$p_val), ]
  rownames(de_results) <- de_results$gene

  .message_if(
    verbose,
    "  Found ",
    sum(de_results$p_val_adj < 0.05, na.rm = TRUE),
    " significant genes (FDR < 0.05)"
  )

  de_results
}
