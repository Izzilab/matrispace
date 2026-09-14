#' @title Matrisome Pair Co-expression Functions
#' @name matrisome-pairs
#' @description Functions for spatial co-expression analysis of MatriComDB
#'   pairs in which at least one member is a matrisome gene.
NULL

.prepare_matrisome_pairs <- function(pair_db) {
  required_columns <- c("Ligand", "Receptor")
  missing_columns <- setdiff(required_columns, colnames(pair_db))
  if (length(missing_columns) > 0) {
    stop(
      "The pair database is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      "."
    )
  }

  gene_1 <- trimws(as.character(pair_db$Ligand))
  gene_2 <- trimws(as.character(pair_db$Receptor))
  keep <- !is.na(gene_1) & !is.na(gene_2) &
    nzchar(gene_1) & nzchar(gene_2) & gene_1 != gene_2
  pair_db <- pair_db[keep, , drop = FALSE]
  gene_1 <- gene_1[keep]
  gene_2 <- gene_2[keep]

  pair_db$Ligand <- pmin(gene_1, gene_2)
  pair_db$Receptor <- pmax(gene_1, gene_2)
  pair_db <- pair_db[!duplicated(pair_db[, required_columns, drop = FALSE]), , drop = FALSE]
  rownames(pair_db) <- NULL
  pair_db
}

#' Calculate spatial matrisome-pair co-expression scores
#'
#' Computes spatially weighted co-expression scores for unique heterotypic
#' MatriComDB pairs at each spot. The method measures whether one member of a
#' pair is expressed in a spot while its partner is expressed in neighboring
#' spots; it does not infer signaling direction or communication activity.
#'
#' @param seurat_obj A Seurat object
#' @param lr_db MatriComDB pair data with legacy `Ligand` and `Receptor`
#'   columns identifying the two pair members. Defaults to package data.
#' @param adj_matrix Sparse adjacency matrix for spatial neighbors.
#'   If NULL, computed from coordinates.
#' @param assay Assay to use (default: "SCT")
#' @param layer Layer within assay (default: "data")
#' @param verbose Print progress messages
#'
#' @return Seurat object with a "LRACTIVITY" assay added. The assay and
#'   argument names are retained for backward compatibility.
#'
#' @details
#' ## Scoring Formula
#'
#' \deqn{Score_{spot} = \sqrt{Gene2_{spot}} \times \sum_{neighbors} \sqrt{Gene1_{neighbor}}}
#'
#' Pair members are alphabetically ordered before scoring so reciprocal rows
#' collapse to one stable feature name. Homomeric rows are excluded.
#'
#' **Product structure:** A score is non-zero only when both pair members are
#' detected in the relevant spot and neighborhood.
#'
#' **Square Root Transformation:**
#' \itemize{
#'   \item Variance stabilization for skewed expression data
#'   \item Reduced influence of highly expressed genes
#'   \item Improved robustness to outliers
#' }
#'
#' @examples
#' \dontrun{
#' seurat_obj <- score_lr_activity(seurat_obj)
#' }
#'
#' @export
score_lr_activity <- function(seurat_obj,
                              lr_db = NULL,
                              adj_matrix = NULL,
                              assay = NULL,
                              layer = "data",
                              verbose = TRUE) {
  .assert_seurat_object(seurat_obj, "seurat_obj")

  if (is.null(lr_db)) {
    lr_db <- lr_database
  }
  lr_db <- .prepare_matrisome_pairs(lr_db)

  assay <- .resolve_assay(seurat_obj, assay)
  expr_matrix_sparse <- .safe_get_assay_data(seurat_obj, assay = assay, slot = layer)
  if (is.null(expr_matrix_sparse)) {
    stop(
      "Could not retrieve expression data from assay '", assay,
      "' using layer '", layer, "'."
    )
  }

  # Build adjacency matrix if not provided
  if (is.null(adj_matrix)) {
    .message_if(verbose, "Building spatial adjacency matrix...")
    coords <- .safe_get_tissue_coordinates(seurat_obj)

    # Use k-nearest neighbors approach
    dist_matrix <- as.matrix(stats::dist(coords[, c("imagerow", "imagecol")]))
    adj_matrix <- Matrix::Matrix(0, nrow = nrow(dist_matrix), ncol = ncol(dist_matrix),
                                 dimnames = list(rownames(coords), rownames(coords)),
                                 sparse = TRUE)

    # Connect 6 nearest neighbors (Visium hexagonal grid)
    for (i in seq_len(nrow(dist_matrix))) {
      neighbors <- order(dist_matrix[i, ])[2:7]  # Skip self (index 1)
      adj_matrix[i, neighbors] <- 1
    }
  }

  # Ensure adjacency matrix matches expression matrix
  if (!all(colnames(expr_matrix_sparse) == colnames(adj_matrix))) {
    adj_matrix <- adj_matrix[colnames(expr_matrix_sparse), colnames(expr_matrix_sparse)]
  }

  # Keep pairs whose two genes are present in the expression data
  genes_in_data <- rownames(expr_matrix_sparse)
  lr_db_filtered <- lr_db[lr_db$Ligand %in% genes_in_data &
                            lr_db$Receptor %in% genes_in_data, ]

  if (nrow(lr_db_filtered) == 0) {
    stop("No valid matrisome pairs were found in the expression data.")
  }

  .message_if(verbose, sprintf("Scoring %d valid matrisome pairs...", nrow(lr_db_filtered)))

  # Convert to dense for loop efficiency
  expr_matrix_dense <- as.matrix(expr_matrix_sparse)
  adj_matrix_dense <- as.matrix(adj_matrix)

  # Core calculation
  list_of_scores <- vector("list", nrow(lr_db_filtered))

  for (i in seq_len(nrow(lr_db_filtered))) {
    pair_row <- lr_db_filtered[i, ]
    gene_1_expr <- expr_matrix_dense[pair_row$Ligand, ]
    gene_2_expr <- expr_matrix_dense[pair_row$Receptor, ]

    # Square root transformation
    sqrt_gene_1_expr <- sqrt(pmax(gene_1_expr, 0))
    sqrt_gene_2_expr <- sqrt(pmax(gene_2_expr, 0))

    # Sum weighted expression of the first pair member across neighbors
    sum_sqrt_gene_1_neighbors <- as.vector(adj_matrix_dense %*% sqrt_gene_1_expr)

    # Final score
    list_of_scores[[i]] <- sqrt_gene_2_expr * sum_sqrt_gene_1_neighbors

    if (isTRUE(verbose) && i %% 500 == 0) {
      message(sprintf("  Processed %d / %d matrisome pairs", i, nrow(lr_db_filtered)))
    }
  }

  # Assemble results
  activity_matrix <- do.call(rbind, list_of_scores)
  dimnames(activity_matrix) <- list(
    paste(lr_db_filtered$Ligand, lr_db_filtered$Receptor, sep = "-"),
    colnames(expr_matrix_dense)
  )

  sparse_activity_matrix <- methods::as(activity_matrix, "sparseMatrix")

  seurat_obj[["LRACTIVITY"]] <- Seurat::CreateAssayObject(counts = sparse_activity_matrix)
  seurat_obj <- Seurat::SetAssayData(
    seurat_obj,
    assay = "LRACTIVITY",
    layer = "data",
    new.data = sparse_activity_matrix
  )

  .message_if(verbose, "Done!")
  return(seurat_obj)
}

#' Find enriched matrisome-pair co-expression per group
#'
#' Calculates enrichment statistics (avg_log2FC and percentage difference)
#' for each matrisome pair across cell or spot groups.
#'
#' @param seurat_obj A Seurat object with LRACTIVITY assay
#' @param assay_name Assay name (default: "LRACTIVITY")
#' @param group_by Metadata column for grouping (e.g., "ecm_niche")
#'
#' @return Data frame with columns: feature, cluster, avg_log2FC, perc_difference, pct.1, pct.2
#'
#' @examples
#' \dontrun{
#' lr_stats <- find_lr_enrichment(seurat_obj, group_by = "ecm_niche")
#' }
#'
#' @export
find_lr_enrichment <- function(seurat_obj, assay_name = "LRACTIVITY", group_by) {
  .assert_seurat_object(seurat_obj, "seurat_obj")

  if (!group_by %in% colnames(seurat_obj@meta.data)) {
    stop("group_by '", group_by, "' not found in metadata")
  }

  score_matrix <- Seurat::GetAssayData(seurat_obj, assay = assay_name, layer = "data")
  groups <- seurat_obj[[group_by, drop = TRUE]]
  cell_groups <- split(colnames(score_matrix), groups)

  all_cluster_stats <- lapply(names(cell_groups), function(cluster_name) {
    group_1_cells <- cell_groups[[cluster_name]]
    group_2_cells <- unlist(cell_groups[names(cell_groups) != cluster_name], use.names = FALSE)

    epsilon <- 1e-9
    group_1_avg <- Matrix::rowMeans(score_matrix[, group_1_cells, drop = FALSE])
    group_2_avg <- Matrix::rowMeans(score_matrix[, group_2_cells, drop = FALSE])
    avg_log2FC <- log2((group_1_avg + epsilon) / (group_2_avg + epsilon))

    pct_1 <- Matrix::rowMeans(score_matrix[, group_1_cells, drop = FALSE] > 0)
    pct_2 <- Matrix::rowMeans(score_matrix[, group_2_cells, drop = FALSE] > 0)
    perc_difference <- pct_1 - pct_2

    data.frame(
      feature = rownames(score_matrix),
      cluster = cluster_name,
      avg_log2FC = avg_log2FC,
      perc_difference = perc_difference,
      pct.1 = pct_1,
      pct.2 = pct_2
    )
  })

  dplyr::bind_rows(all_cluster_stats)
}

#' Aggregate reciprocal rows into unique matrisome pairs
#'
#' Identifies reciprocal rows (e.g., A-B and B-A) and combines them into one
#' unordered heterotypic pair. The legacy function name is retained for
#' backward compatibility.
#'
#' @param stats_df Long-format stats from \code{find_lr_enrichment}
#' @param means_matrix Wide matrix with features as rows, clusters as columns
#'
#' @return List with coalesced \code{stats} data frame and \code{means} matrix
#'
#' @examples
#' \dontrun{
#' means_matrix <- AverageExpression(seurat_obj, assays = "LRACTIVITY")$LRACTIVITY
#' coalesced <- aggregate_lr_axes(lr_stats, means_matrix)
#' }
#'
#' @export
aggregate_lr_axes <- function(stats_df, means_matrix) {

  if (nrow(stats_df) == 0 || nrow(means_matrix) == 0) {
    warning("Empty stats_df or means_matrix")
    return(list(
      stats = data.frame(feature = character(), cluster = character(),
                         avg_log2FC = numeric(), perc_difference = numeric()),
      means = matrix(nrow = 0, ncol = ncol(means_matrix),
                     dimnames = list(NULL, colnames(means_matrix)))
    ))
  }

  # Standardize cluster names
  stats_df$cluster <- gsub(" |\\-", "_", stats_df$cluster)
  colnames(means_matrix) <- gsub(" |\\-", "_", colnames(means_matrix))

  # Create canonical feature names (sorted alphabetically)
  original_features <- rownames(means_matrix)
  canonical_names <- sapply(original_features, function(name) {
    paste(sort(strsplit(name, "-")[[1]]), collapse = "-")
  })

  # Use canonical names without implying signaling direction or communication
  new_feature_names <- canonical_names

  name_map <- data.frame(
    original_feature = original_features,
    new_feature = new_feature_names
  )

  # Coalesce means matrix
  coalesced_means_matrix <- rowsum(as.matrix(means_matrix), group = new_feature_names)

  # Coalesce statistics
  coalesced_stats <- stats_df %>%
    dplyr::left_join(name_map, by = c("feature" = "original_feature"))

  if (!"new_feature" %in% colnames(coalesced_stats) ||
      all(is.na(coalesced_stats$new_feature))) {
    coalesced_stats$new_feature <- coalesced_stats$feature
  }

  coalesced_stats <- coalesced_stats %>%
    dplyr::group_by(.data$new_feature, .data$cluster) %>%
    dplyr::summarise(
      avg_log2FC = mean(.data$avg_log2FC, na.rm = TRUE),
      perc_difference = mean(.data$perc_difference, na.rm = TRUE),
      .groups = "drop"
    )

  colnames(coalesced_stats)[colnames(coalesced_stats) == "new_feature"] <- "feature"

  list(stats = coalesced_stats, means = coalesced_means_matrix)
}
