# Tests for matrisome-pair co-expression functions

test_that("matrisome_pairs loads correctly", {
  data(matrisome_pairs)
  expect_s3_class(matrisome_pairs, "data.frame")
  expect_true("Gene1" %in% colnames(matrisome_pairs))
  expect_true("Gene2" %in% colnames(matrisome_pairs))
  expect_gt(nrow(matrisome_pairs), 10000)
  expect_false(any(matrisome_pairs$Gene1 == matrisome_pairs$Gene2))
  canonical_pairs <- paste(
    pmin(matrisome_pairs$Gene1, matrisome_pairs$Gene2),
    pmax(matrisome_pairs$Gene1, matrisome_pairs$Gene2),
    sep = "-"
  )
  expect_false(anyDuplicated(canonical_pairs) > 0)
})

test_that("pair preparation removes homomeric and reciprocal rows", {
  pair_db <- data.frame(
    Gene1 = c("GeneB", "GeneA", "GeneC"),
    Gene2 = c("GeneA", "GeneB", "GeneC")
  )

  result <- .prepare_matrisome_pairs(pair_db)

  expect_equal(nrow(result), 1)
  expect_equal(result$Gene1, "GeneA")
  expect_equal(result$Gene2, "GeneB")
})

test_that("matrisome-pair scoring creates the renamed assay", {
  counts <- Matrix::Matrix(
    matrix(
      c(4, 0, 1, 0, 9, 1),
      nrow = 2,
      dimnames = list(c("GeneA", "GeneB"), c("Spot1", "Spot2", "Spot3"))
    ),
    sparse = TRUE
  )
  seurat_obj <- Seurat::CreateSeuratObject(counts = counts)
  seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)
  adjacency <- Matrix::Matrix(
    matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3),
    dimnames = list(colnames(counts), colnames(counts)),
    sparse = TRUE
  )

  result <- score_matrisome_pair_coexpression(
    seurat_obj,
    pair_db = data.frame(Gene1 = "GeneA", Gene2 = "GeneB"),
    adj_matrix = adjacency,
    assay = "RNA",
    verbose = FALSE
  )

  expect_true("MATRISOMEPAIR" %in% SeuratObject::Assays(result))
  expect_equal(rownames(result[["MATRISOMEPAIR"]]), "GeneA-GeneB")
})

test_that("renamed pair-analysis functions are exported", {
  pair_functions <- c(
    "score_matrisome_pair_coexpression",
    "find_matrisome_pair_enrichment",
    "aggregate_matrisome_pairs"
  )

  expect_true(all(pair_functions %in% getNamespaceExports("matrispace")))
})

test_that("aggregate_matrisome_pairs handles empty input", {
  empty_stats <- data.frame(
    feature = character(),
    cluster = character(),
    avg_log2FC = numeric(),
    perc_difference = numeric()
  )
  empty_means <- matrix(nrow = 0, ncol = 3,
                        dimnames = list(NULL, c("A", "B", "C")))

  result <- aggregate_matrisome_pairs(empty_stats, empty_means)
  expect_type(result, "list")
  expect_equal(nrow(result$stats), 0)
})

test_that("aggregate_matrisome_pairs coalesces reciprocal rows", {
  # Create test data with reciprocal pair
  stats_df <- data.frame(
    feature = c("GeneA-GeneB", "GeneB-GeneA", "GeneC-GeneD"),
    cluster = c("Cluster1", "Cluster1", "Cluster1"),
    avg_log2FC = c(1.0, 0.8, 0.5),
    perc_difference = c(0.2, 0.15, 0.1)
  )

  means_matrix <- matrix(
    c(1, 0.5, 0.3),
    nrow = 3, ncol = 1,
    dimnames = list(c("GeneA-GeneB", "GeneB-GeneA", "GeneC-GeneD"), "Cluster1")
  )

  result <- aggregate_matrisome_pairs(stats_df, means_matrix)

  # Should have two unique unordered pairs
  expect_equal(nrow(result$stats), 2)
  expect_true("GeneA-GeneB" %in% result$stats$feature)
  expect_false(any(grepl("_AXIS", result$stats$feature)))
})
