# Tests for matrisome-pair co-expression functions

test_that("lr_database loads correctly", {
  data(lr_database)
  expect_s3_class(lr_database, "data.frame")
  expect_true("Ligand" %in% colnames(lr_database))
  expect_true("Receptor" %in% colnames(lr_database))
  expect_gt(nrow(lr_database), 10000)
  expect_false(any(lr_database$Ligand == lr_database$Receptor))
  canonical_pairs <- paste(
    pmin(lr_database$Ligand, lr_database$Receptor),
    pmax(lr_database$Ligand, lr_database$Receptor),
    sep = "-"
  )
  expect_false(anyDuplicated(canonical_pairs) > 0)
})

test_that("pair preparation removes homomeric and reciprocal rows", {
  pair_db <- data.frame(
    Ligand = c("GeneB", "GeneA", "GeneC"),
    Receptor = c("GeneA", "GeneB", "GeneC")
  )

  result <- .prepare_matrisome_pairs(pair_db)

  expect_equal(nrow(result), 1)
  expect_equal(result$Ligand, "GeneA")
  expect_equal(result$Receptor, "GeneB")
})

test_that("aggregate_lr_axes handles empty input", {
  empty_stats <- data.frame(
    feature = character(),
    cluster = character(),
    avg_log2FC = numeric(),
    perc_difference = numeric()
  )
  empty_means <- matrix(nrow = 0, ncol = 3,
                        dimnames = list(NULL, c("A", "B", "C")))

  result <- aggregate_lr_axes(empty_stats, empty_means)
  expect_type(result, "list")
  expect_equal(nrow(result$stats), 0)
})

test_that("aggregate_lr_axes coalesces reciprocal rows without axis terminology", {
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

  result <- aggregate_lr_axes(stats_df, means_matrix)

  # Should have two unique unordered pairs
  expect_equal(nrow(result$stats), 2)
  expect_true("GeneA-GeneB" %in% result$stats$feature)
  expect_false(any(grepl("_AXIS", result$stats$feature)))
})
