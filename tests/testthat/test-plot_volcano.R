# Tests for plot_volcano()

# ─── Helper ───────────────────────────────────────────────────────────────────

.make_volcano_results <- function() {
    data.frame(
        chrom         = rep("chr_sim", 20L),
        position      = seq(1000L, 20000L, by = 1000L),
        strand        = rep("+", 20L),
        mod_type      = rep("6mA", 20L),
        mod_context   = rep("GATC", 20L),
        dm_pvalue     = c(runif(15L, 0, 0.1), runif(5L, 0.1, 1)),
        dm_padj       = c(runif(5L, 0, 0.05), runif(10L, 0.05, 0.5), runif(5L, 0.5, 1)),
        dm_delta_beta = c(runif(5L, 0.2, 0.8), runif(5L, -0.8, -0.2),
                          runif(10L, -0.1, 0.1)),
        stringsAsFactors = FALSE
    )
}

.make_volcano_results_multicontext <- function() {
    res <- rbind(.make_volcano_results(), .make_volcano_results())
    res$mod_context <- rep(c("GATC", "CCWGG"), each = 20L)
    res
}

# ─── Data mapping ─────────────────────────────────────────────────────────────

test_that("plot_volcano: p$data contains exact non-NA padj rows with computed neg_log10_padj", {
    res <- .make_volcano_results()
    p <- plot_volcano(res)
    expect_s3_class(p, "ggplot")
    # p$data should have exactly the non-NA padj rows
    n_nonNA <- sum(!is.na(res$dm_padj))
    expect_equal(nrow(p$data), n_nonNA)
    # neg_log10_padj should be -log10(padj) clamped
    expected_y <- -log10(pmax(res$dm_padj[!is.na(res$dm_padj)], .Machine$double.xmin))
    expect_equal(p$data$neg_log10_padj, expected_y)
    # dm_delta_beta should match input exactly
    expect_equal(p$data$dm_delta_beta, res$dm_delta_beta[!is.na(res$dm_padj)])
})

test_that("plot_volcano: significance categories match threshold rules", {
    res <- .make_volcano_results()
    p <- plot_volcano(res, delta_beta_threshold = 0.3, padj_threshold = 0.01)
    # Verify significance assignments
    pd <- p$data
    # Rows with padj <= 0.01 and delta_beta >= 0.3 should be "Hypermethylated"
    hyper <- pd$dm_padj <= 0.01 & pd$dm_delta_beta >= 0.3
    if (any(hyper)) {
        expect_equal(as.character(pd$significance[hyper]), rep("Hypermethylated", sum(hyper)))
    }
    # Rows with padj <= 0.01 and delta_beta <= -0.3 should be "Hypomethylated"
    hypo <- pd$dm_padj <= 0.01 & pd$dm_delta_beta <= -0.3
    if (any(hypo)) {
        expect_equal(as.character(pd$significance[hypo]), rep("Hypomethylated", sum(hypo)))
    }
})

test_that("plot_volcano: vlines at exact +/- delta_beta_threshold", {
    res <- .make_volcano_results()
    p <- plot_volcano(res, delta_beta_threshold = 0.3, padj_threshold = 0.01)
    layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    vline_idx <- which(layer_classes == "GeomVline")
    expect_gte(length(vline_idx), 2L)
    # Extract exact intercept values
    intercepts <- vapply(vline_idx, function(i) {
        layer_data <- p$layers[[i]]$data
        if (is.function(layer_data)) NA_real_ else layer_data$xintercept[1]
    }, numeric(1))
    intercepts <- intercepts[!is.na(intercepts)]
    # Both +0.3 and -0.3 must be present
    expect_true(any(abs(intercepts - 0.3) < 1e-10))
    expect_true(any(abs(intercepts + 0.3) < 1e-10))
})

# ─── Threshold lines ──────────────────────────────────────────────────────────

test_that("plot_volcano: NULL delta_beta_threshold has 2 layers (points + hline)", {
    res <- .make_volcano_results()
    p <- plot_volcano(res)
    expect_gte(length(p$layers), 2L)
})

test_that("plot_volcano: numeric delta_beta_threshold adds vlines (>= 4 layers)", {
    res <- .make_volcano_results()
    p <- plot_volcano(res, delta_beta_threshold = 0.2)
    expect_gte(length(p$layers), 4L)
})

# ─── Faceting ─────────────────────────────────────────────────────────────────

test_that("plot_volcano: single mod_context does not facet", {
    res <- .make_volcano_results()
    p <- plot_volcano(res)
    expect_false(inherits(p$facet, "FacetWrap"))
})

test_that("plot_volcano: multiple mod_context levels facet when facet = TRUE", {
    res <- .make_volcano_results_multicontext()
    p <- plot_volcano(res, facet = TRUE)
    expect_true(inherits(p$facet, "FacetWrap"))
})

test_that("plot_volcano: facet = FALSE suppresses faceting for multi-context results", {
    res <- .make_volcano_results_multicontext()
    p <- plot_volcano(res, facet = FALSE)
    expect_false(inherits(p$facet, "FacetWrap"))
})

test_that("plot_volcano: error on invalid facet argument", {
    res <- .make_volcano_results()
    expect_error(plot_volcano(res, facet = "yes"), "facet")
    expect_error(plot_volcano(res, facet = NA), "facet")
})

# ─── NA handling ─────────────────────────────────────────────────────────────

test_that("plot_volcano: rows with NA padj are excluded from p$data exactly", {
    res <- .make_volcano_results()
    res$dm_padj[1:5] <- NA
    p <- plot_volcano(res)
    expect_s3_class(p, "ggplot")
    # p$data should have exactly the non-NA padj rows
    n_nonNA <- sum(!is.na(res$dm_padj))
    expect_equal(nrow(p$data), n_nonNA)
    # No NA padj values should appear in p$data
    expect_false(any(is.na(p$data$dm_padj)))
})

test_that("plot_volcano: rows with NA delta_beta are kept in p$data", {
    res <- .make_volcano_results()
    res$dm_delta_beta[1:3] <- NA
    p <- plot_volcano(res)
    expect_s3_class(p, "ggplot")
    # NA delta_beta rows are still plotted (only NA padj excluded)
    n_nonNA_padj <- sum(!is.na(res$dm_padj))
    expect_equal(nrow(p$data), n_nonNA_padj)
    # NA delta_beta values should be present in p$data
    expect_true(any(is.na(p$data$dm_delta_beta)))
})

# ─── Error conditions ─────────────────────────────────────────────────────────

test_that("plot_volcano: error on non-data.frame input", {
    expect_error(plot_volcano(list(dm_delta_beta = 1, dm_padj = 0.01)),
                 "data.frame")
})

test_that("plot_volcano: error when dm_delta_beta column is missing", {
    res <- .make_volcano_results()
    res$dm_delta_beta <- NULL
    expect_error(plot_volcano(res), "dm_delta_beta")
})

test_that("plot_volcano: error when dm_padj column is missing", {
    res <- .make_volcano_results()
    res$dm_padj <- NULL
    expect_error(plot_volcano(res), "dm_padj")
})

test_that("plot_volcano: error when all padj values are NA", {
    res <- .make_volcano_results()
    res$dm_padj <- NA_real_
    expect_error(plot_volcano(res))
})

test_that("plot_volcano: NULL delta_beta_threshold colors by padj alone", {
    set.seed(42L)
    res <- .make_volcano_results()
    p <- plot_volcano(res, delta_beta_threshold = NULL)
    # With NULL threshold, significance is based on padj only + sign of delta_beta
    pd <- p$data
    sig_pos <- pd$dm_padj <= 0.05 & pd$dm_delta_beta > 0
    if (any(sig_pos)) {
        expect_equal(as.character(pd$significance[sig_pos]), rep("Hypermethylated", sum(sig_pos)))
    }
    sig_neg <- pd$dm_padj <= 0.05 & pd$dm_delta_beta < 0
    if (any(sig_neg)) {
        expect_equal(as.character(pd$significance[sig_neg]), rep("Hypomethylated", sum(sig_neg)))
    }
})

test_that("plot_volcano: error on invalid delta_beta_threshold", {
    res <- .make_volcano_results()
    expect_error(plot_volcano(res, delta_beta_threshold = 0), "delta_beta_threshold")
    expect_error(plot_volcano(res, delta_beta_threshold = 1), "delta_beta_threshold")
    expect_error(plot_volcano(res, delta_beta_threshold = -0.1), "delta_beta_threshold")
})

test_that("plot_volcano: error on invalid padj_threshold", {
    res <- .make_volcano_results()
    expect_error(plot_volcano(res, padj_threshold = 0), "padj_threshold")
    expect_error(plot_volcano(res, padj_threshold = 1), "padj_threshold")
})

# ─── Comma example data ───────────────────────────────────────────────────────

test_that("plot_volcano: works with results() output from comma_example_data", {
    data(comma_example_data)
    cd_dm <- diffMethyl(comma_example_data, ~ condition, mod_type = "6mA")
    res <- results(cd_dm)
    p <- plot_volcano(res)
    expect_s3_class(p, "ggplot")
    # p$data row count should match non-NA padj rows in results
    n_nonNA <- sum(!is.na(res$dm_padj))
    expect_equal(nrow(p$data), n_nonNA)
})
