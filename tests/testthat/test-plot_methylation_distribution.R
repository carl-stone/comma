# Tests for plot_methylation_distribution()

# ─── Helper ───────────────────────────────────────────────────────────────────

.make_dist_data <- function() {
    n_sites   <- 10L
    positions <- seq(1000L, 10000L, by = 1000L)
    set.seed(1L)
    betas <- matrix(
        runif(n_sites * 3L, 0.1, 0.9),
        nrow = n_sites, ncol = 3L,
        dimnames = list(NULL, c("ctrl_1", "ctrl_2", "treat_1"))
    )
    cov_mat <- matrix(20L, nrow = n_sites, ncol = 3L,
                      dimnames = dimnames(betas))
    site_gr <- GenomicRanges::GRanges(
        seqnames = rep("chr_sim", n_sites),
        ranges   = IRanges::IRanges(start = positions, width = 1L),
        strand   = rep("+", n_sites),
        mod_type    = factor(rep("6mA", n_sites), levels = c("4mC", "5mC", "6mA")),
        motif       = rep("GATC", n_sites)
    )
    GenomeInfoDb::seqinfo(site_gr) <- GenomeInfoDb::Seqinfo(
        seqnames = "chr_sim",
        seqlengths = 100000L,
        isCircular = FALSE
    )
    cd <- S4Vectors::DataFrame(
        sample_name = c("ctrl_1", "ctrl_2", "treat_1"),
        condition   = c("control", "control", "treatment"),
        replicate   = 1:3,
        row.names   = c("ctrl_1", "ctrl_2", "treat_1")
    )
    rse <- SummarizedExperiment::SummarizedExperiment(
        assays     = list(methylation = betas, coverage = cov_mat),
        rowRanges  = site_gr,
        colData    = cd
    )
    new("commaData", rse)
}

## Object with two modification types
.make_dist_data_two_mods <- function() {
    n_6ma <- 8L; n_5mc <- 4L
    n_sites <- n_6ma + n_5mc
    positions <- seq(1000L, n_sites * 1000L, by = 1000L)
    mod_types  <- factor(c(rep("6mA", n_6ma), rep("5mC", n_5mc)),
                         levels = c("4mC", "5mC", "6mA"))
    motif_vals <- c(rep("GATC", n_6ma), rep("CCWGG", n_5mc))
    set.seed(2L)
    betas <- matrix(
        runif(n_sites * 2L, 0.1, 0.9),
        nrow = n_sites, ncol = 2L,
        dimnames = list(NULL, c("samp1", "samp2"))
    )
    cov_mat <- matrix(20L, nrow = n_sites, ncol = 2L,
                      dimnames = dimnames(betas))
    site_gr <- GenomicRanges::GRanges(
        seqnames = rep("chr_sim", n_sites),
        ranges   = IRanges::IRanges(start = positions, width = 1L),
        strand   = rep("+", n_sites),
        mod_type    = mod_types,
        motif       = motif_vals
    )
    GenomeInfoDb::seqinfo(site_gr) <- GenomeInfoDb::Seqinfo(
        seqnames = "chr_sim",
        seqlengths = 100000L,
        isCircular = FALSE
    )
    cd <- S4Vectors::DataFrame(
        sample_name = c("samp1", "samp2"),
        condition   = c("ctrl", "treat"),
        replicate   = 1:2,
        row.names   = c("samp1", "samp2")
    )
    rse <- SummarizedExperiment::SummarizedExperiment(
        assays     = list(methylation = betas, coverage = cov_mat),
        rowRanges  = site_gr,
        colData    = cd
    )
    new("commaData", rse)
}

# ─── Basic return type ────────────────────────────────────────────────────────

test_that("plot_methylation_distribution: returns ggplot for valid input", {
    obj <- .make_dist_data()
    p <- plot_methylation_distribution(obj)
    expect_s3_class(p, "ggplot")
    # Verify density layer has data points (beta values)
    bd <- ggplot2::ggplot_build(p)$data[[1]]
    expect_true(nrow(bd) > 0L)
    # Each density group should have observations (n > 0)
    expect_true(all(bd$n > 0L))
    # x values should be in [0, 1] (beta range)
    expect_true(all(bd$x >= 0 & bd$x <= 1))
})

test_that("plot_methylation_distribution: mod_type filter returns ggplot", {
    obj <- .make_dist_data_two_mods()
    p_unfiltered <- plot_methylation_distribution(obj)
    p <- plot_methylation_distribution(obj, mod_type = "6mA")
    expect_s3_class(p, "ggplot")
    # Verify filtered data has different content than unfiltered
    bd_filtered   <- ggplot2::ggplot_build(p)$data[[1]]
    bd_unfiltered <- ggplot2::ggplot_build(p_unfiltered)$data[[1]]
    n_filtered   <- sum(bd_filtered$n[!duplicated(paste(bd_filtered$group,
                                                        bd_filtered$PANEL))])
    n_unfiltered <- sum(bd_unfiltered$n[!duplicated(paste(bd_unfiltered$group,
                                                          bd_unfiltered$PANEL))])
    expect_true(n_filtered < n_unfiltered)
})

# ─── Faceting ─────────────────────────────────────────────────────────────────

test_that("plot_methylation_distribution: multi-mod object produces facets", {
    obj <- .make_dist_data_two_mods()
    p <- plot_methylation_distribution(obj)
    expect_s3_class(p, "ggplot")
    # facet_wrap wraps in a FacetWrap layer, not FacetNull
    expect_false(inherits(p$facet, "FacetNull"))
})

test_that("plot_methylation_distribution: single-mod object has no facets", {
    obj <- .make_dist_data()
    p <- plot_methylation_distribution(obj)
    expect_true(inherits(p$facet, "FacetNull"))
})

# ─── NA handling ─────────────────────────────────────────────────────────────

test_that("plot_methylation_distribution: NAs in beta values are silently excluded", {
    obj <- .make_dist_data()
    p_orig <- plot_methylation_distribution(obj)
    # Inject NAs into the methylation matrix
    methyl_mat <- methylation(obj)
    methyl_mat[1:3, "ctrl_1"] <- NA
    SummarizedExperiment::assay(obj, "methylation") <- methyl_mat
    p <- plot_methylation_distribution(obj)
    expect_s3_class(p, "ggplot")
    # Verify fewer data points after NA injection
    bd_orig <- ggplot2::ggplot_build(p_orig)$data[[1]]
    bd_na   <- ggplot2::ggplot_build(p)$data[[1]]
    n_orig <- sum(bd_orig$n[!duplicated(paste(bd_orig$group,
                                              bd_orig$PANEL))])
    n_na   <- sum(bd_na$n[!duplicated(paste(bd_na$group,
                                            bd_na$PANEL))])
    expect_true(n_na < n_orig)
})

# ─── Error conditions ─────────────────────────────────────────────────────────

test_that("plot_methylation_distribution: error on non-commaData input", {
    expect_error(plot_methylation_distribution(data.frame(x = 1)),
                 "commaData")
})

test_that("plot_methylation_distribution: error on invalid mod_type", {
    obj <- .make_dist_data()
    expect_error(plot_methylation_distribution(obj, mod_type = "4mC"),
                 "not found")
})

test_that("plot_methylation_distribution: error on invalid per_sample", {
    obj <- .make_dist_data()
    expect_error(plot_methylation_distribution(obj, per_sample = "yes"),
                 "per_sample")
})

# ─── Comma example data ───────────────────────────────────────────────────────

test_that("plot_methylation_distribution: works with comma_example_data", {
    data(comma_example_data)
    p <- plot_methylation_distribution(comma_example_data)
    expect_s3_class(p, "ggplot")
})


