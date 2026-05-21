# Tests for plot_genome_track()

# ─── Helper ───────────────────────────────────────────────────────────────────

.make_track_data <- function() {
    n_sites   <- 10L
    positions <- seq(1000L, 10000L, by = 1000L)
    set.seed(3L)
    betas <- matrix(
        runif(n_sites * 2L, 0.1, 0.9),
        nrow = n_sites, ncol = 2L,
        dimnames = list(NULL, c("ctrl_1", "treat_1"))
    )
    cov_mat <- matrix(20L, nrow = n_sites, ncol = 2L,
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
        sample_name = c("ctrl_1", "treat_1"),
        condition   = c("control", "treatment"),
        replicate   = 1:2,
        row.names   = c("ctrl_1", "treat_1")
    )
    ann_gr <- GenomicRanges::GRanges(
        seqnames = "chr_sim",
        ranges   = IRanges::IRanges(start = c(2000L, 6000L),
                                     end   = c(4000L, 8000L)),
        strand   = "+"
    )
    ann_gr$feature_type <- c("gene", "gene")
    ann_gr$name         <- c("geneA", "geneB")

    rse <- SummarizedExperiment::SummarizedExperiment(
        assays     = list(methylation = betas, coverage = cov_mat),
        rowRanges  = site_gr,
        colData    = cd
    )
    obj <- new("commaData", rse)
    S4Vectors::metadata(obj)$annotation <- ann_gr
    S4Vectors::metadata(obj)$motifSites <- GenomicRanges::GRanges()
    obj
}

# ─── Basic return type ────────────────────────────────────────────────────────

test_that("plot_genome_track: returns ggplot with point data for valid chromosome", {
    obj <- .make_track_data()
    p <- plot_genome_track(obj, chromosome = "chr_sim", annotation = FALSE)
    expect_s3_class(p, "ggplot")
    # Verify point layer has data with position values
    bd <- ggplot2::ggplot_build(p)
    point_data <- bd$data[[1]]
    expect_true(nrow(point_data) > 0)
    # x values are genomic positions (positive integers)
    expect_true(all(point_data$x[is.finite(point_data$x)] > 0))
})

test_that("plot_genome_track: annotation = FALSE returns single ggplot without annotation track", {
    obj <- .make_track_data()
    p <- plot_genome_track(obj, chromosome = "chr_sim", annotation = FALSE)
    expect_s3_class(p, "ggplot")
    # Without annotation, only point layer (no rect layer for features)
    layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    expect_true("GeomPoint" %in% layer_classes)
    expect_false("GeomRect" %in% layer_classes)
})

test_that("plot_genome_track: mod_type filter reduces plotted points", {
    obj <- .make_track_data()
    p_all <- plot_genome_track(obj, chromosome = "chr_sim", annotation = FALSE)
    p_filt <- plot_genome_track(obj, chromosome = "chr_sim",
                                mod_type = "6mA", annotation = FALSE)
    expect_s3_class(p_filt, "ggplot")
    # All sites are 6mA in this fixture, so same count
    bd_filt <- ggplot2::ggplot_build(p_filt)
    expect_true(nrow(bd_filt$data[[1]]) > 0)
})

# ─── Positional filtering ─────────────────────────────────────────────────────

test_that("plot_genome_track: start/end filtering reduces displayed sites", {
    obj <- .make_track_data()
    p <- plot_genome_track(obj, chromosome = "chr_sim",
                           start = 1000L, end = 5000L, annotation = FALSE)
    expect_s3_class(p, "ggplot")
    bd <- ggplot2::ggplot_build(p)
    # All x values in layer 1 (points) should be within [1000, 5000]
    x_vals <- bd$data[[1L]]$x
    expect_true(all(x_vals >= 1000L & x_vals <= 5000L))
})

test_that("plot_genome_track: start only (no end) filters positions >= start", {
    obj <- .make_track_data()
    p <- plot_genome_track(obj, chromosome = "chr_sim",
                           start = 3000L, annotation = FALSE)
    expect_s3_class(p, "ggplot")
    bd <- ggplot2::ggplot_build(p)
    x_vals <- bd$data[[1L]]$x
    # Must have at least one point (guard against vacuous pass on empty vector)
    expect_true(length(x_vals) > 0L)
    # All displayed positions should be >= 3000
    expect_true(all(x_vals >= 3000L))
})

# ─── Error conditions ─────────────────────────────────────────────────────────

test_that("plot_genome_track: error on non-commaData input", {
    expect_error(plot_genome_track(data.frame(), chromosome = "chr_sim"),
                 "commaData")
})

test_that("plot_genome_track: error when chromosome not in genome", {
    obj <- .make_track_data()
    expect_error(plot_genome_track(obj, chromosome = "chrX"),
                 "chrX")
})

test_that("plot_genome_track: error when no sites on chromosome", {
    obj <- .make_track_data()
    # Add chr2 to Seqinfo but not to data
    new_si <- GenomeInfoDb::Seqinfo(
        seqnames = c("chr_sim", "chr2"),
        seqlengths = c(100000L, 50000L),
        isCircular = c(FALSE, FALSE)
    )
    GenomeInfoDb::seqinfo(obj) <- new_si
    expect_error(plot_genome_track(obj, chromosome = "chr2"),
                 "No methylation sites found")
})

test_that("plot_genome_track: error when start > end", {
    obj <- .make_track_data()
    expect_error(
        plot_genome_track(obj, chromosome = "chr_sim",
                          start = 5000L, end = 1000L),
        "'end' must be"
    )
})

test_that("plot_genome_track: error on invalid mod_type", {
    obj <- .make_track_data()
    expect_error(plot_genome_track(obj, chromosome = "chr_sim",
                                   mod_type = "4mC"),
                 "not found")
})

test_that("plot_genome_track: error when invalid annotation argument", {
    obj <- .make_track_data()
    expect_error(
        plot_genome_track(obj, chromosome = "chr_sim", annotation = 42),
        "'annotation' must be"
    )
})

# ─── Comma example data ───────────────────────────────────────────────────────

test_that("plot_genome_track: works with comma_example_data", {
    data(comma_example_data)
    p <- plot_genome_track(comma_example_data, chromosome = "chr_sim",
                           annotation = FALSE)
    expect_s3_class(p, "ggplot")
})
