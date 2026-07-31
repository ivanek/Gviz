test_that("BiomartGeneRegionTrack works", {
    ## expect_s4_class(BiomartGeneRegionTrack(), "BiomartGeneRegionTrack")

    biomartMapping <- list(
        gene_id = "ensembl_gene_id",
        transcript_id = "ensembl_transcript_id",
        exon_id = "ensembl_exon_id",
        start = "exon_chrom_start",
        end = "exon_chrom_end",
        rank = "rank",
        strand = "strand",
        symbol = c("external_gene_name", "external_gene_id"),
        feature = "gene_biotype",
        chromosome = "chromosome_name",
        u5s = "5_utr_start",
        u5e = "5_utr_end",
        u3s = "3_utr_start",
        u3e = "3_utr_end",
        cdsl = c("cds_length", "cds_start"),
        phase = "phase"
    )
    expect_identical(.getBMFeatureMap(), biomartMapping)
})
