#' @include GeneRegionTrack-class.R
NULL

## BiomartGeneRegionTrack Class ----------------------------------------------

#' BiomartGeneRegionTrack class and methods
#'
#'
#' A class to hold gene model data for a genomic region fetched dynamically
#' from EBI's Biomart Ensembl data source.
#'
#'
#' A track containing all gene models in a particular region as fetched from
#' EBI's Biomart service. Usually the user does not have to take care of the
#' Biomart connection, which will be established automatically based on the
#' provided genome and chromosome information. However, for full flexibility a
#' valid \code{\linkS4class{Mart}} object may be passed on to the constructor.
#' Please note that this assumes a connection to one of the Ensembl gene data
#' sources, mapping the available query data back to the internal object slots.
#'
#' @template BiomartGeneRegionTrack-class_param
#'
#' @name BiomartGeneRegionTrack-class
#'
#' @return
#'
#' The return value of the constructor function is a new object of class
#' \code{BiomartGeneRegionTrack}.
#' @section Objects from the class:
#'
#' Objects can be created using the constructor function
#' \code{BiomartGeneRegionTrack}.
#'
#' @author Florian Hahne
#'
#' @inherit GdObject-class seealso
#'
#' @references EBI Biomart webservice at \url{http://www.biomart.org}.
#'
#' @examples
#' \dontshow{
#' ## Load some sample data
#' data(bmTrack)
#' }
#'
#' ## Construct the object
#' \dontrun{
#' bmTrack <- BiomartGeneRegionTrack(
#'     start = 26682683, end = 26711643,
#'     chromosome = 7, genome = "mm9"
#' )
#' }
#'
#' \dontshow{
#' ## For some annoying reason the postscript device does not know about
#' ## the sans font
#' if (!interactive()) {
#'     font <- ps.options()$family
#'     displayPars(bmTrack) <- list(fontfamily = font, fontfamily.title = font)
#' }
#' }
#'
#' ## Plotting
#' plotTracks(bmTrack)
#'
#' ## Track names
#' names(bmTrack)
#' names(bmTrack) <- "foo"
#' plotTracks(bmTrack)
#'
#' ## Subsetting and splitting
#' subTrack <- subset(bmTrack, from = 26700000, to = 26705000)
#' length(subTrack)
#' subTrack <- bmTrack[transcript(bmTrack) == "ENSMUST00000144140"]
#' split(bmTrack, transcript(bmTrack))
#'
#' ## Accessors
#' start(bmTrack)
#' end(bmTrack)
#' width(bmTrack)
#' position(bmTrack)
#' width(subTrack) <- width(subTrack) + 100
#'
#' strand(bmTrack)
#' strand(subTrack) <- "-"
#'
#' chromosome(bmTrack)
#' chromosome(subTrack) <- "chrX"
#'
#' genome(bmTrack)
#' genome(subTrack) <- "hg19"
#'
#' range(bmTrack)
#' ranges(bmTrack)
#'
#' ## Annotation
#' identifier(bmTrack)
#' identifier(bmTrack, "lowest")
#' identifier(subTrack) <- "bar"
#'
#' feature(bmTrack)
#' feature(subTrack) <- "foo"
#'
#' exon(bmTrack)
#' exon(subTrack) <- letters[1:2]
#'
#' gene(bmTrack)
#' gene(subTrack) <- "bar"
#'
#' symbol(bmTrack)
#' symbol(subTrack) <- "foo"
#'
#' transcript(bmTrack)
#' transcript(subTrack) <- c("foo", "bar")
#' chromosome(subTrack) <- "chr7"
#' plotTracks(subTrack)
#'
#' values(bmTrack)
#'
#' ## Grouping
#' group(bmTrack)
#' group(subTrack) <- "Group 1"
#' transcript(subTrack)
#' plotTracks(subTrack)
#'
#' ## Stacking
#' stacking(bmTrack)
#' stacking(bmTrack) <- "dense"
#' plotTracks(bmTrack)
#'
#' ## coercion
#' as(bmTrack, "data.frame")
#' as(bmTrack, "UCSCData")
#'
#' ## HTML image map
#' coords(bmTrack)
#' tags(bmTrack)
#' bmTrack <- plotTracks(bmTrack)$foo
#' coords(bmTrack)
#' tags(bmTrack)
#' @importFrom biomaRt getBM useEnsembl useMart listDatasets listAttributes listFilters
#'
#' @exportClass BiomartGeneRegionTrack
setClass("BiomartGeneRegionTrack",
    contains = "GeneRegionTrack",
    representation = representation(biomart = "MartOrNULL", filter = "list"),
    prototype = prototype(
        biomart = NULL,
        filter = list(),
        columns = c("feature", "transcript", "symbol", "gene", "rank"),
        name = "BiomartGeneRegionTrack",
        dp = DisplayPars(
            C_segment = "burlywood4",
            D_segment = "lightblue",
            J_segment = "dodgerblue2",
            Mt_rRNA = "yellow",
            Mt_tRNA = "darkgoldenrod",
            Mt_tRNA_pseudogene = "darkgoldenrod1",
            V_segment = "aquamarine",
            miRNA = "cornflowerblue",
            miRNA_pseudogene = "cornsilk",
            misc_RNA = "cornsilk3",
            misc_RNA_pseudogene = "cornsilk4",
            protein_coding = "#FFD58A",
            pseudogene = "brown1",
            rRNA = "darkolivegreen1",
            rRNA_pseudogene = "darkolivegreen",
            retrotransposed = "blueviolet",
            scRNA = "gold4",
            scRNA_pseudogene = "darkorange2",
            snRNA = "coral",
            snRNA_pseudogene = "coral3",
            snoRNA = "cyan",
            snoRNA_pseudogene = "cyan2",
            tRNA_pseudogene = "antiquewhite3",
            utr3 = "#FFD58A",
            utr5 = "#FFD58A",
            verbose = FALSE
        )
    )
)

## Helper to return the default biomart to feature mapping
.getBMFeatureMap <- function() {
    return(list(
        gene_id = "ensembl_gene_id", transcript_id = "ensembl_transcript_id", exon_id = "ensembl_exon_id",
        start = "exon_chrom_start", end = "exon_chrom_end", rank = "rank", strand = "strand",
        symbol = c("external_gene_name", "external_gene_id"), feature = "gene_biotype", chromosome = "chromosome_name",
        u5s = "5_utr_start", u5e = "5_utr_end", u3s = "3_utr_start", u3e = "3_utr_end", cdsl = c("cds_length", "cds_start"),
        phase = "phase"
    ))
}

## Helper to do the actual fetching of data from Biomart

#' @importFrom utils modifyList
.fetchBMData <- function(object, chromosome_name = NULL, staged = FALSE) {
    if (!is.null(chromosome_name)) {
        chromosome_name <- gsub("^chr", "", chromosome_name)
    }
    ## The map between Biomart DB fields and annotation features. The individual values can be vectors for cases where there is
    ## ambiguity between different marts. This will be dynamically evaluated against available filters.
    origFeatureMap <- .getBMFeatureMap()
    featureMap <- modifyList(origFeatureMap, as.list(.dpOrDefault(object, ".__featureMap", list())))
    needed <- c(
        "gene_id", "transcript_id", "exon_id", "start", "end", "rank", "strand", "symbol", "feature",
        "chromosome", "u5s", "u5e", "u3s", "u3e", "cdsl", "phase"
    )
    if (!all(needed %in% names(featureMap))) {
        stop("'featureMap' needs to include items '", paste(setdiff(needed, names(featureMap)), collapse = ", "), "'")
    }
    avail <- listAttributes(object@biomart)[, 1]
    ambig <- names(featureMap)[listLen(featureMap) > 1]
    for (i in ambig) {
        mt <- match(featureMap[[i]], avail)
        if (!all(is.na(mt))) {
            featureMap[[i]] <- featureMap[[i]][min(which(!is.na(mt)))]
        } else {
            featureMap[[i]] <- NA
        }
    }
    featureMap <- unlist(featureMap)
    ## Deal with the filters
    filterValues <- as.list(object@filter)
    start <- object@start
    end <- object@end
    for (i in c("start", "end", "chromosome_name")) {
        if (!is.null(get(i)) && length(get(i)) > 0) {
            filterValues[[i]] <- get(i)
        }
    }
    ens <- getBM(as.vector(featureMap),
        filters = names(filterValues),
        values = filterValues, bmHeader = FALSE,
        mart = object@biomart, uniqueRows = TRUE
    )
    colnames(ens) <- names(featureMap)
    if (staged && nrow(ens) > 0) {
        filterValues <- list(start = min(ens$start), end = max(ens$end), chromosome_name = ens[1, "chromosome"])
        ens <- getBM(as.vector(featureMap),
            filters = names(filterValues),
            values = filterValues, bmHeader = FALSE,
            mart = object@biomart, uniqueRows = TRUE
        )
        colnames(ens) <- names(featureMap)
    }
    ## Only those transcripts that have a CDS length will be considered protein_coding
    ens$feature <- ifelse(is.na(ens$cdsl) & ens$feature == "protein_coding", "non_coding", ens$feature)
    ## We may have to split exons if they contain UTRs
    hasUtr <- !is.na(ens$u5s) | !is.na(ens$u3s)
    ensUtr <- ens[hasUtr, , drop = FALSE]
    ensUtr$ffeature <- ifelse(is.na(ensUtr$u5s), "utr3", "utr5")
    ensUtr$us <- ifelse(ensUtr$ffeature == "utr3", ensUtr$u3s, ensUtr$u5s)
    ensUtr$ue <- ifelse(ensUtr$ffeature == "utr3", ensUtr$u3e, ensUtr$u5e)
    ensUtr$u5e <- ensUtr$u5s <- ensUtr$u3e <- ensUtr$u3s <- NULL
    allUtr <- ensUtr$us == ensUtr$start & ensUtr$ue == ensUtr$end
    utrFinal <- ensUtr[allUtr, , drop = FALSE]
    ensUtr <- ensUtr[!allUtr, , drop = FALSE]
    ensUtrS <- split(ensUtr, ifelse(ensUtr$start == ensUtr$us, "left", "right"))
    utrFinal <- rbind(utrFinal, do.call(rbind, lapply(names(ensUtrS), function(i) {
        y <- ensUtrS[[i]]
        if (nrow(y) == 0) {
            return(NULL)
        }
        yy <- y[rep(seq_len(nrow(y)), each = 2), ]
        sel <- seq(1, nrow(yy), by = 2)
        yy[sel, "end"] <- if (i == "left") yy[sel, "ue"] else yy[sel, "us"] - 1
        yy[sel, "ffeature"] <- yy[sel, ifelse(i == "left", "ffeature", "feature")]
        yy[sel, "phase"] <- if (i == "left") -1 else 0
        sel <- seq(2, nrow(yy), by = 2)
        yy[sel, "start"] <- if (i == "left") yy[sel, "ue"] + 1 else yy[sel, "us"]
        yy[sel, "ffeature"] <- yy[sel, ifelse(i == "left", "feature", "ffeature")]
        yy[sel, "phase"] <- if (i == "left") yy[sel, "phase"] else -1
        yy
    })))
    utrFinal$feature <- utrFinal$ffeature
    keep <- c(
        "gene_id", "transcript_id", "exon_id", "start",
        "end", "rank", "strand", "symbol", "feature",
        "chromosome", "phase"
    )
    ens <- rbind(ens[!hasUtr, keep, drop = FALSE], utrFinal[, keep])
    ens$chromosome <- .chrName(ens$chromosome, force = TRUE)
    range <- GRanges(
        seqnames = ens$chromosome,
        ranges = IRanges(start = ens$start, end = ens$end),
        strand = ens$strand, feature = as.character(ens$feature),
        gene = as.character(ens$gene_id), exon = as.character(ens$exon_id),
        transcript = as.character(ens$transcript_id), symbol = as.character(ens$symbol),
        rank = as.numeric(ens$rank), phase = as.integer(ens$phase)
    )
    suppressWarnings(genome(range) <- unname(genome(object)[1]))
    range <- sort(range)
    return(range)
}


## Create an MD5 hash for a BiomartGeneRegion track taking into account the Biomart details for caching

#' @importFrom digest digest
.bmGuid <- function(bmtrack) {
    digest(list(
        genome = genome(bmtrack), host = bmtrack@biomart@host, mart = bmtrack@biomart@biomart, schema = bmtrack@biomart@vschema,
        dataset = bmtrack@biomart@dataset, filters = bmtrack@filter
    ))
}

## Initialize ----------------------------------------------------------------

# #' @noRd
# #' @keywords internal
#' @describeIn BiomartGeneRegionTrack-class Initialize.
#' @export
setMethod("initialize", "BiomartGeneRegionTrack", function(.Object, start = NULL, end = NULL, biomart, filter = list(), range, genome = NULL, chromosome = NULL, strand = NULL,
                                                           featureMap = NULL, symbol = NULL, gene = NULL, transcript = NULL, entrez = NULL, ...) {
    if ((missing(range) || is.null(range)) && is.null(genome) && is.null(chromosome)) {
        return(.Object)
    }
    ## the display parameter defaults
    .makeParMapping()
    .Object <- .updatePars(.Object, "BiomartGeneRegionTrack")
    verb <- list(...)$verbose
    displayPars(.Object) <- list(!is.null(verb) && verb == TRUE)
    ## preparing filters
    strand <- .strandName(strand, extended = TRUE)
    if (strand %in% 0:1) {
        filter$strand <- c(1, -1)[strand + 1]
    }
    filterOrig <- filter
    idFilters <- list(
        symbol = c("external_gene_name", "external_gene_id", "hgnc_symbol", "wikigene_name", "dbass3_name"),
        gene = "ensembl_gene_id",
        transcript = "ensembl_transcript_id",
        entrez = "entrezgene"
    )
    staged <- FALSE
    for (i in names(idFilters)) {
        if (!is.null(get(i))) {
            mt <- match(idFilters[[i]], listFilters(biomart)[, 1])
            sfilt <- listFilters(biomart)[mt[!is.na(mt)], 1]
            if (length(sfilt) > 0) {
                staged <- TRUE
                filter[[sfilt[1]]] <- get(i)
                if (!is.null(filter$strand)) {
                    warning(sprintf("Cannot combine %s filter with a strand filter. Strand filtering is ignored.", i))
                    filter$strand <- NULL
                }
                if (!is.null(start) || !is.null(end) || !is.null(chromosome)) {
                    warning(sprintf("Cannot combine %s filter with a range restriction. Ignoring start and end coordinates.", i))
                    start <- end <- chromosome <- NULL
                }
            } else {
                stop(sprintf("Unable to automatically map the %s filter. Manually provide adequate filter list.", i))
            }
        }
    }
    ## We can't have only one in start or end
    if (!is.null(end) && is.null(start)) {
        start <- 1
    } else if (!is.null(start) && is.null(end)) {
        end <- 10e8
    }
    ## filling slots
    .Object@filter <- filter
    .Object@biomart <- biomart
    ## Extending start and end positions to capture genes on the edges. We also need both coordinates set, or none.
    extend <- if (is.null(start) || is.null(end)) 10000 else max(10000, abs(diff(c(end, start))))
    .Object@start <- if (!is.null(start)) max(1, start - extend) else NULL
    .Object@end <- if (!is.null(end)) end + extend else NULL
    displayPars(.Object) <- list(".__featureMap" = featureMap)
    genome(.Object) <- ifelse(is.null(genome), "ANY", genome)
    ## fetching data from Biomart
    range <- if (!is.null(.Object@biomart) && (!is.null(.Object@start) || !is.null(.Object@end) || length(.Object@filter) != 0)) {
        .cacheMartData(.Object, chromosome, staged)
    } else {
        tmp <- GRanges()
        values(tmp) <- DataFrame(feature = "a", gene = "a", exon = "a", transcript = "a", symbol = "a", rank = 0, phase = as.integer(1))[0, ]
        suppressWarnings(genome(tmp) <- unname(genome[1]))
        displayPars(.Object) <- list(".__streamOnly" = TRUE)
        tmp
    }
    .Object@filter <- filterOrig
    if (length(range) == 0) {
        .Object <- setPar(.Object, "size", 0, interactive = FALSE)
    } else {
        chromosome <- if (is.null(chromosome)) seqlevels(range)[1] else chromosome
        rr <- range(range, ignore.strand = TRUE)
        s <- start(rr[seqnames(rr) == chromosome])
        e <- end(rr[seqnames(rr) == chromosome])
        start <- min(s, start)
        end <- max(e, end)
    }
    .Object <- callNextMethod(
        .Object = .Object, range = range, start = start, end = end,
        genome = genome, chromosome = chromosome, strand = strand, ...
    )
    ## We want to warn if searching was performed based on an identifier and now values have been returned
    if ((!is.null(symbol) || !is.null(gene) || !is.null(transcript) || !is.null(entrez)) && length(range) == 0) {
        warning("Search by identifier did not yield any values", call. = FALSE)
    }
    return(.Object)
})

## Constructor ---------------------------------------------------------------

## Constructor. The following arguments are supported:
##    o start, end: numeric vectors of the item start and end coordinates
##    o biomart: a biomaRt object used to to query for gene annotations
##    o genome, chromosome: the reference genome and active chromosome for the track.
##    o strand: character, search for gene models on the plus strand ("+" or 0), the
##       minus strand ("-" or 1) or both strands ("+-" or "-+" or 2)
##    o stacking: character controlling the stacking of overlapping items. One in 'hide',
##       'dense', 'squish', 'pack' or 'full'.
##    o filters: list of additional filters for the biomaRt query, where the item names
##          are the filter identifiers and the item values are the filter values.
##    o name: the name of the track. This will be used for the title panel.
## All additional items in ... are being treated as DisplayParameters

#' @describeIn BiomartGeneRegionTrack-class Constructor function for
#' `BiomartGeneRegionTrack-class`.
#' @export
BiomartGeneRegionTrack <- function(start = NULL, end = NULL, biomart, chromosome = NULL, strand, genome = NULL,
                                   stacking = "squish", filters = list(), featureMap = NULL, name = "BiomartGeneRegionTrack",
                                   symbol = NULL, gene = NULL, entrez = NULL, transcript = NULL, ...) {
    ## Some default checking
    .missingToNull(c("genome"))
    if (missing(strand)) {
        strand <- "*"
    }
    if ((!is.null(start) || !is.null(end)) && is.null(chromosome)) {
        stop("Also need to specify a chromsome when initializing a BiomartGeneRegionTrack with start or end coordinates")
    }
    if (!is.null(chromosome)) {
        chromosome <- .chrName(chromosome)[1]
    }
    if (missing(biomart)) {
        if (is.null(genome)) {
            stop("Need either a valid Mart connection object as 'biomart' argument or a UCSC genome identifier as the 'genome' argument.")
        }
        biomart <- .genome2Dataset(genome)
    } else if (is.null(genome)) {
        genome <- biomart@dataset
    }
    new("BiomartGeneRegionTrack",
        start = start, end = end, chromosome = chromosome, strand = strand,
        biomart = biomart, name = name, genome = genome, stacking = stacking, filter = filters, featureMap = featureMap,
        symbol = symbol, gene = gene, transcript = transcript, entrez = entrez, ...
    )
}


## This filters for genes with refseq IDs only
.refseqFilter <- list("with_refseq_dna" = TRUE)

## General accessors ---------------------------------------------------------
## Annotation Accessors ------------------------------------------------------
## Stacking ------------------------------------------------------------------
## Consolidate ---------------------------------------------------------------
## Collapse  -----------------------------------------------------------------
## Subset --------------------------------------------------------------------

#' @describeIn BiomartGeneRegionTrack-class subset a `BiomartGeneRegionTrack`
#' by coordinates and sort if necessary.
#' @export
setMethod("subset", signature(x = "BiomartGeneRegionTrack"), function(x, from, to, chromosome, use.defaults = TRUE, ...) {
    granges <- unlist(range(split(ranges(x), group(x))))
    ranges <- if (use.defaults) {
        .defaultRange(x, from = from, to = to)
    } else {
        c(
            from = ifelse(is.null(from), min(start(granges)) - 1, from),
            to = ifelse(is.null(to), max(end(granges)) + 1, to)
        )
    }
    if (ranges["from"] < x@start || ranges["to"] > x@end) {
        x@start <- ranges["from"] - 10000
        x@end <- ranges["to"] + 10000
        ranges(x) <- .cacheMartData(x, .chrName(chromosome))
    }
    return(callNextMethod(x = x, from = ranges["from"], to = ranges["to"], use.defaults = FALSE, ...))
})

## Position ------------------------------------------------------------------
## DrawGrid ------------------------------------------------------------------
## DrawAxis ------------------------------------------------------------------
## DrawGD --------------------------------------------------------------------
## SetAs ---------------------------------------------------------------------
## Show ----------------------------------------------------------------------
