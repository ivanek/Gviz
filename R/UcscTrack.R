## UcscTrack -----------------------------------------------------------------
##
## Strictly speaking this is not a class, but rather some sort of meta-constructor for several of the previously
## defined track types directly from UCSC data. It will fetch online data from a particular track (or a sub-table of
## a track) and feed it to one of the original constructors with user-provided argument mappings.
## ---------------------------------------------------------------------------
## Constructor. The following arguments are supported:
##    o track: a character of one of the available UCSC tracks
##    o table: a character of one one of the sub-tables of the track, or NULL to fetch all
##    o trackType: a character giving the name of the constructor to pass the data to, one in
##       c("AnnotationTrack", "GeneRegionTrack", "DataTrack", "GenomeAxisTrack")
##    o genome, chromosome: the reference genome and active chromosome for the track.
##    o from: the starting and end coordinates of the track data
##    o name: the name of the track. This will be used for the title panel.
## All additional items in ... are being treated as DisplayParameters

## A simple caching mechanism for UCSC session information. The overhead for establishing a connection to UCSC is
## quite significant and we can shave off 5 to 10 seconds here by caching sessions and associated information
## for a particular genome and chromosome.
.ucscCache <- new.env()
.ensemblCache <- new.env()
.martCache <- new.env()

## rtracklayer's tableNames() (the REST-API-era replacement for the now
## defunct trackNames()) does not always return a simple named character
## vector: it can come back as a list, with NULL entries for
## access-protected tracks and vector-valued entries for tracks that have
## more than one sub-table. sort()/match.arg() need an atomic vector, so we
## flatten it here, dropping protected (NULL) entries and replicating the
## track name across each of its sub-tables.
.flattenTableNames <- function(tn) {
    if (!is.list(tn)) {
        return(tn)
    }
    tn <- tn[!vapply(tn, is.null, logical(1))]
    if (!length(tn)) {
        return(character(0))
    }
    trackIds <- rep(names(tn), lengths(tn))
    tables <- unlist(tn, use.names = FALSE)
    names(tables) <- trackIds
    tables
}

## Some UCSC tracks (notably "knownGene" on hg38 and other recently-updated
## assemblies) are no longer served as the classic genePred table with
## absolute-coordinate, comma-separated "exonStarts"/"exonEnds" columns.
## They now come back in the bigGenePred/BED12+ schema instead, where exon
## (block) coordinates are given as "chromStarts" (offsets *relative* to
## "chromStart") together with "blockSizes", following the ordinary BED12
## convention. If we don't bridge this, rstarts/rends arguments such as the
## commonly-documented rstarts="exonStarts", rends="exonEnds" silently fail
## to resolve to real columns (they get passed through as literal strings)
## and blow up deep inside GeneRegionTrack's range-building code with a
## confusing "Number of elements ... is invalid" error.
##
## This synthesizes legacy-style absolute exonStarts/exonEnds (and a few
## other commonly-used legacy aliases) from the new columns whenever they're
## present and the legacy ones are missing, so existing calls keep working
## unchanged. It's a no-op when the table already looks like a classic
## genePred (or doesn't have enough bigGenePred columns to convert).
.bigGenePredToGenePredCompat <- function(tableDat) {
    bigGenePredCols <- c("chromStart", "chromEnd", "blockSizes", "chromStarts")
    if (
        !is.data.frame(tableDat) ||
            !all(bigGenePredCols %in% colnames(tableDat))
    ) {
        return(tableDat)
    }
    if (all(c("exonStarts", "exonEnds") %in% colnames(tableDat))) {
        return(tableDat)
    }
    parseCsv <- function(x) as.integer(strsplit(sub(",+$", "", x), ",")[[1]])
    toCsv <- function(x) paste0(paste(x, collapse = ","), ",")
    exonCoords <- Map(
        function(chromStart, blockStarts, blockSizes) {
            starts <- chromStart + parseCsv(blockStarts)
            ends <- starts + parseCsv(blockSizes)
            list(starts = toCsv(starts), ends = toCsv(ends))
        },
        tableDat$chromStart,
        tableDat$chromStarts,
        tableDat$blockSizes
    )
    tableDat$exonStarts <- vapply(exonCoords, `[[`, character(1), "starts")
    tableDat$exonEnds <- vapply(exonCoords, `[[`, character(1), "ends")
    if (
        !"exonCount" %in% colnames(tableDat) &&
            "blockCount" %in% colnames(tableDat)
    ) {
        tableDat$exonCount <- tableDat$blockCount
    }
    if (!"txStart" %in% colnames(tableDat)) {
        tableDat$txStart <- tableDat$chromStart
    }
    if (!"txEnd" %in% colnames(tableDat)) {
        tableDat$txEnd <- tableDat$chromEnd
    }
    if (
        !"cdsStart" %in% colnames(tableDat) &&
            "thickStart" %in% colnames(tableDat)
    ) {
        tableDat$cdsStart <- tableDat$thickStart
    }
    if (
        !"cdsEnd" %in% colnames(tableDat) && "thickEnd" %in% colnames(tableDat)
    ) {
        tableDat$cdsEnd <- tableDat$thickEnd
    }
    tableDat
}

## Default UCSC base URL. UCSC's table-browser backend now sits behind
## api.genome.ucsc.edu / the hubApi REST interface. Some rtracklayer
## versions mishandle the plain "http://" redirect to "https://" for that
## backend (see https://github.com/lawremi/rtracklayer/issues/148), so we
## default to https explicitly unless the user overrides Gviz.ucscUrl.
.gvizUcscUrl <- function() {
    myUcscUrl <- getOption("Gviz.ucscUrl")
    if (is.null(myUcscUrl)) "https://genome.ucsc.edu/cgi-bin/" else myUcscUrl
}

## Newer rtracklayer releases query UCSC via the REST-based hubApi, which
## flattened the old track/table distinction: track information now has to
## be requested through the 'table' argument of ucscTableQuery() instead of
## the (now deprecated/unreliable) positional 'track' argument. Older
## rtracklayer releases still expect 'track'. This helper tries the modern
## calling convention first and falls back to the legacy one so that
## UcscTrack() keeps working across rtracklayer versions.
.ucscTableQueryCompat <- function(session, track, range = NULL) {
    modernArgs <- c(
        list(session, table = track),
        if (!is.null(range)) list(range = range)
    )
    res <- tryCatch(do.call(ucscTableQuery, modernArgs), error = function(e) e)
    if (inherits(res, "error")) {
        legacyArgs <- c(
            list(session, track),
            if (!is.null(range)) list(range = range)
        )
        res <- tryCatch(
            do.call(ucscTableQuery, legacyArgs),
            error = function(e) e
        )
    }
    if (inherits(res, "error")) {
        stop(
            "Unable to query UCSC track/table '",
            track,
            "'. This may be caused by an incompatible ",
            "rtracklayer version or a change in the UCSC REST API. Cause: ",
            conditionMessage(res)
        )
    }
    res
}

#' @importFrom rtracklayer ucscGenomes browserSession
#' @importMethodsFrom rtracklayer chrom close getTable "tableName<-" track
#' @importMethodsFrom rtracklayer ucscTableQuery trackNames tableNames import
#' @importMethodsFrom rtracklayer import.gff import.gff1 import.gff2
#' @importMethodsFrom rtracklayer import.gff3 import.2bit import.bed15
#' @importMethodsFrom rtracklayer import.bw import.ucsc import.bed
#' @importMethodsFrom rtracklayer import.bedGraph import.chain import.wig
#' @importMethodsFrom rtracklayer seqinfo
#'
#' @export
.doCache <- function(token, expression, env, callEnv = environment()) {
    if (!token %in% base::ls(env)) {
        res <- eval(expression, envir = callEnv)
        assign(x = token, value = res, envir = env)
        res
    } else {
        env[[token]]
    }
}
.cacheTracks <- function(genome, chromosome, track, env = .ucscCache) {
    genomes <- .doCache(
        "availableGenomes",
        expression(rtracklayer::ucscGenomes()),
        env
    )
    if (!genome %in% as.character(genomes[, "db"])) {
        stop("'", genome, "' is not a valid UCSC genome.")
    }
    sessionToken <- paste("session", genome, sep = "_")
    tracksToken <- paste("tracks", genome, sep = "_")
    tablesToken <- paste("tables", track, genome, sep = "_")
    cenv <- environment()
    session <- .doCache(
        sessionToken,
        expression({
            tmp <- browserSession(url = .gvizUcscUrl())
            genome(tmp) <- genome
            tmp
        }),
        env,
        cenv
    )
    availTracks <- .doCache(
        tracksToken,
        expression(.flattenTableNames(tableNames(ucscTableQuery(session)))),
        env,
        cenv
    )
    track <- match.arg(track, sort(c(availTracks, names(availTracks))))
    if (!is.na(availTracks[track])) {
        track <- names(availTracks[track])
    }
    availTables <- .doCache(
        tablesToken,
        expression({
            query <- .ucscTableQueryCompat(session, track)
            sort(.flattenTableNames(tableNames(query)))
        }),
        env,
        cenv
    )
    chrInfo <- seqlengths(session)
    return(list(
        session = session,
        availTracks = availTracks,
        availTables = availTables,
        track = track,
        chrInfo = chrInfo
    ))
}
.cacheGenomes <- function(genome = NULL, env = .ucscCache) {
    availToken <- "availableGenomes"
    genomesToken <- paste("genomeBands", genome, sep = "_")
    genomes <- .doCache(availToken, expression(rtracklayer::ucscGenomes()), env)
    bands <- NULL
    if (!is.null(genome)) {
        cenv <- environment()
        bands <- .doCache(
            genomesToken,
            expression({
                if (!genome %in% as.character(genomes[, "db"])) {
                    stop("'", genome, "' is not a valid UCSC genome.")
                }
                sessionToken <- paste("session", genome, sep = "_")
                session <- .doCache(
                    sessionToken,
                    expression({
                        tmp <- browserSession(url = .gvizUcscUrl())
                        genome(tmp) <- genome
                        tmp
                    }),
                    env,
                    cenv
                )
                query <- tryCatch(
                    ucscTableQuery(session, table = "cytoBandIdeo"),
                    error = function(e) {
                        warning(
                            "There doesn't seem to be any cytoband data available for genome '",
                            genome,
                            "' at UCSC or the service is temporarily down. Trying to fetch the chromosome length data."
                        )
                        tryCatch(
                            ucscTableQuery(session, table = "chromInfo"),
                            error = function(e) {
                                stop(
                                    "There doesn't seem to be any chromosome length data available for genome '",
                                    genome,
                                    "' at UCSC or the service is temporarily down."
                                )
                            }
                        )
                    }
                )
                out <- getTable(query)
                if (all(c("chrom", "size") %in% colnames(out))) {
                    out <- data.frame(
                        chrom = out$chrom,
                        chromStart = 0,
                        chromEnd = out$size,
                        name = "",
                        gieStain = "gneg",
                        stringsAsFactors = FALSE
                    )
                }
                out
            }),
            env,
            cenv
        )
    }

    return(list(availableGenomes = genomes, bands = bands))
}
.cacheMartData <- function(bmtrack, chromosome = NULL, staged = FALSE) {
    uid <- .bmGuid(bmtrack)
    req <- if (!is.null(bmtrack@start) && !is.null(bmtrack@end)) {
        GRanges(
            seqnames = chromosome[1],
            IRanges(start = bmtrack@start, bmtrack@end)
        )
    } else {
        NULL
    }
    if (is.null(chromosome) || is.null(.martCache[[uid]])) {
        data <- .fetchBMData(bmtrack, chromosome, staged)
        if (!is.null(req)) {
            .martCache[[uid]] <- list(data = data, ranges = req)
        } else {
            req <- range(data)
        }
        if (length(data) && .dpOrDefault(bmtrack, "verbose", FALSE)) {
            message(
                "Loaded data from Biomart for region ",
                paste(
                    sprintf(
                        "%s:%i-%i(%s)",
                        seqnames(req),
                        start(req),
                        end(req),
                        strand(req)
                    ),
                    collapse = " and "
                )
            )
        }
    } else {
        rr <- .martCache[[uid]][["ranges"]]
        dd <- .martCache[[uid]][["data"]]
        if (!is.null(req) && suppressWarnings(req %within% rr)) {
            genes <- unique(subsetByOverlaps(dd, req)$gene)
            data <- dd[seqnames(dd) == chromosome[1] & dd$gene %in% genes]
            if (.dpOrDefault(bmtrack, "verbose", FALSE)) {
                message(sprintf(
                    "Retrieved data from cache for region %s:%i-%i(%s)",
                    chromosome,
                    start(req),
                    end(req),
                    strand(req)
                ))
            }
        } else {
            data <- .fetchBMData(bmtrack, chromosome, staged)
            if (is.null(req)) {
                req <- range(data)
            }
            .martCache[[uid]][["data"]] <- suppressWarnings(c(
                .martCache[[uid]][["data"]],
                data[
                    !(seqnames(data) == chromosome[1] & data$gene %in% dd$gene)
                ]
            ))
            .martCache[[uid]][["ranges"]] <- suppressWarnings(union(rr, req))
            if (length(req) && .dpOrDefault(bmtrack, "verbose", FALSE)) {
                message(
                    "Loaded data from Biomart for region ",
                    paste(
                        sprintf(
                            "%s:%i-%i(%s)",
                            seqnames(req),
                            start(req),
                            end(req),
                            strand(req)
                        ),
                        collapse = " and "
                    )
                )
            }
        }
    }
    return(data)
}

## empty the session cache
#' @export
clearSessionCache <- function() {
    assignInNamespace(".ucscCache", new.env(), ns = "Gviz")
    assignInNamespace(".ensemblCache", new.env(), ns = "Gviz")
    assignInNamespace(".martCache", new.env(), ns = "Gviz")
}


## Constructor

#' Meta-constructor for Gviz tracks fetched directly from the various
#' UCSC data sources
#'
#'
#' The UCSC data base provides a wealth of annotation information. This
#' function can be used to access UCSC, to retrieve the data available there
#' and to return it as an annotation track object amenable to plotting with
#' [`plotTracks`].
#'
#' `clearSessionCache` can be called to remove all cached items from
#' the session which are generated when connecting with the UCSC data base.
#'
#'
#' The data stored at the UCSC data bases can be of different formats: gene or
#' transcript model data, simple annotation features like CpG Island locations
#' or SNPs, or numeric data like conservation or mappability. This function
#' presents a unified API to download all kinds of data and to map them back to
#' one of the annotation track objects defined in this package. The type of
#' object to hold the data has to be given in the `trackType` argument, and
#' subsequently the function passes all data on to the respective object
#' constructor. All additional named arguments are considered to be relevant for
#' the constructor of choice, and single character scalars are replaced by the
#' respective data columns in the downloaded UCSC tables if available. For
#' instance, assuming the table for track 'foo' contains the columns 'id',
#' 'type', 'fromLoc' and 'toLoc', giving the feature identifier, type, start end
#' end location. In order to create an
#' [`AnnotationTrack`][AnnotationTrack-class] object from that data, we have to
#' pass the additional named arguments `id="id"`, `feature="type"`,
#' `start="fromLoc"` and `end="toLoc"` to the `UcscTrack` function. The complete
#' function call could look like this:
#'
#' `UcscTrack(track="foo", genome="mm39", chromosome=3, from=1000,
#' to=10000, trackType="AnnotationTrack", id="id", feature="type",
#' start="from", end="to")`
#'
#' To reduce the bandwidth, some caching of the UCSC connection takes place. In
#' order to remove these cached session items, call `clearSessionCache`.
#'
#' The `Gviz.ucscUrl` option controls which URL is being used to connect
#' to UCSC. For instance, one could switch to the European UCSC mirror by
#' calling `options(Gviz.ucscUrl = "http://genome-euro.ucsc.edu/cgi-bin/")`.
#'
#' @aliases UcscTrack clearSessionCache
#' @param track Character, the name of the track to fetch from UCSC. To find
#' out about available tracks please consult the online table browser at
#' <http://genome.ucsc.edu/cgi-bin/hgTables?command=start>.
#' @param table Character, the name of the table to fetch from UCSC, or
#' `NULL`, in which case the default selection of tables is used. To find
#' out about available tables for a given track please consult the online table
#' browser at <http://genome.ucsc.edu/cgi-bin/hgTables?command=start>.
#' @param trackType Character, one of `c("AnnotationTrack",
#' "GeneRegionTrack", "DataTrack", "GenomeAxisTrack")`. The function will try
#' to coerce the downloaded data in an object of this class. See below for
#' details.
#' @param genome Character, a valid UCSC genome identifier for which to fetch
#' the data.
#' @param chromosome Character, a valid UCSC character identifier for which to
#' fetch the data.
#' @param name Character, the name to use for the resulting track object.
#' @param from,to A range of genomic locations for which to fetch data.
#' @param ... All additional named arguments are expected to be either
#' display parameters for the resulting objects, or character scalars of column
#' names in the downloaded UCSC data tables that are matched by name to
#' available arguments in the respective constructor functions as defined by
#' the `trackType` argument. See Details section for more information.
#' @return
#'
#' An annotation track object as determined by `trackType`.
#' @author Florian Hahne
#' @seealso
#'
#' [AnnotationTrack-class]
#'
#' [DataTrack-class]
#'
#' [GeneRegionTrack-class]
#'
#' [GenomeAxisTrack-class]
#'
#' [plotTracks]
#' @examples
#' \donttest{
#'
#' ## Create UcscTrack for Known Genes from mm39 genome
#' from <- 65921878
#' to <- 65980988
#' knownGenes <- UcscTrack(
#'     genome = "mm39", chromosome = "chrX", track = "knownGene",
#'     from = from, to = to, trackType = "GeneRegionTrack",
#'     rstarts = "exonStarts", rends = "exonEnds", gene = "name",
#'     symbol = "name", transcript = "name", strand = "strand",
#'     fill = "#8282d2", name = "UCSC Genes"
#' )
#' }
#'
#' ## if the UCSC is not accessible load prepared object
#' data(ucscItems)
#'
#' ## knownGenes is essentially GeneRegionTrack
#' knownGenes
#'
#' ## plotting
#' plotTracks(knownGenes, chromosome = "chrX", from = 65920688, to = 65960068)
#' @export
UcscTrack <- function(
  track,
  table = NULL,
  trackType = c(
      "AnnotationTrack",
      "GeneRegionTrack",
      "DataTrack",
      "GenomeAxisTrack"
  ),
  genome,
  chromosome,
  name = NULL,
  from,
  to,
  ...
) {
    trackType <- match.arg(trackType)
    if (missing(genome) || !isSingleString(genome)) {
        stop("Need to specify genome for creating a UcscTrack")
    }
    if (missing(chromosome)) {
        stop("Need to specify chromosome for creating a UcscTrack")
    }
    chromosome <- .chrName(chromosome)[1]
    sessionInfo <- .cacheTracks(
        genome = genome,
        chromosome = chromosome,
        track = track,
        env = .ucscCache
    )
    if (missing(from)) {
        from <- 1
    }
    if (missing(to)) {
        to <- sessionInfo$chrInfo[chromosome]
    }
    gr <- GRanges(
        ranges = IRanges(start = from, end = to),
        seqnames = chromosome
    )
    suppressWarnings(genome(gr) <- unname(genome))[1]
    query <- .ucscTableQueryCompat(
        sessionInfo$session,
        sessionInfo$track,
        range = gr
    )
    if (!is.null(table)) {
        table <- match.arg(table, sessionInfo$availTables)
        tableName(query) <- table
    }
    if (is.null(name)) {
        name <- if (is.null(table)) track else paste(sessionInfo$track, table)
    }
    tableDat <- if (trackType == "DataTrack") {
        tmp <- try(track(query), silent = TRUE)
        if (is(tmp, "try-error")) {
            warning(tmp)
            data.frame()
        } else {
            as.data.frame(tmp)
        }
    } else {
        tmp <- try(getTable(query), silent = TRUE)
        if (is(tmp, "try-error")) {
            warning(tmp)
            data.frame()
        } else {
            tmp
        }
    }
    if (is(tmp, "try-error") && nrow(tableDat) == 0) {
        stop("Unable to fetch data from UCSC")
    }
    tableDat <- .bigGenePredToGenePredCompat(tableDat)
    if (trackType == "GeneRegionTrack") {
        dots <- list(...)
        for (colArg in c("rstarts", "rends")) {
            val <- dots[[colArg]]
            if (
                is.character(val) &&
                    length(val) == 1 &&
                    !val %in% colnames(tableDat)
            ) {
                stop(
                    "Column '",
                    val,
                    "' (",
                    colArg,
                    ") was not found in the data fetched from UCSC ",
                    "for track '",
                    track,
                    "'. This can happen when UCSC changes a track's table schema. ",
                    "Available columns are: ",
                    paste(colnames(tableDat), collapse = ", ")
                )
            }
        }
    }
    args <- lapply(list(...), function(x) {
        if (is.character(x) && length(x) == 1) {
            if (!x %in% colnames(tableDat)) x else unlist(tableDat[, x])
        } else {
            x
        }
    })
    if (trackType == "GeneRegionTrack") {
        args$start <- from
        args$end <- to
    }
    args <- lapply(args, function(x) if (!length(x)) NULL else x)
    trackObject <- do.call(
        trackType,
        args = c(
            list(chromosome = chromosome, genome = genome, name = name),
            args
        )
    )
    return(trackObject)
}
