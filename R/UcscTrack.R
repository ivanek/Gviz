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

#' @importFrom rtracklayer ucscGenomes browserSession
#' @importMethodsFrom rtracklayer chrom close getTable "tableName<-" track
#' ucscTableQuery trackNames tableNames import import.gff import.gff1
#' import.gff2 import.gff3 import.2bit import.bed15 import.bw import.ucsc
#' import.bed import.bedGraph import.chain import.wig seqinfo
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
    genomes <- .doCache("availableGenomes", expression(rtracklayer::ucscGenomes()), env)
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
            myUcscUrl <- getOption("Gviz.ucscUrl")
            tmp <- if (is.null(myUcscUrl)) browserSession() else browserSession(url = myUcscUrl)
            genome(tmp) <- genome
            tmp
        }), env, cenv
    )
    availTracks <- .doCache(tracksToken, expression(trackNames(ucscTableQuery(session))), env, cenv)
    track <- match.arg(track, sort(c(availTracks, names(availTracks))))
    if (!is.na(availTracks[track])) {
        track <- names(availTracks[track])
    }
    availTables <- .doCache(tablesToken, expression({
        query <- ucscTableQuery(session, track)
        sort(tableNames(query))
    }), env, cenv)
    chrInfo <- seqlengths(session)
    return(list(
        session = session, availTracks = availTracks, availTables = availTables, track = track,
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
        bands <- .doCache(genomesToken, expression({
            if (!genome %in% as.character(genomes[, "db"])) {
                stop("'", genome, "' is not a valid UCSC genome.")
            }
            sessionToken <- paste("session", genome, sep = "_")
            session <- .doCache(
                sessionToken,
                expression({
                    myUcscUrl <- getOption("Gviz.ucscUrl")
                    tmp <- if (is.null(myUcscUrl)) browserSession() else browserSession(url = myUcscUrl)
                    genome(tmp) <- genome
                    tmp
                }), env, cenv
            )
            query <- tryCatch(ucscTableQuery(session, table = "cytoBandIdeo"), error = function(e) {
                warning(
                    "There doesn't seem to be any cytoband data available for genome '", genome,
                    "' at UCSC or the service is temporarily down. Trying to fetch the chromosome length data."
                )
                tryCatch(ucscTableQuery(session, table = "chromInfo"), error = function(e) {
                    stop(
                        "There doesn't seem to be any chromosome length data available for genome '", genome,
                        "' at UCSC or the service is temporarily down."
                    )
                })
            })
            out <- getTable(query)
            if (all(c("chrom", "size") %in% colnames(out))) {
                out <- data.frame(chrom = out$chrom, chromStart = 0, chromEnd = out$size, name = "", gieStain = "gneg", stringsAsFactors = FALSE)
            }
            out
        }), env, cenv)
    }

    return(list(availableGenomes = genomes, bands = bands))
}
.cacheMartData <- function(bmtrack, chromosome = NULL, staged = FALSE) {
    uid <- .bmGuid(bmtrack)
    req <- if (!is.null(bmtrack@start) && !is.null(bmtrack@end)) GRanges(seqnames = chromosome[1], IRanges(start = bmtrack@start, bmtrack@end)) else NULL
    if (is.null(chromosome) || is.null(.martCache[[uid]])) {
        data <- .fetchBMData(bmtrack, chromosome, staged)
        if (!is.null(req)) {
            .martCache[[uid]] <- list(data = data, ranges = req)
        } else {
            req <- range(data)
        }
        if (length(data) && .dpOrDefault(bmtrack, "verbose", FALSE)) {
            message("Loaded data from Biomart for region ", paste(sprintf("%s:%i-%i(%s)", seqnames(req), start(req), end(req), strand(req)), collapse = " and "))
        }
    } else {
        rr <- .martCache[[uid]][["ranges"]]
        dd <- .martCache[[uid]][["data"]]
        if (!is.null(req) && suppressWarnings(req %within% rr)) {
            genes <- unique(subsetByOverlaps(dd, req)$gene)
            data <- dd[seqnames(dd) == chromosome[1] & dd$gene %in% genes]
            if (.dpOrDefault(bmtrack, "verbose", FALSE)) {
                message(sprintf("Retrieved data from cache for region %s:%i-%i(%s)", chromosome, start(req), end(req), strand(req)))
            }
        } else {
            data <- .fetchBMData(bmtrack, chromosome, staged)
            if (is.null(req)) {
                req <- range(data)
            }
            .martCache[[uid]][["data"]] <- suppressWarnings(c(.martCache[[uid]][["data"]], data[!(seqnames(data) == chromosome[1] & data$gene %in% dd$gene)]))
            .martCache[[uid]][["ranges"]] <- suppressWarnings(union(rr, req))
            if (length(req) && .dpOrDefault(bmtrack, "verbose", FALSE)) {
                message("Loaded data from Biomart for region ", paste(sprintf("%s:%i-%i(%s)", seqnames(req), start(req), end(req), strand(req)), collapse = " and "))
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
#' UCSC data sources.
#'
#'
#' The UCSC data base provides a wealth of annotation information. This
#' function can be used to access UCSC, to retrieve the data available there
#' and to return it as an annotation track object amenable to plotting with
#' \code{\link{plotTracks}}.
#'
#' \code{clearSessionCache} is can be called to remove all cached items from
#' the session which are generated when connecting with the UCSC data base.
#'
#'
#' The data stored at the UCSC data bases can be of different formats: gene or
#' transcript model data, simple annotation features like CpG Island locations
#' or SNPs, or numeric data like conservation or mapability. This function
#' presents a unified API to download all kinds of data and to map them back to
#' one of the annotation track objects defined in this package. The type of
#' object to hold the data has to be given in the \code{trackType} argument,
#' and subsequently the function passes all data on to the respective object
#' constructor. All additional named arguments are considered to be relevant
#' for the constructor of choice, and single character scalars are replaced by
#' the respective data columns in the downloaded UCSC tables if available. For
#' instance, assuming the table for track 'foo' contains the columns 'id',
#' 'type', 'fromLoc' and 'toLoc', giving the feature identifier, type, start
#' end end location. In order to create an \code{\linkS4class{AnnotationTrack}}
#' object from that data, we have to pass the additional named arguments
#' \code{id="id"}, \code{feature="type"}, \code{start="fromLoc"} and
#' codeend="toLoc" to the \code{UcscTrack} function. The complete function call
#' could look like this:
#'
#' \code{UcscTrack(track="foo", genome="mm9", chromosome=3, from=1000,
#' to=10000, trackType="AnnotationTrack", id="id", feature="type",
#' start="from", end="to")}
#'
#' To reduce the bandwidth, some caching of the UCSC connection takes place. In
#' order to remove these cached session items, call \code{clearSessionCache}.
#'
#' The \code{Gviz.ucscUrl} option controls which URL is being used to connect
#' to UCSC. For instance, one could switch to the European UCSC mirror by
#' calling \code{options(Gviz.ucscUrl="http://genome-euro.ucsc.edu/cgi-bin/"}.
#'
#' @aliases UcscTrack clearSessionCache
#' @param track Character, the name of the track to fetch from UCSC. To find
#' out about available tracks please consult the online table browser at
#' \url{http://genome.ucsc.edu/cgi-bin/hgTables?command=start}.
#' @param table Character, the name of the table to fetch from UCSC, or
#' \code{NULL}, in which case the default selection of tables is used. To find
#' out about available tables for a given track please consult the online table
#' browser at \url{http://genome.ucsc.edu/cgi-bin/hgTables?command=start}.
#' @param trackType Character, one in \code{c("AnnotationTrack",
#' "GeneRegionTrack", "DataTrack", "GenomeAxisTrack")}. The function will try
#' to coerce the downloaded data in an object of this class. See below for
#' details.
#' @param genome Character, a valid USCS genome identifier for which to fetch
#' the data.
#' @param chromosome Character, a valid USCS character identifier for which to
#' fetch the data.
#' @param name Character, the name to use for the resulting track object.
#' @param from,to A range of genomic locations for which to fetch data.
#' @param \dots All additional named arguments are expected to be either
#' display parameters for the resulting objects, or character scalars of column
#' names in the downloaded UCSC data tables that are matched by name to
#' available arguments in the respective constructor functions as defined by
#' the \code{trackType} argument. See Details section for more information.
#' @return
#'
#' An annotation track object as determined by \code{trackType}.
#' @author Florian Hahne
#' @seealso
#'
#' \code{\linkS4class{AnnotationTrack}}
#'
#' \code{\linkS4class{DataTrack}}
#'
#' \code{\linkS4class{GeneRegionTrack}}
#'
#' \code{\linkS4class{GenomeAxisTrack}}
#'
#' \code{\link{plotTracks}}
#' @examples
#' \dontrun{
#'
#' ## Create UcscTrack for Known Genes from mm9 genome
#' from <- 65921878
#' to <- 65980988
#' knownGenes <- UcscTrack(
#'     genome = "mm9", chromosome = "chrX", track = "knownGene",
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
UcscTrack <- function(track, table = NULL,
                      trackType = c(
                          "AnnotationTrack", "GeneRegionTrack",
                          "DataTrack", "GenomeAxisTrack"
                      ),
                      genome, chromosome, name = NULL, from, to, ...) {
    trackType <- match.arg(trackType)
    if (missing(genome) || !isSingleString(genome)) stop("Need to specify genome for creating a UcscTrack")
    if (missing(chromosome)) stop("Need to specify chromosome for creating a UcscTrack")
    chromosome <- .chrName(chromosome)[1]
    sessionInfo <- .cacheTracks(genome = genome, chromosome = chromosome, track = track, env = .ucscCache)
    if (missing(from)) {
        from <- 1
    }
    if (missing(to)) {
        to <- sessionInfo$chrInfo[chromosome]
    }
    gr <- GRanges(ranges = IRanges(start = from, end = to), seqnames = chromosome)
    suppressWarnings(genome(gr) <- unname(genome))[1]
    query <- ucscTableQuery(sessionInfo$session, sessionInfo$track, gr)
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
        stop("Error fetching data from UCSC")
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
    trackObject <- do.call(trackType, args = c(list(chromosome = chromosome, genome = genome, name = name), args))
    return(trackObject)
}
