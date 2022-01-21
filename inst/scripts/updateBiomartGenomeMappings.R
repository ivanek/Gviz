## This script will try to identify genome versions in all available
## ENSEBML archives as well as the current release, and match these
## versions back to UCSC genome identifiers.

library(biomaRt)

parseArchiveDate <- function(x) {
    datePart <- sub("\\.archive.*", "", sub("^https://", "", x))
    return(as.Date(paste(datePart, "1", sep = ""), format = "%b%Y%d"))
}


getBmFromDate <- function(date) {
    archive <- sprintf("https://%s%s.archive.ensembl.org", format(date, "%b"), format(date, "%Y"))
    res <- try(capture.output(suppressMessages(bm <- useMart(
        host = archive,
        biomart = "ENSEMBL_MART_ENSEMBL"
    ))),
    silent = TRUE
    )
    if (is(res, "try-error")) {
        bm <- NA
    }
    return(bm)
}

year <- function(x) {
    if (missing(x)) {
        x <- Sys.Date()
    }
    return(as.integer(format(x, "%Y")))
}

month <- function(x) {
    if (missing(x)) {
        x <- Sys.Date()
    }
    return(as.integer(format(x, "%m")))
}

day <- function(x) {
    if (missing(x)) {
        x <- Sys.Date()
    }
    return(as.integer(format(x, "%d")))
}

knownEnsemblArchives <- c(
    `Ensembl 105: Dec 2021` = "https://Dec2021.archive.ensembl.org",
    `Ensembl 104: May 2021` = "https://May2021.archive.ensembl.org",
    `Ensembl 103: Feb 2021` = "https://Feb2021.archive.ensembl.org",
    `Ensembl 102: Nov 2020` = "https://Nov2020.archive.ensembl.org",
    `Ensembl 101: Aug 2020` = "https://Aug2020.archive.ensembl.org",
    `Ensembl 100: Apr 2020` = "https://Apr2020.archive.ensembl.org",
    `Ensembl 99: Jan 2020` = "https://Jan2020.archive.ensembl.org",
    `Ensembl 98: Sep 2019` = "https://Sep2019.archive.ensembl.org",
    `Ensembl 97: Jul 2019` = "https://Jul2019.archive.ensembl.org",
    `Ensembl 96: Apr 2019` = "https://Apr2019.archive.ensembl.org",
    `Ensembl 95: Jan 2019` = "https://Jan2019.archive.ensembl.org",
    `Ensembl 94: Oct 2018` = "https://Oct2018.archive.ensembl.org",
    `Ensembl 93: Jul 2018` = "https://Jul2018.archive.ensembl.org",
    `Ensembl 92: Apr 2018` = "https://Apr2018.archive.ensembl.org",
    `Ensembl 91: Dec 2017` = "https://Dec2017.archive.ensembl.org",
    `Ensembl 90: Aug 2017` = "https://Aug2017.archive.ensembl.org",
    `Ensembl 89: May 2017` = "https://May2017.archive.ensembl.org",
    `Ensembl 88: Mar 2017` = "https://Mar2017.archive.ensembl.org",
    `Ensembl 87: Dec 2016` = "https://Dec2016.archive.ensembl.org",
    `Ensembl 86: Oct 2016` = "https://Oct2016.archive.ensembl.org",
    `Ensembl 80: May 2015` = "https://May2015.archive.ensembl.org",
    `Ensembl 77: Oct 2014` = "https://Oct2014.archive.ensembl.org",
    `Ensembl 75: Feb 2014` = "https://Feb2014.archive.ensembl.org",
    `Ensembl 54: May 2009` = "https://May2009.archive.ensembl.org"
)

parseOneArchive <- function(date) {
    bm <- getBmFromDate(date)
    datasets <- listDatasets(bm)
    return(datasets)
}

parseEnsemblArchives <- function(known) {
    archiveDates <- sort(parseArchiveDate(known))
    ## We may have to extend our search date range in case that the known archives vector is
    ## out of date and later archives exist
    latest <- max(archiveDates)
    daysBack <- as.integer(Sys.Date() - latest)
    if (daysBack > day()) {
        daysSince <- Sys.Date() - (1:daysBack)
        archiveDates <- sort(unique(c(
            archiveDates,
            as.Date(paste(unique(format(daysSince, "%b%Y")), "1", sep = ""),
                format = "%b%Y%d"
            )
        )))
    }
}
