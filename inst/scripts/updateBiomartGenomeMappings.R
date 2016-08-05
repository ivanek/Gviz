## This script will try to identify genome versions in all available
## ENSEBML archives as well as the current release, and match these
## versions back to UCSC genome identifiers.

library(biomaRt)

parseArchiveDate <- function(x){
    datePart <- sub("\\.archive.*", "", sub("^http://", "", x))
    return(as.Date(paste(datePart, "1", sep=""), format="%b%Y%d"))
}


getBmFromDate <- function(date){
    archive <- sprintf("%s%s.archive.ensembl.org", format(date, "%b"), format(date, "%Y"))
    res <- try(capture.output(suppressMessages(bm <- useMart(host=archive,
                                                             biomart="ENSEMBL_MART_ENSEMBL"))),
              silent=TRUE)
    if(is(res, "try-error")){
        bm <- NA
    }
    return(bm)
}

year <- function(x){
    if(missing(x)){
        x <- Sys.Date()
    }
    return(as.integer(format(x, "%Y")))
}

month <- function(x){
    if(missing(x)){
        x <- Sys.Date()
    }
    return(as.integer(format(x, "%m")))
}

day <- function(x){
    if(missing(x)){
        x <- Sys.Date()
    }
    return(as.integer(format(x, "%d")))
}

knownEnsemblArchives <- c("http://Jul2016.archive.ensembl.org",
                          "http://Mar2016.archive.ensembl.org",
                          "http://Dec2015.archive.ensembl.org",
                          "http://Sep2015.archive.ensembl.org",
                          "http://Jul2015.archive.ensembl.org",
                          "http://May2015.archive.ensembl.org",
                          "http://Mar2015.archive.ensembl.org",
                          "http://Dec2014.archive.ensembl.org",
                          "http://Oct2014.archive.ensembl.org",
                          "http://Aug2014.archive.ensembl.org",
                          "http://Feb2014.archive.ensembl.org",
                          "http://Dec2013.archive.ensembl.org",
                          "http://Sep2013.archive.ensembl.org",
                          "http://May2012.archive.ensembl.org",
                          "http://May2009.archive.ensembl.org")


parseOneArchive <- function(date){
    bm <- getBmFromDate(date)
    datasets <- listDatasets(bm)
}

parseEnsemblArchives <- function(known){
    archiveDates <- sort(parseArchiveDate(known))
    ## We may have to extend our search date range in case that the known archives vector is
    ## out of date and later archives exist
    latest <- max(archiveDates)
    daysBack <- as.integer(Sys.Date() - latest)
    if(daysBack > day()){
        daysSince <- Sys.Date() - (1:daysBack)
        archiveDates <- sort(unique(c(archiveDates,
                                      as.Date(paste(unique(format(daysSince, "%b%Y")), "1", sep=""),
                                              format="%b%Y%d"))))
    }



}








