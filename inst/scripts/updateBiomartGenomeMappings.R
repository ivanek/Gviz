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

knownEnsemblArchives <- c(`Ensembl GRCh37`="http://grch37.ensembl.org",
                          `Ensembl 99: Jan 2020`="http://Jan2020.archive.ensembl.org",
                          `Ensembl 98: Sep 2019`="http://Sep2019.archive.ensembl.org",
                          `Ensembl 97: Jul 2019`="http://Jul2019.archive.ensembl.org",
                          `Ensembl 96: Apr 2019`="http://Apr2019.archive.ensembl.org",
                          `Ensembl 95: Jan 2019`="http://Jan2019.archive.ensembl.org",
                          `Ensembl 94: Oct 2018`="http://Oct2018.archive.ensembl.org",
                          `Ensembl 93: Jul 2018`="http://Jul2018.archive.ensembl.org",
                          `Ensembl 92: Apr 2018`="http://Apr2018.archive.ensembl.org",
                          `Ensembl 91: Dec 2017`="http://Dec2017.archive.ensembl.org",
                          `Ensembl 90: Aug 2017`="http://Aug2017.archive.ensembl.org",
                          `Ensembl 89: May 2017`="http://May2017.archive.ensembl.org",
                          `Ensembl 88: Mar 2017`="http://Mar2017.archive.ensembl.org",
                          `Ensembl 87: Dec 2016`="http://Dec2016.archive.ensembl.org",
                          `Ensembl 86: Oct 2016`="http://Oct2016.archive.ensembl.org",
                          `Ensembl 85: Jul 2016`="http://Jul2016.archive.ensembl.org",
                          `Ensembl 84: Mar 2016`="http://Mar2016.archive.ensembl.org",
                          `Ensembl 83: Dec 2015`="http://Dec2015.archive.ensembl.org",
                          `Ensembl 82: Sep 2015`="http://Sep2015.archive.ensembl.org",
                          `Ensembl 81: Jul 2015`="http://Jul2015.archive.ensembl.org",
                          `Ensembl 80: May 2015`="http://May2015.archive.ensembl.org",
                          `Ensembl 79: Mar 2015`="http://Mar2015.archive.ensembl.org",
                          `Ensembl 78: Dec 2014`="http://Dec2014.archive.ensembl.org",
                          `Ensembl 77: Oct 2014`="http://Oct2014.archive.ensembl.org",
                          `Ensembl 75: Feb 2014`="http://Feb2014.archive.ensembl.org",
                          `Ensembl 67: May 2012`="http://May2012.archive.ensembl.org",
                          `Ensembl 54: May 2009`="http://May2009.archive.ensembl.org")


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


