library(Gviz)
library(parallel)
library(biomaRt)
library(Biobase)
library(rtracklayer)
library(RMySQL)

fetchAllGenomesFromUcsc <- function() {
    mydb <- dbConnect(MySQL(),
        user = "genome", dbname = "hgcentral",
        host = "genome-mysql.cse.ucsc.edu"
    )
    rs <- dbSendQuery(mydb, "select * from dbDb")
    out <- fetch(rs, n = -1)
    out <- out[out$active != 0, ] # only active
    cols <- c(
        "name", "description", "organism", "scientificName",
        "sourceName", "taxId", "orderKey"
    ) # smallest orderKey == rescent
    out <- out[, cols]
    out <- out[order(out$organism, out$orderKey), ]
    return(out)
}


testBiomartVersion <- function() {
    allGenomes <- union(
        read.delim(system.file(package = "Gviz", "extdata/biomartVersionsLatest.txt"))$ucscId,
        read.delim(system.file(package = "Gviz", "extdata/biomartVersionsNow.txt"))$ucscId
    )
    res <- mclapply(allGenomes, function(genome) {
        return(try(
            {
                map <- Gviz:::.ucsc2Ensembl(genome)
                res <- list()
                if (map$date == "head") {
                    bm <- useMart("ensembl", dataset = map$dataset)
                    ds <- listDatasets(bm)
                    mt <- ds[match(map$dataset, ds$dataset), "version"]
                    if (is.na(mt)) {
                        res <- list(current = genome, set = map$dataset, type = "head", cause = "error")
                    }
                    if (mt != map$value) {
                        res <- list(current = genome, set = map$value, setCurrent = mt, type = "head", cause = "mismatch")
                    }
                } else {
                    bm <- useMart(
                        host = sprintf("https://%s.archive.ensembl.org", tolower(sub(".", "", map$date, fixed = TRUE))),
                        biomart = "ENSEMBL_MART_ENSEMBL", dataset = map$dataset
                    )
                    ds <- listDatasets(bm)
                    mt <- ds[match(map$dataset, ds$dataset), "version"]
                    if (is.na(mt)) {
                        res <- list(current = genome, set = map$dataset, type = "archive", cause = "error")
                    }
                    if (mt != map$value) {
                        res <- list(
                            current = genome, set = map$value, setArchive = sub(".", " ", map$date, fixed = TRUE),
                            setVersion = map$version, setCurrent = mt, type = "archive", cause = "mismatch"
                        )
                    }
                }
                res
            },
            silent = TRUE
        ))
    })
    res <- res[listLen(res) > 0]
    res[sapply(res, is, "try-error")] <- list(type = "error")
    allFields <- unique(unlist(lapply(res, names)))
    dt <- as.data.frame(do.call(rbind, lapply(res, function(x) t(as.data.frame(unlist(x)[allFields])))),
        stringsAsFactors = FALSE
    )
    rownames(dt) <- NULL
    gens <- ucscGenomes()
    gens$version <- as.numeric(gsub("[a-zA-Z]*", "", gens$db))
    gensS <- split(gens, gsub("[0-9]*$", "", gens$db))
    highestVersion <- as.data.frame(t(sapply(gensS, function(x) unlist(x[which.max(x$version), ]))),
        stringsAsFactors = FALSE
    )
    rownames(gens) <- gens$db
    mt <- match(gsub("[0-9]*$", "", dt$current), rownames(highestVersion))
    tmp <- highestVersion[mt, c("db", "date", "name", "version")]
    colnames(tmp) <- paste("UCSCLatest", colnames(tmp), sep = "_")
    dt2 <- cbind(dt, tmp)
    tmp2 <- gens[dt2$current, c("db", "date", "name", "version")]
    colnames(tmp2) <- paste("UCSC", colnames(tmp2), sep = "_")
    dt3 <- cbind(dt2, tmp2)
    dt3 <- dt3[rowSums(is.na(dt3)) < ncol(dt3), ]
    # if(require(writexl)) {
    #     write_xlsx(list("Sheet1"=dt3), path="ensMappings.xlsx")
    # }
    write.table(dt3, file = "ensMappings.txt", sep = "\t", quote = FALSE)
    return(dt3)
}

# ucsc <- fetchAllGenomesFromUcsc()

# library(tidyverse)
# library(rvest)
# url <- "http://Dec2021.archive.ensembl.org/info/website/archives/assembly.html"
# tab <- url %>%
#     read_html() %>%
#     html_table(na.strings = "")
# tab <- tab[[1]] %>%
#     left_join(tab[[2]], by="")
# colnames(tab)[1] <- "scientificName"
# tab %>%
#     pivot_longer(-scientificName, names_to = "dateversion") %>%
#     dplyr::filter(!is.na(value))
#
#
# url <- "http://dec2021.archive.ensembl.org/info/about/species.html"
# species <- url %>%
#          read_html() %>%
#          html_table(na.strings = c("", "-"))
