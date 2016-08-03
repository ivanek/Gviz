library(Gviz)
library(biomaRt)
library(Biobase)
library(rtracklayer)


testBiomartVersion <- function(){
    allGenomes <- union(read.delim(system.file(package="Gviz", "extdata/biomartVersionsLatest.txt"))$ucscId,
                        read.delim(system.file(package="Gviz", "extdata/biomartVersionsNow.txt"))$ucscId)
    res <- mclapply(allGenomes, function(genome){
        return(try({
            map <- Gviz:::.ucsc2Ensembl(genome)
            res <- list()
            if(map$date == "head"){
                bm <- useMart("ensembl", dataset=map$dataset)
                ds <- listDatasets(bm)
                mt <- ds[match(map$dataset, ds$dataset), "version"]
                if(is.na(mt)){
                    res <- list(current=genome, set=map$dataset, type="head", cause="error")
                }
                if(mt != map$value){
                     res <- list(current=genome, set=map$value, setCurrent=mt, type="head", cause="mismatch")
                 }
            }else{
                bm <- useMart(host=sprintf("%s.archive.ensembl.org", tolower(sub(".", "", map$date, fixed=TRUE))),
                              biomart="ENSEMBL_MART_ENSEMBL", dataset=map$dataset)
                ds <- listDatasets(bm)
                mt <- ds[match(map$dataset, ds$dataset), "version"]
                if(is.na(mt)){
                    res <- list(current=genome, set=map$dataset, type="archive", cause="error")
                }
                if(mt != map$value){
                     res <- list(current=genome, set=map$value, setArchive=sub(".", " ", map$date, fixed=TRUE),
                                 setVersion=map$version, setCurrent=mt, type="archive", cause="mismatch")
                 }
            }
            res
        }, silent=TRUE))
    })
    res <- res[listLen(res)>0]
    res[sapply(res, is, "try-error")] <- list(type="error")
    allFields <- unique(unlist(lapply(res, names)))
    dt <- as.data.frame(do.call(rbind, lapply(res, function(x) t(as.data.frame(unlist(x)[allFields])))),
                        stringsAsFactors=FALSE)
    rownames(dt) <- NULL
    gens <- ucscGenomes()
    gens$version <- as.numeric(gsub("[a-zA-Z]*", "", gens$db))
    gensS <- split(gens, gsub("[0-9]*$", "", gens$db))
    highestVersion <- as.data.frame(t(sapply(gensS, function(x) unlist(x[which.max(x$version),]))),
                                    stringsAsFactors=FALSE)
    rownames(gens) <- gens$db
    mt <- match(gsub("[0-9]*$", "", dt$current), rownames(highestVersion))
    tmp <- highestVersion[mt, c("db", "date", "name", "version")]
    colnames(tmp) <- paste("UCSCLatest", colnames(tmp), sep="_")
    dt2 <- cbind(dt, tmp)
    tmp2 <- gens[dt2$current,c("db", "date", "name", "version")]
    colnames(tmp2) <- paste("UCSC", colnames(tmp2), sep="_")
    dt3 <- cbind(dt2, tmp2)
    if(require(xlsx)){
        write.xlsx2(dt3, file="ensMappings.xlsx", sheetName="Sheet1")
    }else{
        write.table(dt3, file="ensMappings.xls", sep="\t", quote=FALSE)
    }
    return(dt3)
}
