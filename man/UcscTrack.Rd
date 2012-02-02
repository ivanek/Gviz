\name{UcscTrack}
\alias{UcscTrack}
\alias{clearSessionCache}

\title{
Meta-constructor for GenomeGraph tracks fetched directly from the
various UCSC data sources.
}

\description{
  
The UCSC data base provides a wealth of annotation information. This
function can be used to access UCSC, to retrieve the data available
there and to return it as an annotation track object ameanable to
plotting with \code{\link{plotTracks}}.

\code{clearSessionCache} is can be called to remove all cached items
from the session which are generated when connecting with the UCSC
data base.

}



\usage{
UcscTrack(track, table=NULL, trackType=c("AnnotationTrack",
"GeneRegionTrack", "DataTrack", "GenomeAxisTrack"), genome, chromosome,
name=NULL, from, to, ...)

clearSessionCache()
}

\arguments{
  
  \item{track}{Character, the name of the track to fetch from UCSC. To
    find out about available tracks please consult the online table
    browser at
    \url{http://genome.ucsc.edu/cgi-bin/hgTables?command=start}.}
  
  \item{table}{Character, the name of the table to fetch from UCSC, or
    \code{NULL}, in which case the default selection of tables is
    used. To find out about available tables for a given track please
    consult the online table browser at
    \url{http://genome.ucsc.edu/cgi-bin/hgTables?command=start}.}
  
  \item{trackType}{Character, one in \code{c("AnnotationTrack",
      "GeneRegionTrack", "DataTrack", "GenomeAxisTrack")}. The function
    will try to coerce the downloaded data in an object of this
    class. See below for details.}
  
  \item{genome}{Character, a valid USCS genome identifier for which to
    fetch the data.}
  
  \item{chromosome}{Character, a valid USCS character identifier for which to
    fetch the data.}
  
  \item{name}{Character, the name to use for the resulting track
    object.}

  \item{from, to}{A range of genomic locations for which to fetch data.}

  \item{\dots}{All additional named arguments are expected to be either
  display parameters for the resulting objects, or character scalars of
  column names in the downloaded UCSC data tables that are matched by
  name to available arguments in the respective constructor functions as
  defined by the \code{trackType} argument. See Details section for more
  information.}
}

\details{
  
  The data stored at the UCSC data bases can be of different formats:
  gene or transcript model data, simple annotation features like CpG
  Island locations or SNPs, or numeric data like conservation or
  mapability. This function presents a unified API to download all kinds
  of data and to map them back to one of the annotation track objects
  defined in this package. The type of object to hold the data has to be
  given in the \code{trackType} argument, and subsequently the function
  passes all data on to the respective object constructor. All
  additional named arguments are considered to be relevant for the
  constructor of choice, and single character scalars are replaced by
  the respective data columns in the dowloaded UCSC tables if
  available. For instance, assuming the table for track 'foo' contains
  the columns 'id', 'type', 'fromLoc' and 'toLoc', giving the featuer
  identifier, type, start end end location. In order to create an
  \code{\linkS4class{AnnotationTrack}} object from that data, we have to
  pass the additional named arguments \code{id="id"},
  \code{feature="type"}, \code{start="fromLoc"} and code{end="toLoc"} to the
  \code{UcscTrack} function. The complete function call could look like
  this:

  \code{UcscTrack(track="foo", genome="mm9", chromosome=3,
    from=1000, to=10000, trackType="AnnotationTrack", id="id",
    feature="type", start="from", end="to")}

  To reduce the bandwidth, some caching of the UCSC connection takes
  place. In order to remove these cached session items, call
  \code{clearSessionCache}.
  
}


\value{
  
  An annotation track object as determined by \code{trackType}.
  
}

\author{
Florian Hahne
}


\seealso{

  \code{\linkS4class{AnnotationTrack}}

  \code{\linkS4class{DataTrack}}

  \code{\linkS4class{GeneRegionTrack}}

  \code{\linkS4class{GenomeAxisTrack}}

  \code{\link{plotTracks}}

}






