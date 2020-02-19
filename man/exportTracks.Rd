\name{exportTracks}
\alias{exportTracks}

\title{
  Export GenomeGraph tracks to a annotation file representation.
}

\description{

  This function is still a bit experimental. So far only BED export is
  supported.
  
}

\usage{

exportTracks(tracks, range, chromosome, file)

}


\arguments{

  \item{tracks}{A list of annotation track objects to be exported into a
  single BED file.}

  \item{range}{A numeric vector or length 2. The genomic range to
    display when opening the file in a browser.}
  
  \item{chromosome}{The chromosome to display when opening
    the file in a browser.}

  \item{file}{Character, the path to the file to write into.}
  
}

\details{
FIXME: Need to support WIG exports as well...
}

\value{
  
  The function is called for its side effect of writing to a file.

}

\examples{

## export AnnotationTrack to BED file
at <- AnnotationTrack(start=seq(1,10), width=1, chromosome="chr1")
exportTracks(list(at), range=c(1, 10), chromosome="chr1", 
             file=paste0(tempfile(),".bed"))

}


\author{
  Florian Hahne  
}
