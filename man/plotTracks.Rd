\name{plotTracks}
\alias{plotTracks}

\title{
The main plotting function for one or several GenomeGraph tracks.
}

\description{
  
  \code{plotTracks} is the main interface when plotting single track
  objects, or lists of tracks linked together across the same genomic
  coordinates. Essentially, the resulting plots are very similar to the
  graphical output of the UCSC Genome Browser, execpt for all of the
  interactivity.

}

\usage{

plotTracks(trackList, from=NULL, to=NULL, ..., sizes=NULL,
panel.only=FALSE, extend.right=0, extend.left=0, title.width=NULL,
add=FALSE, main, cex.main=2, fontface.main=2, col.main="black",
margin=6, chromosome=NULL)

}

\arguments{
  
  \item{trackList}{A list of GenomeGraph track objects, all inheriting from class
    \code{\linkS4class{GdObject}}. The tracks will all be drawn to the
    same genomic coordinates, either as defined by the \code{from} and
    \code{to} arguments if supplied, or by the maximum range across all
    individual items in the list.}

  \item{from, to}{Charactar scalar, giving the range of genomic
    coordinates to draw the tracks in. Note that \code{to} cannot be
    larger than \code{from}. If \code{NULL}, the plotting ranges are
    derived from the individual tracks. See \code{extend.left} and
    \code{extend.right} below for the definition of the final plotting
    ranges.}

  \item{\dots}{Additional arguments which are all interpreted as display
    parameters to tweak the appearance of the plot. These parameters are
    global, meaning that they will be used for all tracks in the list
    where they actually make sense, and they override the track-internal
    settings. See \code{\link{settings}} for details on display
    parameters.}

  \item{sizes}{A numeric vector of relative vertical sizes for the
    individual tracks of lenght
    equal to the number of tracks in \code{trackList}, or \code{NULL} to
    auto-detect the most appropriate vertical size proportions.}
  
  \item{panel.only}{Logical flag, causing the tracks to be plotted as
    lattice-like panel functions without resetting the plotting canvas
    and omitting the title pane. This allows to embed tracks into a
    trellis layout. Usually the function is called for a single track
    only when \code{panel.only==TRUE}.}

  \item{extend.right, extend.left}{Numeric scalar, extend the plotting
    range to the right or to the left by a fixed number of bases. The
    final plotting range is defined as \code{from-extend.left} to
    \code{to+extend.right}.}

  \item{title.width}{A expansion factor for the width of the title
    panels. This can be used to make more space, e.g. to accomodate for
    more detailed data axes. The default is to use as much space as
    needed to fit all the annotation text.}

  \item{add}{Logical flag, add the plot to an existing plotting canvas
    without re-initialising.}

  \item{main}{Character scalar, the plots main header.}
  
  \item{cex.main, fontface.main,col.main }{The fontface, color and
    expansion factor settings for the main header.}
  
  \item{margin}{The margin width to add to the plot in pixels.}
  
  \item{chromosome}{Set the chromosome for all the tracks in the track list.}
  
}

\details{
  
  GenomeGraph tracks are plotted in a vertically stacked layout. Each
  track panel is split up into a title section containing the track
  name, as well as an optional axis for tracks containing numeric data,
  and a data section showing the actual data along genomic
  coordinates. In that sense, the output is very similar to the UCSC
  Genome Browser.

  The layout of the individual tracks is highly customizable though so
  called "display parameters". See \code{\link{settings}} for details.

  While plotting a track, the software automatically computes HTML image
  map coordinates based on the current graphics device. These
  coordinates as well as the associated annotation information can later
  be used to embed images of the plots in semi-interactive HTML
  pages. See \code{\linkS4class{ImageMap}} for details. 

}


\value{

  A list of GenomeGraph tracks, each one augmented by the computed image
  map coordinates in the \code{imageMap} slot, along with the additional
  \code{ImageMap} object \code{titles} containing information about the
  title panels.  }

\author{
Florian Hahne
}


\seealso{

  \code{\linkS4class{GdObject}}

  \code{\linkS4class{ImageMap}}

  \code{\linkS4class{RangeTrack}}

  \code{\linkS4class{StackedTrack}}

  \code{\link{settings}}

}

\examples{

## Create some tracks to plot
st <- c(2000000, 2070000, 2100000, 2160000)
ed <- c(2050000, 2130000, 2150000, 2170000)
str <- c("-", "+", "-", "-")
gr <- c("Group1","Group2","Group1", "Group3")
annTrack <- AnnotationTrack(start=st, end=ed, strand=str, chromosome=7,
                            genome="hg19", feature="test", group=gr,
                            id=paste("annTrack item", 1:4),
                            name="annotation track foo",
                            stacking="squish")

ax <- GenomeAxisTrack()

dt <- DataTrack(start=seq(min(st), max(ed), len=10), width=18000,
data=matrix(runif(40), nrow=4), genome="hg19", chromosome=7,
type="histogram", name="data track bar")

\dontshow{
## For some annoying reason the postscript device does not know about
## the sans font
if(!interactive())
{
font <- ps.options()$family
displayPars(annTrack) <- list(fontfamily=font, fontfamily.title=font)
displayPars(ax) <- list(fontfamily=font, fontfamily.title=font)
displayPars(dt) <- list(fontfamily=font, fontfamily.title=font)
}
}


## Now plot the tracks
res <- plotTracks(list(ax, annTrack, dt))

## Plot only a subrange
res <- plotTracks(list(ax, annTrack, dt), from=2080000, to=2156000)

## Extend plotting ranges
res <- plotTracks(list(ax, annTrack, dt), extend.left=200000, extend.right=200000)

## Add a header
res <- plotTracks(list(ax, annTrack, dt), main="A GenomGraphs plot",
col.main="darkgray")

## Change vertical size and title width
res <- plotTracks(list(ax, annTrack, dt), sizes=c(1,1,5))

names(annTrack) <- "foo"
res <- plotTracks(list(ax, annTrack), title.width=0.6)

## Adding and lattice like plots
library(grid)
grid.newpage()
pushViewport(viewport(height=0.5, y=1, just="top"))
grid.rect()
plotTracks(annTrack, add=TRUE)
popViewport(1)
pushViewport(viewport(height=0.5, y=0, just="bottom"))
grid.rect()
plotTracks(dt, add=TRUE)
popViewport(1)

\dontrun{
library(lattice)
myPanel <- function(x, ...) plotTracks(annTrack, panel.only=TRUE,
from=min(x), to=max(x), shape="box")
a <- seq(1900000, 2250000, len=40)
xyplot(b~a|c, data.frame(a=a, b=1, c=cut(a, 4)), panel=myPanel,
scales=list(x="free"))
}
}





