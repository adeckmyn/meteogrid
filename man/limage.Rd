\name{limage}
\alias{limage}
\title{Visualisation of gridded data.}
\description{This function merges the possibilities of image and filled.contour, and adds a flexible legend. }
\usage{
limage(x=1:dim(z)[1],y=1:dim(z)[2],z,smooth=FALSE,
           nlevels=15,levels=pretty(zlim,nlevels),
           color.palette=colorRampPalette(c("blue","white","red")),
           col=color.palette(length(levels)-1), na.col=par("bg"),
           legend=FALSE,legend.cex=.7,
           legend.width=1/12,legend.sep=c(1/4,1/2),
           legend.skip=1,legend.digits=5,
           plot.title,title.adjust=TRUE,
           legend.title,
           asp=1, useRaster=TRUE, ...) 
}
\arguments{
\item{x,y}{Vectors with regularly spaced co-ordinates.}
\item{z}{Matrix of data values to be plotted.}
\item{smooth}{If FALSE, the image consists of the actual grid boxes (cfr. image), if TRUE the values are smoothed (cfr. filled.contour). Note that in the first case, the domain size is extended by half a grid interval in all directions!}
\item{asp}{Aspect ratio of the axes. Default is 1. In the case of a
  Lat-Lon map, a warning is issued, because Lat and Lon degrees are not
  proportional in general.}
\item{nlevels, levels}{Number of levels (assumed to be evenly distributed) or explicit break points.}
\item{color.palette}{A color palette to use for the plot.}
\item{na.col}{Only used if \code{smooth=FALSE}. Sets an alternative colour for all NA values.}
\item{legend}{If TRUE, a legend is added to the right of the image. Note that in this case the title of the image will be off-centre if \code{title.adjust}=FALSE.} 
\item{legend.width}{The fraction of the plot region occupied by the xlegend. Must be smaller than 1.}
\item{legend.sep}{A vector with 2 elements that describes which fraction of the legend width is occupied by the color strip itself. Default is c(1/4,2/4).}
\item{legend.cex}{The scaling of the legend. Default is 0.6 .}
\item{legend.skip}{An integer >=1. In the legend, numerical values are given explicitely only with interval legend.skip. If there are many levels, you may want to set it to a larger value.}
\item{legend.title, plot.title}{These options define the main title and an optional title for the legend. They are not given simply as a character string, but in the form of a title command, e.g. \code{legend.title=title(main="legend")}. In this way, you can also have a subtitle etc. and chooose where to write the legend title.}
\item{useRaster}{A very good option for 'image'. The postscript output is much better. But it may be a problem with other graphical drivers.}
\item{col, legend.digits, title.adjust}{fine-tuning the plot.}
\item{...}{Other graphical parameters for e.g. image, ... }
}
\details{There are also shortcuts \code{iview, fcview} that go a little further by adding geographical aspects like map boundaries and have different colour defaults.
}
\seealso{
  \code{\link[graphics]{image}},\code{\link[graphics]{contour}},
  \code{\link[graphics]{filled.contour}},\code{\link{vecplot}}
  \code{\link{geofield}}
}
\keyword{file}
