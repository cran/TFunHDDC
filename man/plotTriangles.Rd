\name{plotTriangles}
\alias{plotTriangles}
\title{
plotTriangles
}
\description{
Plot data returned by \code{\link{genTriangles}} as  lines coloured according to the assigned clusters.
}
\usage{
plotTriangles(fdt)
}
\arguments{
  \item{fdt}{
List with an element \eqn{fd} given the functional data, and an element \eqn{groupd} given the classes, as returned from the function \code{\link{genTriangles}}.
}
}
\value{No return value, used for side effects.}
\author{
Cristina Anton and Iain Smith
}
\seealso{
\code{\link{genTriangles}}
}
\examples{
conTrig <- genTriangles()
plotTriangles(conTrig)
}
