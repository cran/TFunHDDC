\name{plotNOx}
\alias{plotNOx}
\title{
plotNOx
}
\description{
Plot data returned by \code{\link{fitNOxBenchmark}} as  lines coloured according to the assigned clusters.
}

\arguments{
  \item{fdn}{
List with an element \eqn{fd} given the functional data, and an element \eqn{groupd} given the classes; usually returned from the function \link{fitNOxBenchmark}.
}
}

\usage{
plotNOx(fdn)
}

\value{No return value, used for side effects.}
\author{
Cristina Anton and Iain Smith
}
\seealso{
\code{\link{fitNOxBenchmark}}
}
\examples{
# Univariate Contaminated Data
data1<-fitNOxBenchmark(15)
plotNOx(data1)
}