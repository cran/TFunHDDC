\name{fitNOxBenchmark}
\alias{fitNOxBenchmark}
\title{
fitNOxBenchmark
}
\description{
Extract NOx data from fda.usc
}
\arguments{
  \item{nbasis}{
The number of basis functions to fit to the simulated data.
}
}

\usage{
fitNOxBenchmark(nbasis=15)
}
\details{
Open NOx data from the poblenou data set of fda.usc. Fit the data to a given number of basis functions and
adjust classes for festive days into just weekdays and weekends.
}
\value{
  \item{fd}{A functional data object representing the fitted NOx data.}
  \item{groupd}{Group classifications for each curve as a curve representing a weekday or weekend/festive day.}
}
\references{
- Febrero-Bande M, Galeano P, Gonz~alez-Manteiga W (2008) Outlier detection
in functional data by depth measures, with application to identify abnormal
nox levels. Environmetrics 19:331-345. <//doi.org/10.1002/env.878>.

- Cristina Anton, Iain Smith Model-based clustering of functional data via mixtures of \eqn{t} distributions. Advances in Data Analysis and Classification (to appear).
}
\author{
Cristina Anton and Iain Smith
}
\seealso{
\code{\link{plotNOx}}
}
\examples{
# Univariate Contaminated Data
data1<-fitNOxBenchmark(15)
plotNOx(data1)
}
\keyword{datagen}