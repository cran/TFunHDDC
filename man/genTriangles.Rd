\name{genTriangles}
\alias{genTriangles}
\title{
genTriangles
}
\description{
Generate contaminated triangle data. Groups 1, 2, 3, and 4 are separable over the two dimensions of functional data. Groups 5 and 6 contain the contaminated curves of groups 1 and 3 respectively.
}
\usage{
genTriangles()
}
\details{
Group 1:
	
	 \eqn{X_1(t) = U + (0.6 - U)H_1(t) + \epsilon_1(t)}
	 
	 \eqn{ X_2(t) = U + (0.5 - U)H_1(t) + \epsilon_1(t)}
	 
		Contaminated \eqn{X_1(t) = \sin(t) + (0.6 - U)H_1(t) + \epsilon_2(t)}
		
		Contaminated
		\eqn{  X_2(t) = \sin(t) + (0.5 - U)H_1(t) + \epsilon_2(t)}

Group 2:
	
		\eqn{X_1(t) = U + (0.6 - U)H_2(t) + \epsilon_1(t)}
		
		\eqn{X_2(t) = U + (0.5 - U)H_2(t) + \epsilon_1(t)}
	
Group 3:
	
		\eqn{X_1(t) = U + (0.5 - U)H_1(t) + \epsilon_1(t)}
		
		\eqn{ X_2(t) = U + (0.6 - U)H_2(t) + \epsilon_1(t)}
		
		 Contaminated \eqn{X_1(t) = \sin(t) + (0.5 - U)H_1(t) + \epsilon_3(t)}
		
		Contaminated \eqn{X_2(t) = \sin(t) + (0.6 - U)H_2(t) + \epsilon_3(t)}
	
Group 4:
	
		\eqn{X_1(t) = U + (0.5 - U)H_2(t) + \epsilon_1(t)}
		
		\eqn{X_2(t) = U + (0.6 - U)H_1(t) + \epsilon_1(t).}
Here \eqn{t\in [1,21]}, \eqn{H_1(t) = (6-\vert t-7\vert)_+}, and \eqn{H_2(t) = (6-\vert t-15\vert)_+},  with \eqn{(\cdot)_+} representing the positive part. \eqn{U \sim \mathcal{U}(0, 0.1)}, and \eqn{\epsilon_1(t)\sim N(0, 0.5)}, \eqn{\epsilon_2(t)\sim N(0, 2)}, \eqn{\epsilon_3(t) \sim Cauchy(0, 4)} are mutually independent white noises and independent of \eqn{U}. We simulate 100 curves for each group, groups 1 and 3 consisting of 80 ordinary curves and 20 contaminated curves. Curves are smoothed using a 25 cubic B-spline basis.
}
\value{
  \item{fd}{List of functional data objects representing the two dimensions of triangle data.}
  \item{groupd}{Group classification for each curve}
}
\references{
- C.Bouveyron and J.Jacques (2011), Model-based Clustering of Time Series in Group-specific Functional Subspaces, Advances in Data Analysis and Classification, vol. 5 (4), pp. 281-300,  <doi:10.1007/s11634-011-0095-6>

- Schmutz A, Jacques J, Bouveyron C, et al (2020) Clustering multivariate functional data in group-specific functional subspaces. Comput Stat
35:1101-1131

- Cristina Anton, Iain Smith Model-based clustering of functional data via mixtures of \eqn{t} distributions. Advances in Data Analysis and Classification (to appear).
}
\author{
Cristina Anton and Iain Smith
}
\seealso{
\code{\link{plotTriangles}}
}
\examples{
# Multivariate Contaminated Triangles
conTrig <- genTriangles()
cls = conTrig$groupd
plotTriangles(conTrig)
}
\keyword{datagen}
