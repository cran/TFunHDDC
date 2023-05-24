\name{genModelFD}
\alias{genModelFD}
\title{
genModelFD
}
\description{
Generate functional data with coefficients distributed according to a finite mixture of contamined normal distributions such that for the \eqn{\textit{k}}th cluster we have the
multivariate contaminated normal distribution with density 
\deqn{
f(\gamma_i;\theta_k)=\alpha_k\phi(\gamma_i;\mu_k,\Sigma_k)+(1-\alpha_k)\phi(\gamma_i;\mu_k,\eta_k\Sigma_k)}
where \eqn{\alpha_k\in (0.5,1)} represents the proportion of uncontaminated data, \eqn{\eta_k>1} is the inflation coefficient due to outliers, and \eqn{\phi(\gamma_i;\mu_k,\Sigma_k)} is the density for the multivariate normal distribution \eqn{N(\mu_k,\Sigma_k)}.
}
\arguments{
  \item{ncurves}{
The number of curves total for the simulation.
}
  \item{nsplines}{
The number of splines to fit to the simulated data.
}
  \item{alpha}{
The proportion of uncontaminated data in each group.
}
  \item{eta}{
The inflation coefficient that measures the increase in variability due to the outliers. 
}
}

\usage{
genModelFD(ncurves=1000, nsplines=35, alpha=c(0.9,0.9,0.9),
           eta=c(10, 5, 15))
}
\details{
The data are generate from the model \eqn{FCLM[a_k, b_k,{\bf{Q}}_k,d_k,\alpha_k,\eta_k]}. 
The number of clusters is fixed to \eqn{K=3} and the mixing proportions are equal \eqn{\pi_1=\pi_2=\pi_3=1/3}. We consider the following values of the parameters

Group 1:\eqn{d=5}, \eqn{a=150}, \eqn{b=5}, \eqn{\mu=(1,0,50,100,0,\ldots,0)}

Group 2: \eqn{d=20}, \eqn{a=15}, \eqn{b=8}, \eqn{\mu=(0,0,80,0,40,2,0,\ldots,0)}

Group 3: \eqn{d=10}, \eqn{a=30}, \eqn{b=10}, \eqn{\mu=(0,\ldots,0,20,0,80,0,0,100)},

where \eqn{d} is the intrinsic dimension of the subgroups, \eqn{\mu} is the mean vector of size 70, \eqn{a} is the values of the \eqn{d}-first diagonal elements of \eqn{\mathbf{D}}, and \eqn{b} the value of the last \eqn{70-d}- elements. Curves as smoothed using 35 Fourier basis functions.
}
\value{
  \item{fd}{A functional data object representing the simulated data.}
  \item{groupd}{Group classifications for each curve.}
}
\references{
- Amovin-Assagba M, Gannaz I, Jacques J (2022) Outlier detection in multivariate
functional data through a contaminated mixture model. Comput Stat
Data Anal 174.
- Cristina Anton, Iain Smith Model-based clustering of functional data via mixtures of \eqn{t} distributions. Advances in Data Analysis and Classification (to appear).
}
\author{
Cristina Anton and Iain Smith
}

\examples{
# Univariate Contaminated Data
data <- genModelFD(ncurves=300, nsplines=35, alpha=c(0.9,0.9,0.9),
                  eta=c(10, 7, 17))
plot(data$fd, col = data$groupd)
clm <- data$groupd
}
\keyword{datagen}