\name{predict.tfunHDDC}
\alias{predict.tfunHDDC}
\title{
predict.tfunHDDC: Predicting Function for t-funHDDC Objects
}
\description{
Provides the matrix of classification probabilities and the classification vector for inputted observations assuming the model provided by the \code{\link{tfunHDDC}} object.
}
\usage{
\method{predict}{tfunHDDC}(object, data=NULL, \dots)
}
\arguments{
  \item{object}{
An object of class \code{\link{tfunHDDC}}
}
\item{data}{
Data frame (univariate funtional data) or a list (multivariate functional data) of new observations on the same variables used in the fitting of the \code{\link{tfunHDDC}} object.  If \code{NULL}, then the observations used in the fitting of the \code{\link{tfunHDDC}} object are inputted.
}

\item{\dots}{
Arguments to be passed to other functions.
}
}
\value{
\item{t}{Matrix of classification probabilities}
\item{class}{Vector of maximum a posteriori classifications}
}
\references{
-Andrews JL, McNicholas PD (2012) Model-based clustering, classication,
and discriminant analysis via mixtures of multivariate t-distributions:
The teigen family. Stat Comput 22:1021{1029. <doi.org/10.1007/
s11222-011-9272-x>

-Andrews JL, Wickins JR, Boers NM, et al (2018) An R package for modelbased
clustering and classication via the multivariate t distribution. Journal
of Statistical Software 83(7):1-32

- Cristina Anton, Iain Smith Model-based clustering of functional data via mixtures of \eqn{t} distributions. Advances in Data Analysis and Classification (to appear).
}
}
\author{
Cristina Anton, Iain Smith
}
\seealso{
\code{\link{tfunHDDC}}
}
\examples{
\donttest{
set.seed(1027)
#simulataed univariate data

data = genModelFD(ncurves=300, nsplines=35, alpha=c(0.9,0.9,0.9),
                  eta=c(10, 7, 17))

plot(data$fd, col = data$groupd)

clm = data$groupd

model1=c("AkjBkQkDk", "AkjBQkDk", "AkBkQkDk", "ABkQkDk", "AkBQkDk", "ABQkDk")

####################classification example with predictions

training=c(1:50,101:150, 201:250)

test=c(51:100,151:200, 251:300)

known1=clm[training]

t4<-tfunHDDC(data$fd[training],K=3,threshold=0.2,init="kmeans",nb.rep=1,
             dfconstr="no", dfupdate="numeric", model=model1[1],known=known1, 
             itermax = 10)

if (!is.null(t4$class)) {
	table(clm[training], t4$class)

	p1<-predict.tfunHDDC(t4,data$fd[test] )

	if (!is.null(p1$class)) table(clm[test], p1$class)
}

###########################NOX data

data1=fitNOxBenchmark(15)

plotNOx(data1)

###example for prediction

training=c(1:50)

test=c(51:115)

known1=data1$groupd[training]

t1<-tfunHDDC(data1$fd[training],K=2,threshold=0.6,init="kmeans",nb.rep=10,
             dfconstr="no", model=c("AkjBkQkDk", "AkjBQkDk", "AkBkQkDk", 
             "ABkQkDk", "AkBQkDk", "ABQkDk"),known=known1) 

if (!is.null(t1$class)) {
	table(data1$groupd[training], t1$class)

	p1<-predict.tfunHDDC(t1,data1$fd[test] )

	if (!is.null(p1$class)) table(data1$groupd[test], p1$class)
}
}
}
\keyword{classif}
