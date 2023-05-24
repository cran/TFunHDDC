#********************************* SIMULATIONS *********************************#
.newContaminatedSimulation <- function(ncurves, p, mu, K=3, prop=NULL,
                                      d=NULL, a=NULL, b=NULL, alpha=NULL,
                                      eta=NULL){
  ####################################################################################
  # Input: ncurves(int) -> number of curves to simulate                              #
  #        p(int) -> number of points to simulate for each curve                     #
  #        mu(matrix) -> matrix that contains mu for each group over all spline coef #
  #        K(int) : default 3 -> number of groups to simulate                        #
  # Purpose: Produce coefficients for simulated data over p splines it produces the  #
  #  sigma matrix.                                                                   #
  # Output: (list) -> $X(matrix) -> the coefficients simulated                       #
  #                   $clx(vector) -> expected classification                        #
  #                   $prms(list) -> parameters used to produce coefficients         #
  ####################################################################################
  #
  if (length(prop)==0) prop<-rep(1/K, K)
  else if (length(prop)!=K) stop("Proportions don't fit with the number of classes.")
  else prop<-prop/sum(prop)

  # Class sizes
  n<-floor(prop*ncurves)
  N<-sum(n)

  #MEANS
  #IS mu<-matrix(0, K, p) mu is now given
  j<-sample(p, K)
  
  # Intrinsic dimensions
  
  if ( length(d)==0 )	d<-sort(ceiling(runif(K, 0, 12*(p>20)+5*(p<=20 && p>=6)+(p<6)*(p-1))), decreasing=TRUE)
  else if ( length(d)!=K || !any(is.numeric(d)) ) stop("Wrong value of d.")

  # Orientation matrices
  Q<-vector(mode='list', length=K)
  for (i in 1:K) Q[[i]]<-qr.Q(qr(mvrnorm(p, mu=rep(0, p), Sigma=diag(1, p))))

  # Variance in the class-specific subspace
  if ( length(a)==0 ) a<-sort(ceiling(runif(K, 30, 350)))
  else if ( length(a)!=K || !any(is.numeric(a)) ) stop("Wrong value of a.")
  if ( length(b)==0 )b<-sort(ceiling(runif(K, 0, 25)))
  else if ( length(b)!=K || !any(is.numeric(b)) ) stop("Wrong value of b.")

  # Simulation
  S1<-vector(mode='list', length=K)
  S2<-vector(mode='list', length=K)
  
  #added eta[i] to define the contamination within groups
  for (i in 1:K)	{
    S1[[i]]<-crossprod(Q[[i]]%*%sqrt(diag(c(rep(a[i], d[i]), rep(b[i], p-d[i])))))
    S2[[i]]<-crossprod(Q[[i]]%*%sqrt(eta[i]*diag(c(rep(a[i], d[i]), rep(b[i], p-d[i])))))
  }


  cls<-X<-NULL
  clo <- rep(1, N)
  #contaminate sigma only on certain rows and all sigma for that
  for (i in 1:K) for(j in 1:n[i]){
    if(runif(0, 1, n=1) < alpha[i]) {
      X<-rbind(X, mvrnorm(1, mu=mu[i, ],Sigma=S1[[i]]))
    }else{
      clo[(i-1)*n[i]+j] <- 0
      X<-rbind(X, mvrnorm(1, mu=mu[i, ],Sigma=S2[[i]]))
    }
  }
  for (i in 1:K) cls<-c(cls, rep(i, n[i]))

  ind <- sample(1:N, N)
  prms <- list(a=a, b=b, prop=prop, d=d, mu=mu, ncurves=sum(n))
  return(list(X = X,
              clx = cls,
              clo = clo,
              prms = prms))

}
# create the mu matrix

genModelFD <- function(ncurves=1000, nsplines=35,
                       alpha=c(0.9,0.9,0.9), eta=c(10, 5, 15)) {
  ####################################################################################
  # Input: ncurves(int) -> number of curves to simulate                              #
  #        nsplines(int) -> number of splines in the produced functional data        #
  #        alpha(vector) -> amount of uncontaminated data in the curves              #
  #        eta(vector) : default 3 -> number of groups to simulate                   #
  # Purpose: Produce coefficients for simulated data over p splines it produces the  #
  #  sigma matrix.                                                                   #
  # Output: (list) -> $X(matrix) -> the coefficients simulated                       #
  #                   $clx(vector) -> expected classification                        #
  #                   $prms(list) -> parameters used to produce coefficients         #
  ####################################################################################
  mu <- rbind(c(1, 0, 50, 100, rep(0, nsplines-4)),
              c(0, 0, 80, 0, 40, 2, rep(0, nsplines-6)),
              c(rep(0, nsplines-6), 20, 0, 80, 0, 0, 100))
  a <- c(150, 15, 30)
  b <- c(5,8,10)
  d <- c(5,20,10)
  #set.seed(5)
  coef <- .newContaminatedSimulation(ncurves, nsplines, mu, d=d, a=a, b=b,
                                     eta = eta, alpha = alpha)
  simdata <- matrix(0, ncol=nsplines, nrow=coef$prms$ncurves)

  #basis <- create.bspline.basis(rangeval = c(0,100), nbasis = nsplines)
  basis <- create.fourier.basis(rangeval = c(0,100), nbasis = nsplines)
  evalbas <- eval.basis(seq(0, 100, length.out = 100), basis)
  finaldata <- coef$X %*% t(evalbas)
  finaldata <- cbind(finaldata, coef$clx)
#  bbasis <- create.bspline.basis(rangeval = c(0,100), nbasis = nsplines)
  bbasis <- create.fourier.basis(rangeval = c(0,100), nbasis = nsplines)
  fds <- smooth.basis(argvals = seq(0, 100, length.out = 100),
                      y = t(finaldata[,1:100]),
                      fdParobj = bbasis)$fd
  return(list(fd = fds,
              groupd = finaldata[, 101]))
}
