#***************************** REQUIRED LIBRARIES *****************************#
# library(fda.usc)
# library(fda)
if(getRversion() >= "2.15.1")  utils::globalVariables(c("poblenou"))

#********************************* FIT DATA ***********************************#
fitNOxBenchmark <- function(nbasis = 15) {

  data( "poblenou", package = "fda.usc", envir = environment() )
  nox_og <- poblenou$nox
  # smooth and fit data on functional basis
  basis <- create.bspline.basis(rangeval = seq(0, 23, length.out=13),
                                nbasis = nbasis)
  nox_fd <- smooth.basis(argvals = seq(0,23, length.out=24),
                         y = t( as.matrix(nox_og$data) ),
                         fdParobj = basis)$fd
  # modify classes for outlying festive days
  nox_class <- ifelse(poblenou$df$day.week %in% c(1,2,3,4,5) , 1, 2)
  nox_class <- ifelse(poblenou$df$day.festive %in% c(1), 2, nox_class)

  return( list(fd = nox_fd,
               groupd = nox_class) )
}

#******************************* PLOTTING DATA ********************************#
plotNOx <- function(fdn) {
  CL = c("1"="red", "2"="blue")
  fda::plot.fd(fdn$fd, col=CL[fdn$groupd])
}
