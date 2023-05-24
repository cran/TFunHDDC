# h1 and h2 functions as defined in A. Schmutz et al. 2020.
.h1 <- function(t) ifelse(6-abs(t-7) > 0, 6-abs(t-7), 0)
.h2 <- function(t) ifelse(6-abs(t-15) > 0, 6-abs(t-15), 0)

#**************************** GENERATE SIMULATION *****************************#

genTriangles <- function(){
  #set.seed(5)
  t <- seq(from=0, to=21, length.out=101)
  h1v <- .h1(t)
  h2v <- .h2(t)

  curves <- matrix(ncol = 101, nrow = 400)
  curves2 <- matrix(ncol = 101, nrow = 400)
  group_1 <- 1:80
  group_2 <- 101:200
  group_3 <- 201:280
  group_4 <- 300:400

  for(i in group_1) {
    U <- runif(min=0, max=0.1, n=101)
    e <- rnorm(101, mean = 0, sd = sqrt(0.25))
    curves[i,] <- U + (0.6-U)*h1v + e

    U <- runif(min=0, max=0.1, n=101)
    e <- rnorm(101, mean = 0, sd = sqrt(0.25))
    curves2[i,] <- U + (0.5-U)*h1v + e
  }

  for(i in group_2) {
    U <- runif(min=0, max=0.1, n=101)
    e <- rnorm(101, mean = 0, sd = sqrt(0.25))
    curves[i,] <- U + (0.6-U)*h2v + e

    U <- runif(min=0, max=0.1, n=101)
    e <- rnorm(101, mean = 0, sd = sqrt(0.25))
    curves2[i,] <- U + (0.5-U)*h2v + e
  }

  for(i in group_3) {
    U <- runif(min=0, max=0.1, n=101)
    e <- rnorm(101, mean = 0, sd = sqrt(0.25))
    curves[i,] <- U + (0.5-U)*h1v + e

    U <- runif(min=0, max=0.1, n=101)
    e <- rnorm(101, mean = 0, sd = sqrt(0.25))
    curves2[i,] <- U + (0.6-U)*h2v + e
  }

  for(i in group_4) {
    U <- runif(min=0, max=0.1, n=101)
    e <- rnorm(101, mean = 0, sd = sqrt(0.25))
    curves[i,] <- U + (0.5-U)*h2v + e

    U <- runif(min=0, max=0.1, n=101)
    e <- rnorm(101, mean = 0, sd = sqrt(0.25))
    curves2[i,] <- U + (0.6-U)*h1v + e
  }

  contam_1 <- 81:100
  contam_3 <- 281:300

  sint <- sin((pi*t)/2)
  for(i in contam_1) {
    U <- runif(min=0, max=0.1, n=101)
    e1 <- rt(101, 4, 0)
    curves[i,] <- (0.6-U)*h1v + sint + e1

    U <- runif(min=0, max=0.1, n=101)
    e1 <- rt(101, 4, 0)
    curves2[i,] <- (0.5-U)*h1v + sint + e1
  }

  for(i in contam_3) {
    U <- runif(min=0, max=0.1, n=101)
    e2 <- rnorm(101, mean=0, sd=2)
    curves[i,] <- (0.5-U)*h1v + sint + e2 # function expects a t

    U <- runif(min=0, max=0.1, n=101)
    e2 <- rnorm(101, mean=0, sd=2)
    curves2[i,] <- (0.6-U)*h2v + sint + e2
  }



  title_name <- "T Distribution Contamination"
  # add more versions of outliers

  fdt <- list()
  basis <- create.bspline.basis(c(0, 21), nbasis=25)
  fdt[[1]] <- smooth.basis(argvals = seq( 0, 21,length.out = ncol(curves) ),
                           y = t(curves), fdParobj = basis)$fd
  fdt[[2]] <- smooth.basis(argvals = seq( 0, 21,length.out = ncol(curves2) ),
                           y = t(curves2), fdParobj = basis)$fd
  return( list( fd=fdt,
               groupd=c(rep(1, 80),
                        rep(5, 20),
                        rep(2, 100),
                        rep(3, 80),
                        rep(6, 20),
                        rep(4, 100)
              ) ) )
}


#******************************* PLOTTING DATA ********************************#

plotTriangles <- function(fdt) {
  # default values
  splits <- NA
  addcol <- TRUE
  clustd <- NA
  CL = c("black", "red", "green", "blue", "purple", "brown")
  
  if(is.list(fdt) && !is.null(fdt$fd)) {
    fd <- fdt$fd
  }else if(is.list(fdt) && is.null(fdt$fd)) {
    stop("plot_triangle_fd: missing $fd in list fdt.")
  }
  if(is.fd(fd)) {
    if(addcol) {
      if(is.na(clustd)) {
        clustd <- fdt$groupd
      }
      plot.fd(fd, type = 'l', col = CL[clustd])
    }else{
      plot.fd(fd, type = 'l', col = 'black')
    }
  }else if(is.list(fd)){
    # return user settings to original
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow = c(1,2))
    if(addcol) {
      if(is.na(clustd)) {
        clustd <- fdt$groupd
      }
      plot.fd(fd[[1]], col = CL[clustd], ylab = "X1(t)", xlab = "t",
           main = "a")
      plot.fd(fd[[2]], col = CL[clustd], ylab = "X2(t)", xlab = "t",
           main = "b")
    }else{
      plot.fd(fd[[1]], col = 'black', ylab = "X1(t)", xlab = "t")
      plot.fd(fd[[2]], col = 'black', ylab = "X2(t)", xlab = "t")
    }
    mtext(fdt$title_name, side = 3, line = - 2, outer = TRUE)
  } else {
    stop("plot_triangle_fd: fdt is not functional data or a list to extract functional data.")
  }
}


