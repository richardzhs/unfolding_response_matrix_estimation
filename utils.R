## -------------------------------------------------------
## Utility functions
## This R file includes the utility functions for simulation, visualization,
## conditional density estimation, and response matrix estimation.
##

library(np)
library(ks)
library(locfit)
library(hrbrthemes)
library(ggplot2)
library(plotly)
library(pracma)
library(hdrcde)
library(parallel)
library(doParallel)
library(grid)
require(gridExtra)



################################################################################
# Functions for simulating data


#' generate simulation data for inclusive jet momentum with
#' the intensity function f0, and C1, C2, C3, delta, sd as the parameters
#' for the normal kernel function. It only works for marginal intensity function
#' that is monotonically decreasing in the range [lb, ub].
#' 
#' @param n number of data points
#' @param f0 density/intensity function on particle level
#' @param lb lower bound for the support of the pdf
#' @param ub upper bound for the support of the pdf
#' @param C1,C2,C3 parameters for the kernel function
#' @param delta the mean shifting parameter; the normal kernel will have mean
#' x + delta
#' @param sd standard deviation of the normal kernel; if it is specified, then
#' the normal kernel will use sd as the standard deviation and don't use C1, C2,
#' C3 for determining the standard deviation
#' 
#' @return a list of simulated x and y
#' 
generateData <- function(n, f0, lb, ub, C1=0, C2=0, C3=0, delta=0, sd=0) {
  # generate data on particle level (rejection sampling)
  M = f0(lb)
  x = rep(NA, n)
  i = 0
  while (i < n) {
    u = runif(n = 1, min = 0, max = 1)
    v = runif(n = 1, min = lb, max = ub)
    if (u <= f0(v) / M) {
      i = i + 1
      x[i] = v
    }
  }
  
  # propagate data into detector level
  y = rep(NA, n)
  for (i in 1:n) {
    if (sd != 0) {
      y[i] = rnorm(n = 1, mean = x[i] + delta, sd = sd)
    } else {
      y[i] = rnorm(n = 1, mean = x[i] + delta, sd = sqrt(x[i]*C1^2+C2^2+x[i]^2*C3^2))
    }
  }
  out = list()
  out[["x"]] = x
  out[["y"]] = y
  return(out)
}

generateData_general_k <- function(n, f0, lb, ub, k) {
  # generate data on particle level (rejection sampling)
  x <- rejection_sampling(n, f0, lb, ub)
  y = rep(NA, n)
  for (i in 1:n) {
    kx <- function(y) {
      k(x[i],y)
    }
    y[i] <- rejection_sampling(1, kx, lb-400, ub+400)
  }
  out = list()
  out[["x"]] = x
  out[["y"]] = y
  return(out)
}

#' Perform rejection sampling to sample from a given density with a bounded support.
#' 
#' @param n number of data points to sample
#' @param f0 density to sample from
#' @param lb,ub lower bound and upper bound of the support of the density
#' 
rejection_sampling <- function(n, f0, lb, ub) {
  M = optimize(f = f0, lower = lb, upper = ub, maximum = TRUE)$objective
  #M = max(f0(seq(lb,ub,by=1)))
  x = rep(NA, n)
  i = 0
  while (i < n) {
    u = runif(n = 1, min = 0, max = 1)
    v = runif(n = 1, min = lb, max = ub)
    if (u <= f0(v) / M) {
      i = i + 1
      x[i] = v
    }
  }
  return(x)
}

#' generate discretized simulation data for a given particle-level mean (lambda)
#' and response matrix (K)
#' 
#' @param lambda particle-level histogram mean (n x 1)
#' @param K response matrix (m x n)
#' 
#' @return a list of simulated x and y
#' 
generateDiscreteData <- function(lambda, K) {
  mu = as.vector(K %*% lambda)
  x = rpois(n = length(lambda), lambda = lambda)
  y = rpois(n = length(mu), lambda = mu)
  out = list()
  out[["x"]] = x
  out[["y"]] = y
  return(out)
}


#' Use Metropolis sampling to sample from a given pdf f
#' 
#' @param f pdf to sample from
#' @param n sample size
#' @param start initial point (default 0)
#' 
metrop_sampling <- function(f, n, start=0) {
  x = rep(NA, n)
  x[1] = start
  for (i in 2:n) {
    xstar = rnorm(1, mean=x[i-1], sd=1)
    alpha = f(xstar) / f(x[i-1])
    if (runif(1) < alpha) {
      x[i] = xstar
    } else {
      x[i] = x[i-1]
    }
  }
  return(x)
}


################################################################################
# Functions for ploting


#' Visualize response matrix K (rows represent true space and columns represent the 
#' smeared space) through heatmap.
#' 
#' @param K matrix
#' @param xbins binx on x
#' @param ybins bins on y
#' @param limits limits on the values from the elements of matrix (default NULL)
#' 
visualize <- function(K, xbins=NULL, ybins=NULL, limits = NULL, xlab="true pT", ylab="smeared pT",
                      legend_name="probability", title="", low_col = "#132B43", high_col = "#56B1F7") {
  if (is.null(xbins)) {
    xbins = 0:nrow(K)
  }
  if (is.null(ybins)) {
    ybins = 0:ncol(K)
  }
  dat = reshape2::melt(K)
  colnames(dat) <- c("Var1", "Var2", "probability")
  true = rep(xbins[-1], each = length(ybins)-1)
  smeared = rep(ybins[-1], times = length(xbins)-1)
  dat = cbind(dat, true, smeared)
  #grDevices::windowsFonts("Arial" = windowsFont("Arial"))
  p <- ggplot(dat, aes(x = true, y = smeared)) + 
    geom_tile(aes(fill = probability), size=1.2) +
    scale_fill_gradient(limits = limits, name = legend_name, low = low_col, high = high_col) +
    #scale_x_continuous(breaks = c(20, 100, 200, 300, 400, 500, 600), expand = 0) +
    #scale_y_continuous(breaks = c(20, 100, 200, 300, 400, 500, 600), expand = 0) +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(title)+
    coord_cartesian(xlim = c(xbins[1], xbins[length(xbins)]), 
                    ylim = c(ybins[1], ybins[length(ybins)]), expand = 0) +
    theme_ipsum(base_family = "Arial", axis_title_size = 12)
  
  ggplotly(p, tooltip="all")
  return(p)
}


#' Plot multiple 2d graphs in one plot
#' 
#' @param data a list of data to plot
#' @param xlim,ylim limits of the x-axis and yaxis
#' @param type type of plots to plot (default is points "p")
#' @param title title of the plot
#' @param legend title of the legend
#' @param xlab,ylab texts on the x-axis and y-axis
#' 
#' 
plot_multiple <- function(data, xlim=NULL, ylim=NULL, type="p", title="",
                          legend="", xlab="", ylab="") {
  for (i in 1:length(data)) {
    plot(data[[i]][[1]], data[[i]][[2]], xlim = xlim, ylim = ylim, 
         xlab = "", ylab = "", type = type, col = i)
    par(new=TRUE)
  }
  title(main = title, 
        xlab = xlab, ylab = ylab)
  legend(x = "topleft",
         legend = legend,
         col = 1:length(data), 
         lwd = 2)
}



################################################################################
## Functions for conditional density estimation

#' Compute the true conditional density evaluated at given (x,y).
#' 
#' @param p conditional denisty function
#' @param x,y points to evaluate
#' @return length(x) x length(y) matrix of evaluated density
compute_cd <- function(p, x, y, rescale=TRUE) {
  cd <- matrix(nrow = length(x), ncol = length(y))
  for (i in 1:nrow(cd)) {
    for (j in 1:ncol(cd)) {
      cd[i,j] = p(x[i],y[j])
    }
  }
  if(rescale)
  {
    delta <- y[2]-y[1]
    for(i in 1:nrow(cd))
    {
      sumz <- sum(cd[i,],na.rm=TRUE)
      if(sumz>0)
        cd[i,] <- cd[i,] / sum(cd[i,],na.rm=TRUE)
    }
    cd <- cd/delta
  }
  return(structure(list(x=x,y=y,z=cd),class="cde"))
}


#' Conditional Density Estimation 
#' (Adapted from hdrcde package: https://github.com/robjhyndman/hdrcde)
#'
#' Calculates kernel conditional density estimate using local polynomial
#' estimation.
#'
#' If bandwidths are omitted, they are computed using normal reference rules
#' described in Bashtannyk and Hyndman (2001) and Hyndman and Yao (2002). Bias
#' adjustment uses the method described in Hyndman, Bashtannyk and Grunwald
#' (1996). If deg>1 then estimation is based on the local parametric estimator
#' of Hyndman and Yao (2002).
#'
#' @param x Numerical vector or matrix: the conditioning variable(s).
#' @param y Numerical vector: the response variable.
#' @param deg Degree of local polynomial used in estimation.
#' @param link Link function used in estimation. Default "identity". The other
#' possibility is "log" which is recommended if degree > 0.
#' @param a Optional bandwidth in x direction.
#' @param b Optional bandwidth in y direction.
#' @param mean Estimated mean of y|x. If present, it will adjust conditional
#' density to have this mean.
#' @param x.margin Values in x-space on which conditional density is
#' calculated. If not specified, an equi-spaced grid of \code{nxmargin} values
#' over the range of x is used.  If x is a matrix, x.margin should be a list of
#' two numerical vectors.
#' @param y.margin Values in y-space on which conditional density is
#' calculated. If not specified, an equi-spaced grid of \code{nymargin} values
#' over the range of y is used.
#' @param x.name Optional name of x variable used in plots.
#' @param y.name Optional name of y variable used in plots.
#' @param use.locfit If TRUE, will use \code{\link[locfit]{locfit}} for
#' estimation. Otherwise \code{\link[stats]{ksmooth}} is used.
#' \code{\link[locfit]{locfit}} is used if degree>0 or link not the identity or
#' the dimension of x is greater than 1 even if \code{use.locfit=FALSE}.
#' @param fw If TRUE (default), will use fixed window width estimation.
#' Otherwise nearest neighbourhood estimation is used. If the dimension of x is
#' greater than 1, nearest neighbourhood must be used.
#' @param rescale If TRUE (default), will rescale the conditional densities to
#' integrate to one.
#' @param nxmargin Number of values used in \code{x.margin} by default.
#' @param nymargin Number of values used in \code{y.margin} by default.
#' @param a.nndefault Default nearest neighbour bandwidth (used only if
#' \code{fw=FALSE} and \code{a} is missing.).
#' @param \dots Additional arguments are passed to locfit.
#' @return A list with the following components: \item{x}{grid in x direction
#' on which density evaluated. Equal to x.margin if specified.} \item{y}{grid
#' in y direction on which density is evaluated. Equal to y.margin if
#' specified. } \item{z}{value of conditional density estimate returned as a
#' matrix. } \item{a}{window width in x direction.} \item{b}{window width in y
#' direction.} \item{x.name}{Name of x variable to be used in plots.}
#' \item{y.name}{Name of y variable to be used in plots.}
#' @author Rob J Hyndman
#' @seealso \code{\link{cde.bandwidths}}
#' @references Hyndman, R.J., Bashtannyk, D.M. and Grunwald, G.K. (1996)
#' "Estimating and visualizing conditional densities". \emph{Journal of
#' Computational and Graphical Statistics}, \bold{5}, 315-336.
#'
#' Bashtannyk, D.M., and Hyndman, R.J. (2001) "Bandwidth selection for kernel
#' conditional density estimation". \emph{Computational statistics and data
#' analysis}, \bold{36}(3), 279-298.
#'
#' Hyndman, R.J. and Yao, Q. (2002) "Nonparametric estimation and symmetry
#' tests for conditional density functions". \emph{Journal of Nonparametric
#' Statistics}, \bold{14}(3), 259-278.
#' @keywords smooth distribution hplot
#' @examples
#' # Old faithful data
#' faithful.cde <- cde(faithful$waiting, faithful$eruptions,
#'   x.name="Waiting time", y.name="Duration time")
#' plot(faithful.cde)
#' plot(faithful.cde, plot.fn="hdr")
#'
#' # Melbourne maximum temperatures with bias adjustment
#' x <- maxtemp[1:3649]
#' y <- maxtemp[2:3650]
#' maxtemp.cde <- cde(x, y,
#'   x.name="Today's max temperature", y.name="Tomorrow's max temperature")
#' # Assume linear mean
#' fit <- lm(y~x)
#' fit.mean <- list(x=6:45,y=fit$coef[1]+fit$coef[2]*(6:45))
#' maxtemp.cde2 <- cde(x, y, mean=fit.mean,
#' 	 x.name="Today's max temperature", y.name="Tomorrow's max temperature")
#' plot(maxtemp.cde)
#' @export cde
cde <- function(x, y, deg=0, link="identity", a, b, mean=NULL,
                x.margin,y.margin, x.name, y.name, use.locfit=FALSE, fw=TRUE, rescale=TRUE,
                nxmargin=15, nymargin=100, a.nndefault=0.3, maxk=100, ...)
{
  xname = deparse(substitute(x))
  yname = deparse(substitute(y))
  miss.xmargin=missing(x.margin)
  miss.ymargin=missing(y.margin)
  miss.a <- missing(a)
  miss.b <- missing(b)
  miss.xname <- missing(x.name)
  miss.yname <- missing(y.name)
  bias.adjust <- !is.null(mean)
  
  x <- as.matrix(x)
  nx <- ncol(x)
  
  if(bias.adjust & nx>1)
    stop("Bias adjustment not implemented for multiple conditioning variables")
  
  use.locfit <- (link!="identity" | deg>0 | use.locfit | nx>1 | !fw)
  fw <- (fw & nx==1)
  rescale <- (rescale | use.locfit)
  
  ## Get x.name
  if(miss.xname)
    x.name <- dimnames(x)[[2]]
  if(is.null(x.name))
  {
    if(nx==1)
      x.name <- xname
    else
      x.name <- paste(xname,"[,",1:nx,"]")
  }
  else if(length(x.name) != nx)
    stop("x.name has wrong length")
  dimnames(x) <- list(NULL,x.name)
  x.df <- as.data.frame(x)
  
  ## Get y.name
  if(miss.yname)
    y.name <- names(y)
  if(is.null(y.name))
    y.name <- yname
  else if(length(y.name)!=1)
    stop("y.name has wrong length")
  
  ##### Choose bandwidths
  if(miss.a | miss.b)
  {
    # Use reference rules
    # Only chooses bandwidth for first column of x.
    bands <- cdeband.rules(x[,1],y,deg=deg,link=link)
    if(miss.b)
      b <- bands$b
    if(miss.a)
    {
      if(fw)
        a <- bands$a ## Fixed width window
      else
        a <- a.nndefault  ## Nearest neighbourhood span
    }
  }
  if(use.locfit)
  {
    if(fw)
      locfit.a <- c(0,2.5*a)  ## For fixed bandwidth of a in locfit
    else
      locfit.a <- a
  }
  
  ##### Find y margin
  if(miss.ymargin)
  {
    yrange <- range(y)
    y.margin <- seq(yrange[1]-4*b,yrange[2]+4*b,l=nymargin)
  }
  else
    y.margin <- sort(y.margin)
  
  #### Find x margin
  if(miss.xmargin)
    x.margin <- NULL
  else if(is.matrix(x.margin)) # turn it into a list
    x.margin <- split(c(x.margin),rep(1:nx,rep(nrow(x.margin),nx)))
  else if(!is.list(x.margin))  #so is a vector
    x.margin <- list(x.margin)
  for(i in 1:nx)
  {
    if(miss.xmargin)
    {
      xrange <- range(x[,i])
      x.margin <- c(x.margin,list(seq(xrange[1],xrange[2],l=nxmargin)))
    }
    else
      x.margin[[i]] <- sort(x.margin[[i]])
  }
  names(x.margin) <- names(x.df)
  x.margin.grid <- expand.grid(x.margin)
  
  ##### Set up
  dim.cde <- c(length(y.margin),unlist(lapply(x.margin,length)))
  cde <- NULL
  GCV <- AIC <- numeric(dim.cde[1])
  n <- length(x)
  
  
  # If bias adjustment
  if(bias.adjust)
  {
    ymean <- mean(y)
    oldwarn <- options(warn=-1)
    approx.mean <- approx(mean$x,mean$y,xout=x)$y
    options(warn=oldwarn$warn)
    if(sum(is.na(approx.mean))>0)
      stop("Missing values in estimated mean")
    y <- y - approx.mean
    y.margin <- y.margin - ymean
  }
  
  
  ##### Do the calculations
  oldwarn <- options(warn=-1)
  xrange <- range(x[,1]) # How to handle multiple x??
  for(i in 1:dim.cde[1])
  {
    newy <- Kernel(y,y.margin[i],b,type="normal")
    if(max(abs(newy)) < 1e-20)
    {
      cde <- c(cde,rep(0,length(x.margin[[1]])))
    }
    else if(!use.locfit)
    {
      junk <- ksmooth(x[,1],newy,bandwidth=2.697959*a,kernel="normal",x.points=x.margin[[1]])$y
      junk[is.na(junk)] <- 0 ## No data in these areas
      cde <- c(cde,list(junk))
    }
    else
    {
      yscale <- mean(newy)
      newy <- newy/yscale
      junk <- locfit(newy ~ lp(x, nn = locfit.a, deg = deg), maxk=maxk)
      #junk <- locfit::locfit.raw(x,newy, alpha=locfit.a, deg=deg,link=link,family="qgauss",
      #kern="gauss",maxit=400,...)
      sum.coef <- sum(abs(junk$eva$coef))
      fits <- try(predict(junk,newdata=as.matrix(x.margin.grid)),silent=TRUE)
      if(class(fits)!="try-error")
      {
        AIC[i] <- -2 * junk$dp["lk"] + 2 * junk$dp["df2"]
        GCV[i] <- (-2 * n * junk$dp["lk"])/(n - junk$dp["df2"])^2
      }
      else
      {
        fits <- rep(NA,nrow(x.margin.grid))
        AIC[i] <- Inf  # or something huge
        GCV[i] <- Inf
      }
      cde <- c(cde,list(array(fits*yscale,dim.cde[-1])))
    }
  }
  options(warn=oldwarn$warn)
  
  AIC[AIC==Inf] <- 1.5*max(AIC[AIC<Inf])
  GCV[GCV==Inf] <- 1.5*max(GCV[GCV<Inf])
  z <- array(unlist(cde),dim=c(dim.cde[-1],dim.cde[1]))
  
  if(rescale)
  {
    delta <- y.margin[2]-y.margin[1]
    if(nx==1)
    {
      for(i in 1:dim.cde[2])
      {
        sumz <- sum(z[i,],na.rm=TRUE)
        if(sumz>0)
          z[i,] <- z[i,] / sum(z[i,],na.rm=TRUE)
      }
    }
    else if(nx==2)
    {
      for(i in 1:dim.cde[2])
        for(j in 1:dim.cde[3])
          z[i,j,] <- z[i,j,] / sum(z[i,j,],na.rm=TRUE)
    }
    z <- z/delta
  }
  z[z<0] <- 0
  
  # Bias adjustment
  if(bias.adjust)
  {
    #   browser()
    oldwarn <- options(warn=-1)
    approx.mean <- approx(mean$x,mean$y,xout=x.margin[[1]],rule=2)$y
    options(warn=oldwarn$warn)
    for(i in 1:dim.cde[2])
    {
      amean <- approx.mean[i] - sum(z[i,]*y.margin)*(y.margin[2]-y.margin[1])
      z[i,] <- approx(y.margin+amean,z[i,],xout=y.margin+ymean)$y
    }
    z[is.na(z)] <- 0
    y.margin <- y.margin + ymean
  }
  
  #    browser()
  
  ## Return the result
  if(nx==1)
    x.margin <- x.margin[[1]]  ## No need to keep it as a list.
  return(structure(list(x=x.margin,y=y.margin,z=z,a=a,b=b,deg=deg,link=link,
                        fn=switch(use.locfit+1,"ksmooth","locfit"),x.name=x.name, y.name=y.name,
                        fixed.width=fw,AIC=mean(AIC),GCV=mean(GCV),call=match.call()),class="cde"))
}


Kernel <- function(y,y0,b,type="epanech")
{
  if(type=="epanech")
    K <- epanech
  else
    K <- dnorm
  t(K(sweep(matrix(y0, nrow=length(y0), ncol=length(y)), 2, y),0,b))
}

"epanech" <- function(x,a,h)
{
  xx <- (x-a)/h
  0.75*(1-xx^2) * as.numeric(abs(xx)<1)
}


#' @export
print.cde <- function(x, ...)
{
  cat("Conditional density estimate:\n")
  cat(paste(x$y.name,"|",x$x.name,"\n\nCall: "))
  print(x$call)
  cat("\n  a=",x$a,"  b=",x$b)
  cat("\n  Degree=",x$deg,"  Link=",x$link,"\n")
}



# Conditional density estimation assuming location scale model
cde_ls <- function(x, y, xgrid, ygrid, rescale=TRUE) {
  cde <- matrix(nrow = length(xgrid), ncol = length(ygrid))
  #mu <- npreg(txdat = x, tydat = y, exdat = xgrid, residuals = TRUE)
  m1 <- smooth.spline(x, y)
  mu <- predict(m1, x=x)$y
  mu.eval <- predict(m1, x=xgrid)$y
  #m1 <- loess(y~x, data = data.frame(x,y), deg=1)
  #mu <- m1$fitted
  #mu.eval <- predict(m1, newdata = data.frame(x=xgrid))
  eps <- y - mu
  #sigma2 <- npreg(txdat = x, tydat = eps^2, exdat = xgrid)
  m2 <- smooth.spline(x, eps^2)
  sigma2 <- predict(m2, x)$y
  sigma2.eval <- predict(m2, xgrid)$y
  #m2 <- loess(eps2~x, data = data.frame(x,eps2=eps^2), deg=1)
  #sigma2 <- m2$fitted
  #sigma2.eval <- predict(m2, newdata = data.frame(x=xgrid))
  # evaluation points for the density estimates of epsilon
  eeps <- numeric(length(xgrid)*length(ygrid))
  k = 1
  for (i in 1:length(xgrid)) {
    for (j in 1:length(ygrid)) {
      eeps[k] = (ygrid[j] - mu.eval[i]) / sqrt(sigma2.eval[i])
      k = k + 1
    }
  }
  
  # remove the points where the estimated sigma^2 is < 0
  eps = eps[which(sigma2>0)]
  sigma2 = sigma2[which(sigma2>0)]
  
  m <- ks::kde(eps / sqrt(sigma2), eval.points = eeps)
  k = 1
  for (i in 1:length(xgrid)) {
    for (j in 1:length(ygrid)) {
      cde[i,j] = m$estimate[k] / sqrt(sigma2.eval[i])
      k = k + 1
    }
  }
  if(rescale)
  {
    delta <- ygrid[2]-ygrid[1]
    for(i in 1:nrow(cde))
    {
      sumz <- sum(cde[i,],na.rm=TRUE)
      if(sumz>0)
        cde[i,] <- cde[i,] / sum(cde[i,],na.rm=TRUE)
    }
    cde <- cde/delta
  }
  return(structure(list(x=xgrid, y=ygrid, z=cde), class="cde"))
}

#' compute kernel conditional density estimate with varying bandwidths depending
#' on x
#' 
#' @param x,y data points (x and y should have same length)
#' @param xgrid,ygrid grid of points along x-axis and y-axis to compute the cde
#' @param xbw pairs of bandwidths along x (should be a matrix of size length(xgrid) x 2)
#' @param rescale whether to rescale the density estimate to integrate to 1
#' 
#' @return a list of xgrid, ygrid, and matrix of density estimates
#' 
cde_xbw <- function(x, y, xgrid, ygrid, xbw, rescale=TRUE) {
  z <- matrix(nrow = length(xgrid), ncol = length(ygrid))
  for (i in 1:nrow(z)) {
    z[i,] <- cde(x = x, y = y, x.margin = xgrid[i], y.margin = ygrid, 
                       a = xbw[i,1], b = xbw[i,2], rescale = rescale)$z
    cat(xgrid[i]," ")
  }
  return(structure(list(x=xgrid, y=ygrid, z=z), class="cde"))
}

#' compute kernel conditional density estimate where only local points (neighborhood
#' of x) are used for estimation
#' 
#' @param x,y data points (x and y should have same length)
#' @param xgrid,ygrid grid of points along x-axis and y-axis to compute the cde
cde_local_x <- function(x, y, xgrid, ygrid, neighdist=NULL, d1=NULL, d2=NULL, 
                        bw_method=2, rescale=TRUE, parallel=FALSE) {
  z <- matrix(nrow = length(xgrid), ncol = length(ygrid))
  my_parallel <- function(i) {
    if (is.null(neighdist)) {
      if (is.null(d1) | is.null(d2)) {
        stop("d1 and d2 must not be null if neighdist is null")
      }
      neighdist = d1*exp(d2*xgrid[i])
    }
    sub_idx = which((x > xgrid[i] - neighdist) & (x < xgrid[i] + neighdist))
    x_sub = x[sub_idx]
    y_sub = y[sub_idx]
    #bw = hdrcde::cde.bandwidths(x = x_sub, y = y_sub, method = bw_method)
    out <- hdrcde::cde(x = x_sub, y = y_sub, x.margin = xgrid[i], y.margin = ygrid, 
                       rescale = rescale)$z
    cat(xgrid[i]," ")
    return(out)
  }
  if (parallel) {
    numCores <- detectCores()
    registerDoParallel(numCores-2)
    result <- foreach (k=1:nrow(z), .export=c("x", "y", "xgrid", "ygrid", "neighdist",
                                              "d1", "d2"), .combine=c) %dopar% {
      my_parallel(k)
    }
    for (i in 1:nrow(z)) {
      z[i,] <- result[(ncol(z)*i-ncol(z)+1):(ncol(z)*i)]
    }
  } else {
    for (i in 1:nrow(z)) {
      z[i,] = my_parallel(i)
    }
  }
  return(structure(list(x=xgrid, y=ygrid, z=z), class="cde"))
}


# a modified version of cde_local_x where we assume we know the true conditional
# density and hence can compute the optimal global bandwidth
cde_local_x_optim_bw <- function(x, y, xgrid, ygrid, neighdist, cd, rescale=TRUE, parallel=FALSE) {
  z <- matrix(nrow = length(xgrid), ncol = length(ygrid))
  my_parallel <- function(i) {
    sub_idx = which((x > xgrid[i] - neighdist) & (x < xgrid[i] + neighdist))
    x_sub = x[sub_idx]
    y_sub = y[sub_idx]
    bw = hdrcde::cde.bandwidths(x = x_sub, y = y_sub, y.margin = ygrid, method = bw_method)
    out <- hdrcde::cde(x = x_sub, y = y_sub, x.margin = xgrid[i], y.margin = ygrid, 
                       a = bw$a, b = bw$b, rescale = rescale)$z
    cat(xgrid[i]," ")
    return(out)
  }
  if (parallel) {
    numCores <- detectCores()
    registerDoParallel(numCores-2)
    result <- foreach (k=1:nrow(z), .export=c("x", "y", "xgrid", "ygrid"), .combine=c) %dopar% {
      my_parallel(k)
    }
    for (i in 1:nrow(z)) {
      z[i,] <- result[(ncol(z)*i-ncol(z)+1):(ncol(z)*i)]
    }
  } else {
    for (i in 1:nrow(z)) {
      z[i,] = my_parallel(i)
    }
  }
  return(structure(list(x=xgrid, y=ygrid, z=z), class="cde"))
}
  




###############################################################################
# Functions for response matrix estimation

#' Estimate the response matrix by binnings given the data and the specified bins
#' 
#' @param x,y data points (x and y should have same length)
#' @param bins binnings (end-points) of both the true space and smeared space
estimate_K_naive <- function(x, y, xbins, ybins) {
  n = length(x)
  nxbins = length(xbins) - 1
  nybins = length(ybins) - 1
  Kest = matrix(data = 0, nrow = nybins, ncol = nxbins)
  
  for (k in 1:n) {
    i = which((y[k] > ybins) == FALSE)[1] - 1
    j = which((x[k] > xbins) == FALSE)[1] - 1
    Kest[i,j] = Kest[i,j] + 1
  }
  xhist = hist(x,breaks = xbins,plot = FALSE)$count
  Kest = sweep(Kest,2,xhist,`/`)
  Kest[is.nan(Kest)] = 0
  return(Kest)
}


#' Estimate the response matrix using the estimated kernel (conditional density), 
#' estimated true spectrum, and the specified binnings.
#' 
#' @param k kernel function (current implementation assumes it is a matrix which 
#' the entries are evaluations k(x,y) at grids of points x,y)
#' @param f particle-level (true) spectrum (current implementation assumes it is a
#' vector which the elements are f(x) evaluated at points x)
#' @param x,y grids of (equi-distance) points x,y to be evaluated 
#' @param xbins binnings (end-points) of the true space (preferably, the binnings 
#' are contained in the grid of points)
#' @param ybins binnings (end-points) of the smeared space
#' @return the response matrix K
estimate_K <- function(k, x, y, f, xbins, ybins) {
  nxbins = length(xbins) - 1
  nybins = length(ybins) - 1
  K = matrix(0, nrow = nybins, ncol = nxbins)
  dx = abs(x[2] - x[1])
  dy = abs(y[2] - y[1])
  for (i in 1:nybins) {
    lby = which((y >= ybins[i]) == TRUE)[1]
    uby = tail(which((y < ybins[i+1]) == TRUE), n=1)
    for (j in 1:nxbins) {
      lbx = which((x >= xbins[j]) == TRUE)[1]
      ubx = tail(which((x < xbins[j+1]) == TRUE), n=1)
      temp = rep(0, (uby-lby+1))
      for (m in lby:uby) {
        for (n in lbx:ubx) {
          temp[m-lby+1] = temp[m-lby+1] + k[m,n]*f[n]*dx
        }
      }
      numerator = sum(temp * dy)
      denominator = 0
      for (n in lbx:ubx) {
        denominator = denominator + f[n]*dx
      }
      if (denominator == 0) {
        K[i,j] = 0
      } else {
        K[i,j] = numerator / denominator
      }
    }
  }
  return(K)
}


#' Compute the response matrix using the exactly known kernel (conditional density),
#' known true spectrum and the specified binnings.
#' 
#' @param k kernel function (an R function)
#' @param f particle-level (true) spectrum (an R function)
#' @param xbins binnings (end-points) of the true space (preferably, the binnings 
#' are contained in the grid of points)
#' @param ybins binnings (end-points) of the smeared space
#' @return the response matrix K
compute_K <- function(k, f, xbins, ybins) {
  nxbins = length(xbins) - 1
  nybins = length(ybins) - 1
  K = matrix(0, nrow = nybins, ncol = nxbins)
  
  for (j in 1:nxbins) {
    for (i in 1:nybins) {
      integrand <- function(x,y) {
        return(k(y,x)*f(x))
      }
      numerator = dblquad(integrand, xa=xbins[j], xb=xbins[j+1],
                          ya=ybins[i], yb=ybins[i+1])
      denominator = quadgk(f, a=xbins[j], b=xbins[j+1])
      if (denominator == 0) {
        K[i,j] = 0
      } else {
        K[i,j] = numerator / denominator
      }
    }
  }
  return(K)
}

#' Compute the true response matrix given the normal kernel function. It uses the
#' histogram estimate but with very large sample size to approximate the truth
#' 
#' @param C1,C2,C3 parameter for the variance of the normal kernel
#' @param delta parameter for the mean of the normal kernel (default 0)
#' @param f the true particle-level distribution
#' @param lb,ub lower bound and upper bound of the true particle-level distribution
#' @param n sample size (default 100000; if greater than 100000, it will break into
#' sub-sample of each size 100000)
#' @param xbins,ybins bins on the x and y space
#' 
#' 
#' @return the true response matrix K
compute_true_K <- function(C1, C2, C3, delta=0, f, n=100000, xbins, ybins) {
  K <- matrix(0, nrow=length(ybins)-1, ncol=length(xbins)-1)
  m <- n %/% 100000
  cat("m=",m,"\n")
  if (m == 0) {
    dat <- generateData(n, f, lb, ub, C1, C2, C3, delta)
    K <- estimate_K_naive(dat$x, dat$y, xbins, ybins)
    return(K)
  } else{
    for (i in 1:(m+1)) {
      dat <- generateData(100000, f, lb, ub, C1, C2, C3, delta)
      K <- K + 1/(m+1)*estimate_K_naive(dat$x, dat$y, xbins, ybins)
      cat(i, " ")
    }
  }
  return(K)
}


#' Compute the response matrix assuming normal kernel (conditional density),
#' known true spectrum and the specified binnings.
#' The normal kernels assumes the following structure:
#‘ k(y,x) = normal(y-x | mu=x+delta, sigma^2(x))
#' where (sigma(x)/x)^2 = (C1/sqrt(x))^2 + (C2/x)^2 + C3^2
#' 
#' @param delta parameters for the mean of the normal kernel
#' @param C1,C2,C3 parameters for the variance of the normal kernel
#' @param f particle-level (true) spectrum (an R function)
#' @param xbins binnings (end-points) of the true space (preferably, the binnings 
#' are contained in the grid of points)
#' @param ybins binnings (end-points) of the smeared space
#' @return the response matrix K
compute_K_normal <- function(C1, C2, C3, delta, f, xbins, ybins) {
  nxbins = length(xbins) - 1
  nybins = length(ybins) - 1
  K = matrix(0, nrow = nybins, ncol = nxbins)
  for (j in 1:nxbins) {
    for (i in 1:nybins) {
      kx <- function(x) {
        pnorm(ybins[i+1]-x,mean = 0,sd = sqrt(x*C1^2+C2^2+x^2*C3^2)) -
          pnorm(ybins[i]-x,mean = 0,sd = sqrt(x*C1^2+C2^2+x^2*C3^2))
      }
      integrand <- function(x) {
        return(kx(x)*f(x))
      }
      numerator = quadgk(integrand, a=xbins[j], b=xbins[j+1])
      denominator = quadgk(f, a=xbins[j], b=xbins[j+1])
      if (denominator == 0) {
        K[i,j] = 0
      } else {
        K[i,j] = numerator / denominator
      }
    }
  }
  return(K)
}
