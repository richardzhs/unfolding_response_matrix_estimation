source("cde/bandwidth/bandwidth.selection.R")
source("utils.R")

## Setup: inclusive jet transverse momentum spectrum ---------------------------

# parameters for the intensity function
L = 5.1
N0 = 1e17
alpha = 5
beta = 10
gamma = 10
s = (7000)^2

# parameters for the smearing kernel
C1 = 1
C2 = 1
C3 = 0.05

# bounds for true and smeared space
lb = 400
ub = 1000

# number of bins (same for true and smeared space)
nbins = 40
bins = seq(lb, ub, length.out = nbins + 1)
nxbins = 40
xbins = seq(lb, ub, length.out = nxbins + 1)
nybins = 40
ybins = seq(lb, ub, length.out = nybins + 1)

# evaluation points for cde
xmargin = seq(400, 1000, length.out = 1000)
short.idx = seq(1,1000,by=10)
xmargin.short = xmargin[short.idx]
ymargin = seq(200, 1200, length.out = 1000)


## Theoretical functions & quantities ------------------------------------------

# particle-level intensity function (pt in units of GeV)
f0 <- function(pt) {
  L*N0*(pt)^(-alpha)*(1-2/sqrt(s)*pt)^beta*exp(-gamma/pt)
}
plot(f0, xlim = c(lb, ub))

# mean vector lambda for particle-level histograms
lambda = rep(NA, nbins)
for (i in 1:nbins) {
  lambda[i] = integrate(f = f0, lower = bins[i], upper = bins[i+1])$value
}


# response kernel
k1 <- function(x,y) {
  dnorm(y, mean = x, sd = sqrt(x*C1^2+C2^2+x^2*C3^2))
}

# Crystal ball function
CB <- function(y, x, mu, sigma, alpha, gamma) {
  m = x-y
  if ((m-mu)/sigma>-alpha){
    z = exp(-(m-mu)^2/(2*sigma^2))
  } else {
    z = (gamma/alpha)^gamma*exp(-alpha^2/2)*(gamma/alpha-alpha-(m-mu)/sigma)^(-gamma)
  }
  return(z)
}
CB <- Vectorize(CB)
k2 <- function(x,y) {
  CB(y, x, mu=0, sigma=x/50, alpha=1, gamma=4)
}

# smearing matrix K
Kmat = matrix(nrow = nybins, ncol = nxbins)
for (i in 1:nybins) {
  fi <- function(t) {
    (pnorm(ybins[i+1], mean = t, sd = sqrt(t*C1^2+C2^2+t^2*C3^2)) - 
       pnorm(ybins[i], mean = t, sd = sqrt(t*C1^2+C2^2+t^2*C3^2)))*f0(t)
  }
  for (j in 1:nxbins) {
    Kmat[i,j] = integrate(f = fi, lower = xbins[j], upper = xbins[j+1])$value/lambda[j]
  }
}
K <- compute_K_normal(C1,C2,C3,0,f0,xbins,ybins)
dat <- generateData(100000000, f0, lb, ub, C1, C2, C3)
xbins = seq(lb, ub, length.out = 10 + 1)
Ktrue10 <- compute_true_K(C1,C2,C3,0,f0,100000000,xbins,ybins)
saveRDS(Ktrue10, "cde/data/dat3/Ktrue_40x10.rds")
xbins = seq(lb, ub, length.out = 20 + 1)
Ktrue20 <- compute_true_K(C1,C2,C3,0,f0,100000000,xbins,ybins)
saveRDS(Ktrue20, "cde/data/dat3/Ktrue_40x20.rds")
xbins = seq(lb, ub, length.out = 40 + 1)
Ktrue40 <- compute_true_K(C1,C2,C3,0,f0,100000000,xbins,ybins)
saveRDS(Ktrue40, "cde/data/dat3/Ktrue_40x40.rds")

# mean vector mu for detector-level histograms
mu = as.vector(Kmat %*% lambda)


# visualizing K
visualize(Kmat,xbins = xbins, ybins = ybins)


## Toy Simulation --------------------------------------------------------------

n = 100000
out <- generateData(n, f0, lb, ub, C1, C2, C3)
x <- out[[1]]
y <- out[[2]]
plot(x,y)

# estimated K
Kest = matrix(data = 0, nrow = nbins, ncol = nbins)

for (k in 1:n) {
  i = which((y[k] > bins) == FALSE)[1] - 1
  j = which((x[k] > bins) == FALSE)[1] - 1
  Kest[i,j] = Kest[i,j] + 1
}
# calculate the proportion from the counts
Kest = sweep(Kest,2,colSums(Kest),`/`)
colSums(Kest)
Kest = estimate_K_naive(x,y,bins,bins)
colSums(Kest)

# visualizing Kest
visualize(Kest, xbins = bins, ybins = bins, limits = c(0,0.5))



## marginal densities ----------------------------------------------------------

f1 <- function(pt) {
  L*N0*(pt)^(-alpha)*(1-2/sqrt(s)*pt)^beta*exp(-gamma/pt)
}

p1 <- function(pt) {
  f1(pt)/integrate(f = f1, lower = lb, upper = ub)$value
}

f2 <- function(pt) {
  L*5.5e19*(pt)^(-6)*(1-2/sqrt(s)*pt)^12*exp(-gamma/pt)
}

f3 <- function(pt) {
  dunif(pt, min = lb, max = ub)
}

plot(f1, xlim = c(lb, ub), ylim = c(0, 15000))
par(new = TRUE)
plot(f2, xlim = c(lb, ub), ylim = c(0, 15000), col = 2)
par(new = TRUE)
plot(f3, xlim = c(lb, ub), ylim = c(0, 15000), col = 3)



## generating data -------------------------------------------------------------

dat1 = generateData(50000, f1, lb, ub, C1, C2, C3)
dat2 = generateData(100000, f1, lb, ub, C1, C2, C3)
dat3 = generateData(1000000, f1, lb, ub, C1, C2, C3)
dat4 = generateData(10000, f1, lb, ub, C1, C2, C3)
saveRDS(dat1, "cde/data/dat3/dat1.rds")
saveRDS(dat2, "cde/data/dat3/dat2.rds")
saveRDS(dat3, "cde/data/dat3/dat3.rds")
saveRDS(dat4, "cde/data/dat3/dat4.rds")



## Different cde methods -------------------------------------------------------
dat = readRDS("cde/data/dat3/dat1.rds")
dat = generateData(100000, f1, lb, ub, C1, C2, C3)
cd <- compute_cd(k1, x = xmargin, y = ymargin, rescale = TRUE)
opt.bw <- find_oracle_global_bw_sample(cd, 5,6,5,6,dat$x,dat$y,numCores = 2)
cde.out1 = cde_ls(dat$x, dat$y, xmargin, ymargin, rescale = TRUE)
cde.out2 = hdrcde::cde(x = dat$x, y = dat$y, x.margin = xmargin, y.margin = ymargin, rescale = TRUE)
#bw = find_oracle_global_bw_sample(cd, amin = 0.5, amax = 30, bmin = 0.5, bmax = 30,
#                                  x = dat$x, y = dat$y)
#saveRDS(bw, "bw.rds")
#cde.out3 = cde(x = dat$x, y = dat$y, x.margin = xmargin.short, y.margin = ymargin, 
#               a = bw$kernel[1], b = bw$kernel[2], rescale = TRUE)

find_optimal_window_width <- function(cd,d1,d2,dat,rescale=TRUE,parallel=FALSE) {
  mse = matrix(nrow=length(d1),ncol=length(d2))
  xmargin = cd$x
  ymargin = cd$y
  for (i in 1:length(d1)) {
    for (j in 1:length(d2)) {
      cde.out <- cde_local_x(dat$x, dat$y, xmargin, ymargin, d1 = d1[i], d2 = d2[j],
                  rescale = rescale, parallel = parallel)
      mse[i,j] = mean((cd$z-cde.out$z)^2)
      cat("(", d1[i], ",", d2[j], ")  ")
    }
  }
  opt = which(mse==min(mse), arr.ind = TRUE)
  return(list("d1"=d1,"d2"=d2,"mse"=mse,"opt_d1"=d1[opt[1]],"opt_d2"=d2[opt[2]]))
}

find_suitable_window_widths <- function(d1,d2,dat,xmargin,ymargin,min_prop=0.05,max_prop=0.2,
                                        tol_fail=round(length(xmargin)*0.05)) {
  n = length(dat$x)
  out <- c()
  for (i in 1:length(d1)) {
    for (j in 1:length(d2)) {
      fail = 0
      for (k in 1:length(xmargin)) {
        neighdist = d1[i]*exp(d2[j]*xmargin[k])
        sub_idx = which((dat$x > xmargin[k] - neighdist) & (dat$x < xmargin[k] + neighdist))
        if (length(sub_idx) < min_prop*n | length(sub_idx) > max_prop*n) {
          fail = fail + 1
        }
      }
      if (fail < tol_fail) {
        out <- rbind(out, c(d1[i], d2[j]))
      }
    }
  }
  return(out)
}

d1.d2 = find_suitable_window_widths(seq(0.5,5,by=0.5),1/seq(50,500,by=50),dat,xmargin,ymargin,
                                    min_prop = 0.05, max_prop = 0.2, tol_fail = 10)
opt.d1.d2 = find_optimal_window_width(cd, seq(0.5,10,by=0.5),1/seq(50,700,by=50), dat, parallel = TRUE)
saveRDS(opt.d1.d2,"opt.d1.d2.rds")
opt.d1.0 = find_optimal_window_width(cd, c(5:100,seq(105,700,by=5)), 0, dat, parallel = TRUE)
saveRDS(opt.d1.0,"opt.d1.0.rds")

opt.d1.d2 = readRDS("opt.d1.d2.rds")
opt.d1.0 = readRDS("opt.d1.0.rds")
plot(opt.d1.0$d1, opt.d1.0$mse, type = "l", xlab = "d1", ylab = "mse",
     main = "MSE for constant window size for local cde")
library(scatterplot3d)
xd1 = rep(opt.d1.d2$d1, length(opt.d1.d2$d2))
yd2 = rep(opt.d1.d2$d2, each=length(opt.d1.d2$d1))
zmse = expand.grid(opt.d1.d2$mse)
scatterplot3d(xd1, yd2, zmse$Var1, color="blue", pch=20, angle=30, 
              xlab="d1", ylab="d2", zlab="mse")

cde.out3 = cde_local_x(dat$x, dat$y, xmargin, ymargin, d1 = 8, d2 = 0,
                       rescale = TRUE, parallel = FALSE)

#bw.x = find_oracle_bw_x_sample(cd = cd, amin = 0.5, amax = 30, bmin = 0.5, bmax = 30,
#                               x = dat$x, y = dat$y, parallel = TRUE)
#bw.x = readRDS("cde/data/dat1/bwx.rds")
#cde.out5 = cde.xbw(dat$x, dat$y, xmargin, ymargin, bw.x$kernel)
#bw.eebs.short = readRDS("bandwidth/data/dat1/bw.eebs.short.rds")
#cde.out6 = cde.xbw(dat$x, dat$y, xmargin.short, ymargin, bw.eebs.short$bw)

#bw.lpcde = lpbwcde(x_data = dat$x, y_data = dat$y, x=800, y_grid = c(600,800),bw_type = "mse-rot")
#cde.out7 = lpcde(x_data = dat$x, y_data = dat$y, y_grid = ymargin, x = xmargin.short[1], bw = 2)


# toy sample for lpcde (seems working)
set.seed(42);
x_data = rnorm(2000)
y_data = rnorm(2000, mean=x_data)
x = 0
# Construct bandwidth
bw1 <- lpbwcde(y_data = y_data, x_data = x_data, x=x, bw_type = "mse-rot")
summary(bw1)
toy.cde = lpcde(x_data = x_data, y_data = y_data, x=x, bw=bw1$BW[,"bw"])
ylim=c(0,0.4)
plot(toy.cde$Estimate[,1],toy.cde$Estimate[,3], ylim=ylim)
par(new=TRUE)
plot(toy.cde$Estimate[,1], dnorm(toy.cde$Estimate[,1]), col='red', ylim=ylim)
lines(toy.cde$Estimate[,1], dnorm(toy.cde$Estimate[,1]), col='red')

plot(cd)
plot(cde.out1)
plot(cde.out2)
plot(cde.out3)
plot(cde.out5)


conditional_plot <- function(cde.out, cd, idx) {
  xlim = c(200, 1200)
  ylim = c(0,0.013)
  ylab = "z"
  plot(ymargin, cde.out$z[idx,], xlim = xlim, ylim = ylim, ylab = "", col = 2)
  par(new = TRUE)
  plot(ymargin, cd$z[idx,], xlim = xlim, ylim = ylim, ylab = ylab)
  legend(x = "topright",
         legend = c("true conditional density","cde"),
         col = c(1,2),
         lwd = 2)
}


conditional_plot(cde.out4, cd, 92)


flash <- function(cde.out, cd) {
  n = length(cde.out$x)
  for (idx in 1:n) {
    conditional_plot(cde.out, cd, idx)
  }
}


mse1.x = rowSums((cd$z-cde.out1$z)^2)*(cd$y[2]-cd$y[1])
mse2.x = rowSums((cd$z-cde.out2$z)^2)*(cd$y[2]-cd$y[1])
mse3.x = rowSums((cd$z-cde.out3$z)^2)*(cd$y[2]-cd$y[1])
mse4.x = rowSums((cd$z-cde.out4$z)^2)*(cd$y[2]-cd$y[1])
mse5.x = rowSums((cd$z-cde.out5$z)^2)*(cd$y[2]-cd$y[1])
mse6.x = rowSums((cd$z-cde.out6$z)^2)*(cd$y[2]-cd$y[1])

plot.data = list(list(cd$x,mse1.x),
                 list(cd$x,mse2.x),
                 list(cd$x,mse3.x))
                 #list(cd$x,mse4.x))
                 #list(cd$x,mse5.x),
                 #list(cd$x,mse6.x))

plot_multiple(plot.data, ylim = c(0,0.002), type = "l",
              legend = c("location scale","global normal reference",
                         "local x"))



## Integrated Mean Squared Error for CDE----------------------------------------
cd <- compute_cd(k1, x = xmargin, y = ymargin, rescale = TRUE)
niter = 200
n = 100000

numCores <- detectCores()
registerDoParallel(50)

mse <- foreach (k=1:niter, .combine = '+') %dopar% {
  dat = generateData(n, f1, lb, ub, C1, C2, C3)
  cde.out1 = cde_ls(dat$x, dat$y, xmargin, ymargin, rescale = TRUE)
  cde.out2 = hdrcde::cde(x = dat$x, y = dat$y, x.margin = xmargin, y.margin = ymargin, rescale = TRUE)
  cde.out3 = cde_local_x(dat$x, dat$y, xmargin, ymargin, d1 = opt.d1.d2$opt_d1, 
                         d2 = opt.d1.d2$opt_d2, rescale = TRUE, parallel = FALSE)
  cde.out4 = cde_local_x(dat$x, dat$y, xmargin, ymargin, d1 = opt.d1.0$opt_d1, d2 = 0,
                         rescale = TRUE, parallel = FALSE)
  #bw = hdrcde::cde.bandwidths(dat$x, dat$y, method=1)
  #cde.out4 = hdrcde::cde(x = dat$x, y = dat$y, x.margin = xmargin.short, y.margin = ymargin, 
  #                       rescale = TRUE, a = bw$a, b = bw$b)
  #cde.out5 = cde_local_x(dat$x, dat$y, xmargin.short, ymargin, bw_method = 1, 
  #                       d1 = 0.5, d2 = 1/100, rescale = TRUE, parallel = FALSE)
  cat(k, " ")
  out = array(NA, c(nrow(cd$z), ncol(cd$z), 4))
  out[,,1] = 1/niter*(cd$z-cde.out1$z)^2*(cd$y[2]-cd$y[1])*(cd$x[2]-cd$x[1])
  out[,,2] = 1/niter*(cd$z-cde.out2$z)^2*(cd$y[2]-cd$y[1])*(cd$x[2]-cd$x[1])
  out[,,3] = 1/niter*(cd$z-cde.out3$z)^2*(cd$y[2]-cd$y[1])*(cd$x[2]-cd$x[1])
  out[,,4] = 1/niter*(cd$z-cde.out4$z)^2*(cd$y[2]-cd$y[1])*(cd$x[2]-cd$x[1])
  return(out)
}

cond.mse <- foreach (k=1:niter, .combine = '+') %dopar% {
  dat = generateData(n, f1, lb, ub, C1, C2, C3)
  cde.out1 = cde_ls(dat$x, dat$y, xmargin, ymargin, rescale = TRUE)
  cde.out2 = hdrcde::cde(x = dat$x, y = dat$y, x.margin = xmargin, y.margin = ymargin, rescale = TRUE)
  cde.out3 = cde_local_x(dat$x, dat$y, xmargin, ymargin, d1 = opt.d1.d2$opt_d1, 
                         d2 = opt.d1.d2$opt_d2, rescale = TRUE, parallel = FALSE)
  cde.out4 = cde_local_x(dat$x, dat$y, xmargin, ymargin, d1 = opt.d1.0$opt_d1, d2 = 0,
                 rescale = TRUE, parallel = FALSE)
  #bw = hdrcde::cde.bandwidths(dat$x, dat$y, method=1)
  #cde.out4 = hdrcde::cde(x = dat$x, y = dat$y, x.margin = xmargin.short, y.margin = ymargin, 
  #                       rescale = TRUE, a = bw$a, b = bw$b)
  #cde.out5 = cde_local_x(dat$x, dat$y, xmargin.short, ymargin, bw_method = 1, 
  #                       d1 = 0.5, d2 = 1/100, rescale = TRUE, parallel = FALSE)
  cat(k, " ")
  return(cbind(1/niter*rowSums((cd$z-cde.out1$z)^2)*(cd$y[2]-cd$y[1]),
    1/niter*rowSums((cd$z-cde.out2$z)^2)*(cd$y[2]-cd$y[1]),
    1/niter*rowSums((cd$z-cde.out3$z)^2)*(cd$y[2]-cd$y[1]),
    1/niter*rowSums((cd$z-cde.out4$z)^2)*(cd$y[2]-cd$y[1])))
    #1/niter*rowSums((cd$z-cde.out5$z)^2)*(cd$y[2]-cd$y[1])))
}
saveRDS(mse, "mse.rds")

mse = readRDS("mse.rds")
cond.mse = readRDS("cond_mse.rds")
cond.mse.list <- list(list(cd$x, cond.mse[,1]), list(cd$x, cond.mse[,2]), list(cd$x, cond.mse[,3]), list(cd$x, cond.mse[,4]))


plot_multiple(cond.mse.list, ylim = c(0,0.001),legend = c("location-scale","normal reference",
                                                      "local normal reference (d1=1, d2=1/400)",
                                                     "local normal reference (d1=8, d2=0)"),
              xlab="x", ylab="mse", title="Conditional MSE for different CDE methods")


visualize(mse[,,1])


## Propagate to response matrix ------------------------------------------------
# one-sample data
dat = generateData(100000, f1, lb, ub, C1, C2, C3)
dat2 = generateData_general_k(100000, f1, lb, ub, k2)

# true conditional density
cd1 <- compute_cd(k1, xmargin, ymargin, rescale = TRUE)
cd2 <- compute_cd(k2, xmargin, ymargin, rescale = TRUE)
# location scale
cde.out1 = cde_ls(dat$x, dat$y, xmargin, ymargin, rescale = TRUE)
# kernel conditional density
cde.out2 = hdrcde::cde(x = dat$x, y = dat$y, x.margin = xmargin.short, y.margin = ymargin, rescale = TRUE)
# local kernel conditional density
cde.out3 = cde_local_x(dat2$x, dat2$y, xmargin.short, ymargin, d1 = 1, d2 = 1/400,
            rescale = TRUE, parallel = FALSE)
#cde.out4 = cde_local_x(dat$x, dat$y, xmargin, ymargin, d1 = 8, d2 = 0,
#                       rescale = TRUE, parallel = FALSE)
# local linear conditional density
cde.out4 = hdrcde::cde(x = dat$x, y = dat$y, x.margin = xmargin.short, y.margin = ymargin, 
                       rescale = TRUE, deg = 1)


xbins = seq(lb, ub, length.out = 40 + 1)
ybins = seq(lb, ub, length.out = 40 + 1)

# true response matrix
K <- compute_K_normal(C1, C2, C3, 0, f1, xbins, ybins)
cd2 <- compute_cd(k2,xmargin,ymargin)
K <- estimate_K(t(cd2$z), x=cd2$x, y=cd2$y, f=f1(cd2$x), xbins, ybins)
K1 <- estimate_K(k = t(cde.out1$z), x = cde.out1$x, y = cde.out1$y, f = f1(cde.out1$x), 
                 xbins = xbins, ybins = ybins)
K2 <- estimate_K(k = t(cde.out2$z), x = cde.out2$x, y = cde.out2$y, f = f1(cde.out2$x), 
                 xbins = xbins, ybins = ybins)
K3 <- estimate_K(k = t(cde.out3$z), x = cde.out3$x, y = cde.out3$y, f =  f1(cde.out3$x), 
                 xbins = xbins, ybins = ybins)
K4 <- estimate_K(k = t(cde.out4$z), x = cde.out4$x, y = cde.out4$y, f =  f1(cde.out4$x), 
                 xbins = xbins, ybins = ybins)
Knaive <- estimate_K_naive(dat2$x, dat2$y, xbins, ybins)
write.csv(K, "unfolding_osb_po_uq-main/data/jet_pt/config80/response_matrices/Ktrue.csv",row.names = FALSE)
write.csv(K1, "unfolding_osb_po_uq-main/data/jet_pt/config75/response_matrices_K=40x40_n=50000/K1.csv", row.names = FALSE)
write.csv(K2, "unfolding_osb_po_uq-main/data/jet_pt/config75/response_matrices_K=40x40_n=50000/K2.csv", row.names = FALSE)
write.csv(K3, "unfolding_osb_po_uq-main/data/jet_pt/config75/response_matrices_K=40x40_n=50000/K3.csv", row.names = FALSE)
write.csv(K4, "unfolding_osb_po_uq-main/data/jet_pt/config75/response_matrices_K=40x40_n=50000/K4.csv", row.names = FALSE)
write.csv(Knaive, "unfolding_osb_po_uq-main/data/jet_pt/config75/response_matrices_K=40x40_n=50000/Knaive.csv", row.names = FALSE)


# visualize K
Ktrue <- read.csv("../unfolding-main/data/jet_pt/normal_kernel/response_matrices/Ktrue.csv")
K1 <- read.csv("../unfolding-main/data/jet_pt/normal_kernel/response_matrices/K1-1.csv")
K2 <- read.csv("../unfolding-main/data/jet_pt/normal_kernel/response_matrices/K2-1.csv")
K3 <- read.csv("../unfolding-main/data/jet_pt/normal_kernel/response_matrices/K3-1.csv")
K4 <- read.csv("../unfolding-main/data/jet_pt/normal_kernel/response_matrices/K4-1.csv")
Knaive <- read.csv("../unfolding-main/data/jet_pt/normal_kernel/response_matrices/Knaive-1.csv")
xbins = seq(lb, ub, length.out = 40 + 1)
limits=c(-0.05,0.25)
p1 <- visualize(as.matrix(K1),xbins,ybins,limits=limits,title="location-scale")
p2 <- visualize(as.matrix(K2),xbins,ybins,limits=limits,title="kernel cde")
p3 <- visualize(as.matrix(K3),xbins,ybins,limits=limits,title="local kernel cde")
p4 <- visualize(as.matrix(K4),xbins,ybins,limits=limits,title="local linear cde")
ptrue <- visualize(as.matrix(K),xbins,ybins,limits=limits,title="true response matrix")
pnaive <- visualize(as.matrix(Knaive),xbins,ybins,limits=limits,title="histogram estimate")
grid.arrange(pnaive, p2, p4, p3, p1, ptrue, ncol=3,
             top = textGrob("Response Matrix Estimation (K=40x40, n=100000)",
                            gp=gpar(fontsize=20,font=3)))


# generate response matrices using different set of monte carlo samples
registerDoParallel(50)
nsamples <- 50
start.time = Sys.time()
foreach (i=1:100) %dopar% {
  dat = generateData_general_k(100000, f1, lb, ub, k2)
  cde.out1 = cde_ls(dat$x, dat$y, xmargin.short, ymargin, rescale = TRUE)
  cde.out2 = hdrcde::cde(x = dat$x, y = dat$y, x.margin = xmargin.short, y.margin = ymargin, rescale = TRUE)
  cde.out3 = cde_local_x(dat$x, dat$y, xmargin.short, ymargin, d1 = 1, d2 = 1/400,
                         rescale = TRUE, parallel = FALSE)
  cde.out4 = hdrcde::cde(x = dat$x, y = dat$y, x.margin = xmargin.short, y.margin = ymargin, 
                         rescale = TRUE, deg = 1)
  
  
  K1 <- estimate_K(k = t(cde.out1$z), x = cde.out1$x, y = cde.out1$y, f = f1(cde.out1$x), 
                   xbins = xbins, ybins = ybins)
  K2 <- estimate_K(k = t(cde.out2$z), x = cde.out2$x, y = cde.out2$y, f = f1(cde.out2$x), 
                   xbins = xbins, ybins = ybins)
  K3 <- estimate_K(k = t(cde.out3$z), x = cde.out3$x, y = cde.out3$y, f =  f1(cde.out3$x), 
                   xbins = xbins, ybins = ybins)
  K4 <- estimate_K(k = t(cde.out4$z), x = cde.out4$x, y = cde.out4$y, f =  f1(cde.out4$x), 
                   xbins = xbins, ybins = ybins)
  Knaive <- estimate_K_naive(dat$x, dat$y, xbins, ybins)
  
  
  write.csv(K1, strcat(c("resp_mat/K1-", i, ".csv")), row.names = FALSE)
  write.csv(K2, strcat(c("resp_mat/K2-", i, ".csv")), row.names = FALSE)
  write.csv(K3, strcat(c("resp_mat/K3-", i, ".csv")), row.names = FALSE)
  write.csv(K4, strcat(c("resp_mat/K4-", i, ".csv")), row.names = FALSE)
  write.csv(Knaive, strcat(c("resp_mat/Knaive-", i, ".csv")), row.names = FALSE)
  cat(i, " ")
}
end.time = Sys.time()
cat("spent:", end.time-start.time)

# true lambda
xhist <- compute_lambda_intensity(p1, 100000000, xbins)
write.csv(xhist, file = "unfolding_osb_po_uq-main/data/jet_pt/config77/true_histograms/x.csv", row.names = FALSE)
# true mu
mu <- K %*% xhist
write.csv(mu, file = "unfolding_osb_po_uq-main/data/jet_pt/config77/smeared_means/y.csv", row.names = FALSE)
# generate detector-level histograms
nsim <- 1000
for (i in 1:nsim) {
  out <- generateDiscreteData(xhist, K)
  y = out[[2]]
  write.csv(y, file = strcat(c("unfolding_osb_po_uq-main/data/jet_pt/config77/smeared_histograms/y", i, ".csv")), row.names = FALSE)
}


## Mean squared/absolute error for response matrices ------------------------------------
# number of iterations
niter = 1000
# number of data points in each iteration
n = 50000
xbins = seq(lb, ub, length.out = 40 + 1)
ybins = seq(lb, ub, length.out = 40 + 1)
# True matrix
K <- compute_K_normal(C1,C2,C3,0,f0,xbins,ybins)
K <- read.csv("unfolding_osb_po_uq-main/data/jet_pt/config76/response_matrices/Ktrue.csv")

mat_mse <- array(0, c(nrow(K), ncol(K), 5))
mat_mae <- array(0, c(nrow(K), ncol(K), 5))
for (i in 1:niter) {
  dat = generateData(n, f1, lb, ub, C1, C2, C3)
  # location scale
  cde.out1 = cde_ls(dat$x, dat$y, xmargin, ymargin, rescale = TRUE)
  # kernel conditional density
  cde.out2 = hdrcde::cde(x = dat$x, y = dat$y, x.margin = xmargin, y.margin = ymargin, rescale = TRUE)
  # local kernel conditional density
  cde.out3 = cde_local_x(dat$x, dat$y, xmargin, ymargin, d1 = 1, d2 = 1/400,
                         rescale = TRUE, parallel = FALSE)
  # local linear conditional density
  cde.out4 = hdrcde::cde(x = dat$x, y = dat$y, x.margin = xmargin, y.margin = ymargin, 
                         rescale = TRUE, deg = 1)
  
  K1 <- estimate_K(k = t(cde.out1$z), x = cde.out1$x, y = cde.out1$y, f = f1(cde.out1$x), xbins = xbins, ybins = ybins)
  K2 <- estimate_K(k = t(cde.out2$z), x = cde.out2$x, y = cde.out2$y, f = f1(cde.out2$x), xbins = xbins, ybins = ybins)
  K3 <- estimate_K(k = t(cde.out3$z), x = cde.out3$x, y = cde.out3$y, f = f1(cde.out3$x), xbins = xbins, ybins = ybins)
  K4 <- estimate_K(k = t(cde.out4$z), x = cde.out4$x, y = cde.out4$y, f = f1(cde.out4$x), xbins = xbins, ybins = ybins)
  Knaive <- estimate_K_naive(dat$x,dat$y,xbins,ybins)
  mat_mse[,,1] <- mat_mse[,,1] + 1/niter * as.matrix((K1-K)^2)
  mat_mse[,,2] <- mat_mse[,,2] + 1/niter * as.matrix((K2-K)^2)
  mat_mse[,,3] <- mat_mse[,,3] + 1/niter * as.matrix((K3-K)^2)
  mat_mse[,,4] <- mat_mse[,,4] + 1/niter * as.matrix((K4-K)^2)
  mat_mse[,,5] <- mat_mse[,,5] + 1/niter * as.matrix((Knaive-K)^2)
  cat(i, " ")
}

# or use the pre-generated response matrices to compute the mse (which is faster)
for (i in 1:niter) {
  K1 <- read.csv(strcat(c("unfolding_osb_po_uq-main/data/jet_pt/config76/response_matrices/K1-", i, ".csv")))
  K2 <- read.csv(strcat(c("unfolding_osb_po_uq-main/data/jet_pt/config76/response_matrices/K2-", i, ".csv")))
  K3 <- read.csv(strcat(c("unfolding_osb_po_uq-main/data/jet_pt/config76/response_matrices/K3-", i, ".csv")))
  K4 <- read.csv(strcat(c("unfolding_osb_po_uq-main/data/jet_pt/config76/response_matrices/K4-", i, ".csv")))
  Knaive <- read.csv(strcat(c("unfolding_osb_po_uq-main/data/jet_pt/config76/response_matrices/Knaive-", i, ".csv")))
  mat_mse[,,1] <- mat_mse[,,1] + 1/niter * as.matrix((K1-K)^2)
  mat_mse[,,2] <- mat_mse[,,2] + 1/niter * as.matrix((K2-K)^2)
  mat_mse[,,3] <- mat_mse[,,3] + 1/niter * as.matrix((K3-K)^2)
  mat_mse[,,4] <- mat_mse[,,4] + 1/niter * as.matrix((K4-K)^2)
  mat_mse[,,5] <- mat_mse[,,5] + 1/niter * as.matrix((Knaive-K)^2)
  mat_mae[,,1] <- mat_mae[,,1] + 1/niter * as.matrix(abs(K1-K))
  mat_mae[,,2] <- mat_mae[,,2] + 1/niter * as.matrix(abs(K2-K))
  mat_mae[,,3] <- mat_mae[,,3] + 1/niter * as.matrix(abs(K3-K))
  mat_mae[,,4] <- mat_mae[,,4] + 1/niter * as.matrix(abs(K4-K))
  mat_mae[,,5] <- mat_mae[,,5] + 1/niter * as.matrix(abs(Knaive-K))
  cat(i, " ")
}


# use parallel computing
numCores <- detectCores()
registerDoParallel(50)
mat_mse <- foreach (k=1:niter, .combine = '+') %dopar% {
  dat = generateData(n, f1, lb, ub, C1, C2, C3)
  # location scale
  cde.out1 = cde_ls(dat$x, dat$y, xmargin, ymargin, rescale = TRUE)
  # kernel conditional density
  cde.out2 = hdrcde::cde(x = dat$x, y = dat$y, x.margin = xmargin, y.margin = ymargin, rescale = TRUE)
  # local kernel conditional density
  cde.out3 = cde_local_x(dat$x, dat$y, xmargin, ymargin, d1 = 1, d2 = 1/400,
                         rescale = TRUE, parallel = FALSE)
  # local linear conditional density
  cde.out4 = hdrcde::cde(x = dat$x, y = dat$y, x.margin = xmargin, y.margin = ymargin, 
                         rescale = TRUE, deg = 1)
  
  K1 <- estimate_K(k = t(cde.out1$z), x = cde.out1$x, y = cde.out1$y, f = f1(cde.out1$x), xbins = xbins, ybins = ybins)
  K2 <- estimate_K(k = t(cde.out2$z), x = cde.out2$x, y = cde.out2$y, f = f1(cde.out2$x), xbins = xbins, ybins = ybins)
  K3 <- estimate_K(k = t(cde.out3$z), x = cde.out3$x, y = cde.out3$y, f = f1(cde.out3$x), xbins = xbins, ybins = ybins)
  K4 <- estimate_K(k = t(cde.out4$z), x = cde.out4$x, y = cde.out4$y, f = f1(cde.out4$x), xbins = xbins, ybins = ybins)
  Knaive <- estimate_K_naive(dat$x,dat$y,xbins,ybins)
  
  out = array(NA, c(nrow(K), ncol(K), 5))
  out[,,1] <- 1/niter*(K1-K)^2
  out[,,2] <- 1/niter*(K2-K)^2
  out[,,3] <- 1/niter*(K3-K)^2
  out[,,4] <- 1/niter*(K4-K)^2
  out[,,5] <- 1/niter*(Knaive-K)^2
  cat(k, " ")
  return(out)
}
saveRDS(mat_mse, "mat_mse.rds")

xbins = seq(lb, ub, length.out = 40 + 1)
ybins = seq(lb, ub, length.out = 40 + 1)
limits = c(0, 0.05)
low = "white"
high = "red"
p1 <- visualize(mat_mae[,,1],xbins=xbins,ybins=ybins, limits=limits,
          title = "location-scale", legend_name = "mae", low_col = low, high_col = high)
p2 <- visualize(mat_mae[,,2],xbins=xbins,ybins=ybins, limits=limits, 
          title = "kernel cde", legend_name = "mae", low_col = low, high_col = high)
p3 <- visualize(mat_mae[,,3],xbins=xbins,ybins=ybins, limits=limits, 
          title = "local kernel cde", legend_name = "mae", low_col = low, high_col = high)
p4 <- visualize(mat_mae[,,4],xbins=xbins,ybins=ybins, limits=limits, 
          title = "local linear cde", legend_name = "mae", low_col = low, high_col = high)
pnaive <- visualize(mat_mae[,,5],xbins=xbins,ybins=ybins, limits=limits,
          title = "histogram estimate", legend_name = "mae", low_col = low, high_col = high)


grid.arrange(p1, p2, p3, p4, pnaive, ncol=3,
             top = textGrob("MAE for Response Matrix Estimation for Crystal Ball Function (K=40x40, n=100000)",
                            gp=gpar(fontsize=20,font=3)))


ylim=c(0,0.2)
plot(bins[-1], colSums(mse1), ylim=ylim)
par(new=TRUE)
plot(bins[-1], colSums(mse2), ylim=ylim, col=2)
par(new=TRUE)
plot(bins[-1], colSums(mse3), ylim=ylim, col=3)
par(new=TRUE)
plot(bins[-1], colSums(msenaive), ylim=ylim, col=4)
legend(x = "topleft",
       legend = c("location-scale", "kernel", "kernel (optimal bw)", "naive"),
       col = c(1,2,3,4),
       lwd = 2)


## Bias calculation under different marginal density ---------------------------

mean1l <- matrix(0, nrow = nbins + 1, ncol = 100)
mean1k <- matrix(0, nrow = nbins + 1, ncol = 100)
mean1e <- matrix(0, nrow = nbins + 1, ncol = 100)
mean2l <- matrix(0, nrow = nbins + 1, ncol = 100)
mean2k <- matrix(0, nrow = nbins + 1, ncol = 100)
mean2e <- matrix(0, nrow = nbins + 1, ncol = 100)
mean3l <- matrix(0, nrow = nbins + 1, ncol = 100)
mean3k <- matrix(0, nrow = nbins + 1, ncol = 100)
mean3e <- matrix(0, nrow = nbins + 1, ncol = 100)

cdelist1l <- list()
cdelist1k <- list()
cdelist1e <- list()
cdelist2l <- list()
cdelist2k <- list()
cdelist2e <- list()
cdelist3l <- list()
cdelist3k <- list()
cdelist3e <- list()

N <- 50
n <- 50000


for (i in 1:N) {
  dat1 = generateData(n, f1, lb, ub, C1, C2, C3)
  cde1l <- cde(dat1$x, dat1$y, deg = 1, x.margin = bins, y.margin = ymargin, use.locfit = TRUE)
  cde1k <- cde(dat1$x, dat1$y, x.margin = bins, y.margin = ymargin)
  cde1e <- cde.eps(dat1$x, dat1$y, xgrid = bins, ygrid = ymargin)
  dat2 = generateData(n, f2, lb, ub, C1, C2, C3)
  cde2l <- cde(dat2$x, dat2$y, deg = 1, x.margin = bins, y.margin = ymargin, use.locfit = TRUE)
  cde2k <- cde(dat2$x, dat2$y, x.margin = bins, y.margin = ymargin)
  cde2e <- cde.eps(dat2$x, dat2$y, xgrid = bins, ygrid = ymargin)
  dat3 = generateData(n, f3, lb, ub, C1, C2, C3)
  cde3l <- cde(dat3$x, dat3$y, deg = 1, x.margin = bins, y.margin = ymargin, use.locfit = TRUE)
  cde3k <- cde(dat3$x, dat3$y, x.margin = bins, y.margin = ymargin)
  cde3e <- cde.eps(dat3$x, dat3$y, xgrid = bins, ygrid = ymargin)
  cdelist1l[[i]] <- cde1l$z
  cdelist1k[[i]] <- cde1k$z
  cdelist1e[[i]] <- cde1e$z
  cdelist2l[[i]] <- cde2l$z
  cdelist2k[[i]] <- cde2k$z
  cdelist2e[[i]] <- cde2e$z
  cdelist3l[[i]] <- cde3l$z
  cdelist3k[[i]] <- cde3k$z
  cdelist3e[[i]] <- cde3e$z
  mean1l = mean1l + 1/N*cde1l$z
  mean1k = mean1k + 1/N*cde1k$z
  mean1e = mean1e + 1/N*cde1e$z
  mean2l = mean2l + 1/N*cde2l$z
  mean2k = mean2k + 1/N*cde2k$z
  mean2e = mean2e + 1/N*cde2e$z
  mean3l = mean3l + 1/N*cde3l$z
  mean3k = mean3k + 1/N*cde3k$z
  mean3e = mean3e + 1/N*cde3e$z
  cat(i, " ")
}

save(cdelist1l, cdelist1k, cdelist2l, cdelist2k, cdelist3l, cdelist3k, 
     mean1l, mean1k, mean2l, mean2k, mean3l, mean3k,
     file = "cde list + mean estimate (no bw specified) (3.12).RData")

save(mean1k, mean1l, mean2k, mean2l, mean3k, mean3l, file = "mean estimate (3.1).RData")
load("cde list + mean estimate (no bw specified) (3.12).RData")
load("mean estimate (2.28 no bw specified).RData")



# Compute integrated bias given x
bias1l <- apply(X = (mean1l - cd$z)^2, 1, mean)
bias1k <- apply(X = (mean1k - cd$z)^2, 1, mean)
bias1e <- apply(X = (mean1e - cd$z)^2, 1, mean)
bias2l <- apply(X = (mean2l - cd$z)^2, 1, mean)
bias2k <- apply(X = (mean2k - cd$z)^2, 1, mean)
bias2e <- apply(X = (mean2e - cd$z)^2, 1, mean)
bias3l <- apply(X = (mean3l - cd$z)^2, 1, mean)
bias3k <- apply(X = (mean3k - cd$z)^2, 1, mean)
bias3e <- apply(X = (mean3e - cd$z)^2, 1, mean)


# bias comparison
ylim = c(0, 1e-5)
type = "p"
plot(bins, bias1l, ylim = ylim, col = 1, xlab = "", ylab = "", type = type, pch = 0)
par(new = TRUE)
plot(bins, bias1k, ylim = ylim, col = 2, xlab = "", ylab = "", type = type, pch = 0)
par(new = TRUE)
plot(bins, bias1e, ylim = ylim, col = 3, xlab = "", ylab = "", type = type, pch = 0)
par(new = TRUE)
plot(bins, bias2l, ylim = ylim, col = 1, xlab = "", ylab = "", type = type, pch = 1)
par(new = TRUE)
plot(bins, bias2k, ylim = ylim, col = 2, xlab = "", ylab = "", type = type, pch = 1)
par(new = TRUE)
plot(bins, bias2e, ylim = ylim, col = 3, xlab = "", ylab = "", type = type, pch = 1)
par(new = TRUE)
plot(bins, bias3l, ylim = ylim, col = 1, xlab = "", ylab = "", type = type, pch = 2)
par(new = TRUE)
plot(bins, bias3k, ylim = ylim, col = 2, xlab = "", ylab = "", type = type, pch = 2)
par(new = TRUE)
plot(bins, bias3e, ylim = ylim, col = 3, xlab = "", ylab = "", type = type, pch = 2)
par(xpd=TRUE)
title(main = "Bias (f1,f2: jet pT spectrum, f3: Uniform)", 
      xlab = "x", ylab = "bias")
legend(x = "topright",
       legend = c("f1-local", "f2-local", "f3-local", "f1-kernel", "f2-kernel",
                  "f3-kernel", "f1-add_error", "f2-add_error", "f3-add_error"),
       col = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
       pch = c(0, 1, 2, 0, 1, 2, 0, 1, 2),
       lty = rep(NA, 9),
       lwd = 2)


# Compute integrated variance given x
var1l <- ise1l - bias1l
var1k <- ise1k - bias1k
var2l <- ise2l - bias2l
var2k <- ise2k - bias2k
var3l <- ise3l - bias3l
var3k <- ise3k - bias3k

# variance comparison
ylim = c(0, 5e-6)
type = "p"
plot(bins, var1l, ylim = ylim, col = 1, xlab = "", ylab = "", type = type, pch = 0)
par(new = TRUE)
plot(bins, var1k, ylim = ylim, col = 1, xlab = "", ylab = "", type = type, pch = 1)
par(new = TRUE)
plot(bins, var2l, ylim = ylim, col = 2, xlab = "", ylab = "", type = type, pch = 0)
par(new = TRUE)
plot(bins, var2k, ylim = ylim, col = 2, xlab = "", ylab = "", type = type, pch = 1)
par(new = TRUE)
plot(bins, var3l, ylim = ylim, col = 3, xlab = "", ylab = "", type = type, pch = 0)
par(new = TRUE)
plot(bins, var3k, ylim = ylim, col = 3, xlab = "", ylab = "", type = type, pch = 1)
title(main = "Variance (f1,f2: jet momentum spectrum, f3: Uniform)", 
      xlab = "x", ylab = "variance")
legend(x = "topleft",
       legend = c("f1-local", "f1-kernel", "f2-local", "f2-kernel",
                  "f3-local", "f3-kernel"),
       col = c(1, 1, 2, 2, 3, 3),
       pch = c(0, 1, 0, 1, 0, 1),
       lty = rep(NA, 6),
       lwd = 2)



# Compute integrated mean squared error given x
ise1l <- numeric(length(bins))
ise1k <- numeric(length(bins))
ise2l <- numeric(length(bins))
ise2k <- numeric(length(bins))
ise3l <- numeric(length(bins))
ise3k <- numeric(length(bins))

for (k in 1:N) {
  cde1l <- cdelist1l[[k]]
  cde1k <- cdelist1k[[k]]
  cde2l <- cdelist2l[[k]]
  cde2k <- cdelist2k[[k]]
  cde3l <- cdelist3l[[k]]
  cde3k <- cdelist3k[[k]]
  
  ise1l = ise1l + 1/N * rowMeans((cde1l-cd)^2)
  ise1k = ise1k + 1/N * rowMeans((cde1k-cd)^2)
  ise2l = ise2l + 1/N * rowMeans((cde2l-cd)^2)
  ise2k = ise2k + 1/N * rowMeans((cde2k-cd)^2)
  ise3l = ise3l + 1/N * rowMeans((cde3l-cd)^2)
  ise3k = ise3k + 1/N * rowMeans((cde3k-cd)^2)
}

# 
# ylim = c(0, 0.012)
# xlim = c(200, 1400)
# type = "l"
# plot(ymargin, mean3l[5,], xlim = xlim, ylim = ylim, type = type)
# par(new=TRUE)
# plot(ymargin, cd[5,], xlim = xlim, ylim = ylim, col = 2, type = type)




## Bandwidth selection for conditional kernel density estimate -----------------

a = matrix(rep(seq(2,20,length.out=2),2), nrow = 2, byrow = TRUE)
b = matrix(rep(seq(2,20,length.out=2),2), nrow = 2, byrow = TRUE)

a = seq(2,20,length.out=2)
b = seq(2,20,length.out=2)
cde.collection <- cde_different_bw(dat$x, dat$y, a = a, b = b, xmargin=c(500,501), ymargin=c(600,601,602))

cd <- compute_cd(k1, xmargin.short, ymargin, rescale = TRUE)
# This is one sample of data
dat = readRDS("bandwidth/data/dat1/dat.rds")


# global optimal bandwidth for the specific sample
bw = readRDS("bandwidth/data/dat1/bw.rds")
# local optimal (conditioned on x) bandwidth for the specific sample 
bw.x = readRDS("bandwidth/data/dat1/bwx.rds")
#bw.x.temp = readRDS("bandwidth/data/dat1/bwx_temp.rds")

start = Sys.time()
bw1 = cde.bandwidths(dat$x, dat$y)
print(Sys.time()-start)
start = Sys.time()
bw2 = cde.bandwidths(dat$x, dat$y, method = 2)
print(Sys.time()-start)
start = Sys.time()
bw3 = cde.bandwidths(dat$x, dat$y, method = 3)
print(Sys.time()-start)
start = Sys.time()
bw4 = cde.bandwidths(dat$x, dat$y, method = 4)
print(Sys.time()-start)
bw.lpcde = lpcde(x_data = as.matrix(dat$x), y_data = as.matrix(dat$y), x = mean(dat$x))
bw.eebs = readRDS("bandwidth/data/dat1/bw_eebs.rds")
bw.eebs = EEBS_x(dat$x,dat$y,2,3,2,3,c(500),parallel = TRUE)


cde.out.1 = cde(x = dat$x, y = dat$y, x.margin = xmargin.short, y.margin = ymargin, 
                a = bw1$a, b = bw1$b, rescale = TRUE)
cde.out = cde(x = dat$x, y = dat$y, x.margin = xmargin.short, y.margin = ymargin, rescale = TRUE)
cde.out.3 = cde(x = dat$x, y = dat$y, x.margin = xmargin.short, y.margin = ymargin, 
                a = bw3$a, b = bw3$b, rescale = TRUE)
cde.out.4 = cde(x = dat$x, y = dat$y, x.margin = xmargin.short, y.margin = ymargin, 
                a = bw4$a, b = bw4$b, rescale = TRUE)
cde.out.opt.bw = cde(x = dat$x, y = dat$y, x.margin = xmargin.short, y.margin = ymargin, 
              a = bw$kernel[1], b = bw$kernel[2], rescale = TRUE)
cde.out.opt.xbw = cde.xbw(dat$x,dat$y,xmargin.short,ymargin,bw.x$kernel[seq(1,100,by=11),])
cde.out.ls = cde.ls(dat$x, dat$y, xmargin.short, ymargin, rescale = TRUE)
cde.out.eebs = cde.xbw(dat$x,dat$y,xmargin.short,ymargin,bw.eebs$bw)
visualize((t(cde.out.eebs$z-cd$z))^2, ybins = c(400,xmargin.short), xbins = c(200, ymargin))

isex = list()
isex[[1]] <- list(cd$x, rowSums((cde.out$z-cd$z)^2*(ymargin[2]-ymargin[1])))
isex[[1]] <- list(cd$x, rowSums((cde.out.opt.bw$z-cd$z)^2*(ymargin[2]-ymargin[1])))
isex[[2]] <- list(cd$x, rowSums((cde.out.opt.xbw$z-cd$z)^2*(ymargin[2]-ymargin[1])))
isex[[3]] <- list(cd$x, rowSums((cde.out.ls$z-cd$z)^2*(ymargin[2]-ymargin[1])))
isex[[3]] <- list(cd$x, rowSums((cde.out.eebs$z-cd$z)^2*(ymargin[2]-ymargin[1])))
isex[[6]] <- list(cd$x, rowSums((cde.out.1$z-cd$z)^2*(ymargin[2]-ymargin[1])))
isex[[7]] <- list(cd$x, rowSums((cde.out.3$z-cd$z)^2*(ymargin[2]-ymargin[1])))
isex[[8]] <- list(cd$x, rowSums((cde.out.4$z-cd$z)^2*(ymargin[2]-ymargin[1])))

plot_multiple(isex, legend = c(#"global normal reference rule", 
                               "global optimal bandwidth", "optimal bandwidth along x", 
                               "location-scale"#, "EEBS"
                               #"Bashtannyk-Hyndman algorithm",
                               #"Bashtannyk-Hyndman regression", "Bashtannyk-Hyndman bootstrap"
                               ),
              ylim=c(0,0.003),type="l",xlab = "x",ylab = "Integrated squared error conditioned on x")


K.out1 = estimate_K(t(cde.out$z),x = xmargin, y = ymargin, f = f1(xmargin), xbins = bins, ybins = bins)
K.out.opt.bw = estimate_K(t(cde.out.opt.bw$z),x = xmargin, y = ymargin, f = f1(xmargin), xbins = bins, ybins = bins)
K.out.opt.xbw = estimate_K(t(cde.out.opt.xbw$z),x = xmargin, y = ymargin, f = f1(xmargin), xbins = bins, ybins = bins)
K.out.ls = estimate_K(t(cde.out.ls$z),x = xmargin, y = ymargin, f = f1(xmargin), xbins = bins, ybins = bins)
K.out.1 = estimate_K(t(cde.out.1$z),x = xmargin, y = ymargin, f = f1(xmargin), xbins = bins, ybins = bins)
K.out.3 = estimate_K(t(cde.out.3$z),x = xmargin, y = ymargin, f = f1(xmargin), xbins = bins, ybins = bins)
K.out.4 = estimate_K(t(cde.out.4$z),x = xmargin, y = ymargin, f = f1(xmargin), xbins = bins, ybins = bins)


visualize(K.out,bins,bins,limits = c(0,0.3))
visualize(K.out.opt.bw,bins,bins,limits = c(0,0.3))
visualize(K.out.opt.xbw,bins,bins,limits = c(0,0.3))
visualize((K-K.out)^2,bins,bins, limits=c(0,0.012))
visualize((K-K.out.opt.bw)^2,bins,bins, limits=c(0,0.012))
visualize((K-K.out.opt.xbw)^2,bins,bins, limits=c(0,0.012))


isex_K = list()
isex_K[[1]] <- list(bins[-1], colSums((K-K.out)^2))
isex_K[[2]] <- list(bins[-1], colSums((K-K.out.opt.bw)^2))
isex_K[[3]] <- list(bins[-1], colSums((K-K.out.opt.xbw)^2))
isex_K[[4]] <- list(bins[-1], colSums((K-K.out.ls)^2))
isex_K[[5]] <- list(bins[-1], colSums((K-K.out.1)^2))
isex_K[[6]] <- list(bins[-1], colSums((K-K.out.3)^2))
isex_K[[7]] <- list(bins[-1], colSums((K-K.out.4)^2))
plot_multiple(isex_K, legend = c("global normal reference rule", "global optimal",
                                 "optimal along x", "location-scale",
                                 "Bashtannyk-Hyndman algorithm",
                                 "Bashtannyk-Hyndman regression", "Bashtannyk-Hyndman bootstrap"),
              ylim=c(0,0.05),type="b")




## conditional bias and variance estimation for EEBS ---------------------------

a = seq(0.5,10,by=0.5)
b = seq(0.5,10,by=0.5)
xlist = seq(400,1000,by=100)
ylist = seq(400,1000,by=100)
true.bias2.list = list()
true.variance.list = list()
bias2.eebs.list = list()
variance.eebs.list = list()


for (i in 1:length(xlist)) {
  for (j in 1:length(ylist)) {
    x0 = xlist[i]
    y0 = ylist[j]
    cat("(",x0,",",y0,") \n")
    result = true_bias2_variance(compute_cd(p = k1,x = c(x0,x0+1), y = c(y0,y0+1),
                                            rescale = FALSE)$z[1,1],a,b,x0,y0,f1,lb,ub,C1,C2,C3)
    true.bias2.list[[paste(x0,y0,sep = ",")]] = result[["bias2"]]
    true.variance.list[[paste(x0,y0,sep = ",")]] = result[["var"]]
    
    cde.eval <- cde_different_bw(dat$x, dat$y, a, b, c(x0,x0+1), c(y0,y0+1), FALSE)
    X <- cde.eval[cde.eval$x == x0 & cde.eval$y == y0,c('a','b','z')]
    data <- cbind(X[,c('a','b')], (X$a)*(X$b), (X$a)^2, (X$b)^2, X$z)
    colnames(data) <- c('s', 'b', 'ab', 'a2', 'b2', 'z')
    bias2.eebs.list[[paste(x0,y0,sep = ",")]] = bias_eebs(data)
    variance.eebs.list[[paste(x0,y0,sep = ",")]] = variance_eebs_boot(dat$x,dat$y,x0,y0,a,b)
  }
}
saveRDS(true.bias2.list, "bandwidth/data/dat1/true.bias2.list.rds")
saveRDS(true.variance.list, "bandwidth/data/dat1/true.variance.list.rds")
saveRDS(bias2.eebs.list, "bandwidth/data/dat1/bias2.eebs.list.rds")
saveRDS(variance.eebs.list, "bandwidth/data/dat1/variance.eebs.list.rds")



true.bias2.list = readRDS("bandwidth/data/dat1/true.bias2.list.rds")
true.variance.list = readRDS("bandwidth/data/dat1/true.variance.list.rds")
bias2.eebs.list = readRDS("bandwidth/data/dat1/bias2.eebs.list.rds")
variance.eebs.list = readRDS("bandwidth/data/dat1/variance.eebs.list.rds")
lb = min(min(bias2.eebs),min(true.bias2))
ub = max(max(bias2.eebs),max(true.bias2))
x0 = 600
y0 = 600
visualize(true.bias2.list[[paste(x0,y0,sep = ",")]], c(0,a), c(0,b), xlab="a", ylab="b", 
          legend_name = "true bias2", title = paste("x=",x0,", y=",y0,sep=""))
visualize(bias2.eebs.list[[paste(x0,y0,sep = ",")]], c(0,a), c(0,b), xlab="a", ylab="b", 
          legend_name = "EEBS bias2", title = paste("x=",x0,", y=",y0,sep=""))
visualize(true.variance.list[[paste(x0,y0,sep = ",")]], c(0,a), c(0,b), xlab="a", ylab="b", 
          legend_name = "true variance", title = paste("x=",x0,", y=",y0,sep=""))
visualize(variance.eebs.list[[paste(x0,y0,sep = ",")]], c(0,a), c(0,b), xlab="a", ylab="b", 
          legend_name = "EEBS variance", title = paste("x=",x0,", y=",y0,sep=""))
which(bias2.eebs==min(bias2.eebs), arr.ind = TRUE)
which(true.bias2==min(true.bias2), arr.ind = TRUE)


ymargin.short = seq(200,1200,length.out=20)
bias.eebs.x = list()
true.bias.x = list()
for (i in 1:length(ymargin.short)) {
  bias.eebs.x[[i]] = bias_eebs(dat$x,dat$y,a,b,x0,ymargin.short[i])
  true.bias.x[[i]] = true_bias(compute_cd(p = k1,x = c(x0,x0+1), y = c(ymargin.short[i],ymargin.short[i]+1))$z[1,1],
                                        dat$x,dat$y,a,b,x0,ymargin.short[i])
  cat(i, " ")
}
idx = 12
lb = min(min(bias.eebs.x[[idx]]),min(true.bias.x[[idx]]))
ub = max(max(bias.eebs.x[[idx]]),max(true.bias.x[[idx]]))
visualize(true.bias.x[[idx]], 1:20, 1:20, limits = c(0,0.25))
visualize(bias.eebs.x[[idx]], 1:20, 1:20, limits = c(0,1))


integ.bias.eebs.x = 0
integ.true.bias.x = 0
for (i in 1:length(ymargin.short)) {
  dy = ymargin.short[2] - ymargin.short[1]
  integ.bias.eebs.x = integ.bias.eebs.x + bias.eebs.x[[i]] * dy
  integ.true.bias.x = integ.true.bias.x + true.bias.x[[i]] * dy
}
lb = min(min(integ.bias.eebs.x),min(integ.true.bias.x))
ub = max(max(integ.bias.eebs.x),max(integ.true.bias.x))
visualize(integ.true.bias.x, 1:20, 1:20, c(0,250))
visualize(integ.bias.eebs.x, 1:20, 1:20, c(0,250))
which(integ.bias.eebs.x==min(integ.bias.eebs.x), arr.ind = TRUE)
which(integ.true.bias.x==min(integ.true.bias.x), arr.ind = TRUE)
