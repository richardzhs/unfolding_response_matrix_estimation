## -------------------------------------------------------
## Response Matrix Estimation for Simulated Jet pT.
##
## This R file contains the codes for applying the cde methods for estimating the
## response matrix to the simulated inclusive jet transverse momentum spectrum.
##
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
nxbins = 40
xbins = seq(lb, ub, length.out = nxbins + 1)
nybins = 40
ybins = seq(lb, ub, length.out = nybins + 1)

# evaluation points for cde
xmargin = seq(400, 1000, length.out = 1000)
ymargin = seq(200, 1200, length.out = 1000)


## Theoretical functions & quantities ------------------------------------------

# particle-level intensity function (pt in units of GeV)
f0 <- function(pt) {
  L*N0*(pt)^(-alpha)*(1-2/sqrt(s)*pt)^beta*exp(-gamma/pt)
}
plot(f0, xlim = c(lb, ub))

# mean vector lambda for particle-level histograms
lambda = rep(NA, nxbins)
for (i in 1:nxbins) {
  lambda[i] = integrate(f = f0, lower = xbins[i], upper = xbins[i+1])$value
}

# response kernel
k <- function(x,y) {
  dnorm(y, mean = x, sd = sqrt(x*C1^2+C2^2+x^2*C3^2))
}


## generate MC data -------------------------------------------------------------
dat_mc = generateData(100000, f0, lb, ub, C1, C2, C3)


## Estimate conditional densities and response matrices ------------------------
# location-scale
cde.out1 = cde_ls(dat_mc$x, dat_mc$y, xmargin, ymargin, rescale = TRUE)
# global kernel
cde.out2 = cde(x = dat_mc$x, y = dat_mc$y, x.margin = xmargin, y.margin = ymargin, rescale = TRUE)
# local kernel
cde.out3 = cde_local_x(dat_mc$x, dat_mc$y, xmargin, ymargin, d1 = 1, d2 = 1/400,
                       rescale = TRUE, parallel = FALSE)
# local linear
cde.out4 = cde(x = dat_mc$x, y = dat_mc$y, x.margin = xmargin, y.margin = ymargin, 
               rescale = TRUE, deg = 1, maxk = 1000)

# estimated MC particle-level density
f1 <- density(x = dat_mc$x, n = length(xmargin), from = xmargin[1], to = xmargin[length(xmargin)])$y


# estimated response matrices 
K1 <- estimate_K(k = t(cde.out1$z), x = cde.out1$x, y = cde.out1$y, f = f1, 
                 xbins = xbins, ybins = ybins)
K2 <- estimate_K(k = t(cde.out2$z), x = cde.out2$x, y = cde.out2$y, f = f1, 
                 xbins = xbins, ybins = ybins)
K3 <- estimate_K(k = t(cde.out3$z), x = cde.out3$x, y = cde.out3$y, f =  f1, 
                 xbins = xbins, ybins = ybins)
K4 <- estimate_K(k = t(cde.out4$z), x = cde.out4$x, y = cde.out4$y, f =  f1, 
                 xbins = xbins, ybins = ybins)
Knaive <- estimate_K_naive(dat_mc$x[dat_mc$x<xbins[length(xbins)] & dat_mc$y<ybins[length(ybins)]], 
                           dat_mc$y[dat_mc$x<xbins[length(xbins)] & dat_mc$y<ybins[length(ybins)]], xbins, ybins)

# true response matrix
K <- compute_K_normal(C1, C2, C3, 0, f1, xbins, ybins)
# true mu
mu <- K %*% lambda
# detector-level histogram
y <- generateDiscreteData(lambda, K)[[2]]

## Save the results ------------------------------------------------------------
#write.csv(K, "data/jet_pt/normal_kernel/response_matrices/Ktrue.csv",row.names = FALSE)
#write.csv(K1, "data/jet_pt/normal_kernel/response_matrices/K1.csv",row.names = FALSE)
#write.csv(K2, "data/jet_pt/normal_kernel/response_matrices/K2.csv",row.names = FALSE)
#write.csv(K3, "data/jet_pt/normal_kernel/response_matrices/K3.csv",row.names = FALSE)
#write.csv(K4, "data/jet_pt/normal_kernel/response_matrices/K4.csv",row.names = FALSE)
#write.csv(Knaive, "data/jet_pt/normal_kernel/response_matrices/Knaive.csv", row.names = FALSE)


#write.csv(lambda, file = "data/jet_pt/normal_kernel/true_means/x.csv", row.names = FALSE)
#write.csv(mu, file = "data/jet_pt/normal_kernel/smeared_means/y.csv", row.names = FALSE)
#write.csv(y, file = "data/jet_pt/normal_kernel/smeared_histograms/y.csv", row.names = FALSE)


## Visualize the estimated response matrices -----------------------------------
limits = c(-0.05,0.7)
ptrue <- visualize(K, xbins, ybins, limits = limits, title = "true response matrix")
p1 <- visualize(K1, xbins, ybins, limits = limits, title = "location-scale")
p2 <- visualize(K2, xbins, ybins, limits = limits, title = "kernel cde")
p3 <- visualize(K3, xbins, ybins, limits = limits, title = "local kernel cde")
p4 <- visualize(K4, xbins, ybins, limits = limits, title = "local linear cde")
pnaive <- visualize(Knaive, xbins, ybins, limits = limits, title = "histogram estimate")

p1