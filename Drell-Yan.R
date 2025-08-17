## -------------------------------------------------------
## Application to jet transverse momentum spectrum in Drell-Yan events + Jets
##
## This R file contains the codes for applying the cde methods for estimating the
## response matrix to the jet transverse momentum spectrum in Drell-Yan event.
##
source("utils.R")

## Read jet data ---------------------------------------------------------------
jet_clean = read.csv("data/Drell-Yan/jet_clean.csv")
# for reproducibility
set.seed(100)
mc_sample <- sample(1:nrow(jet_clean),100000,replace = FALSE)
dat_mc <- data.frame(x = jet_clean$x[mc_sample], y = jet_clean$y[mc_sample])
dat_expr <- data.frame(x = jet_clean$x[-mc_sample], y = jet_clean$y[-mc_sample])

xmargin = seq(20, 400, length.out = 1000)
ymargin = seq(20, 400, length.out = 1000)

xbins = seq(20, 400, length.out = 41)
ybins = seq(20, 400, length.out = 41)


## Estimate conditional densities and response matrices ------------------------
# location-scale
cde.out1 = cde_ls(dat_mc$x, dat_mc$y, xmargin, ymargin, rescale = TRUE)
# global kernel
cde.out2 = cde(x = dat_mc$x, y = dat_mc$y, x.margin = xmargin, y.margin = ymargin, rescale = TRUE)
# local kernel
cde.out3 = cde_local_x(dat_mc$x, dat_mc$y, xmargin, ymargin, d1 = 1, d2 = 1/100,
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
Knaive_full <- estimate_K_naive(jet_clean$x[jet_clean$x<xbins[length(xbins)] & jet_clean$y<ybins[length(ybins)]], 
                                jet_clean$y[jet_clean$x<xbins[length(xbins)] & jet_clean$y<ybins[length(ybins)]], xbins, ybins)

# experimental histograms
x_expr <- hist(dat_expr$x[dat_expr$x<xbins[length(xbins)]&dat_expr$y<ybins[length(ybins)]], breaks = xbins)$count
y_expr <- hist(dat_expr$y[dat_expr$x<xbins[length(xbins)]&dat_expr$y<ybins[length(ybins)]], breaks = xbins)$count


## Save the results ------------------------------------------------------------
#write.csv(K1, "data/Drell-Yan/response_matrices/K1.csv",row.names = FALSE)
#write.csv(K2, "data/Drell-Yan/response_matrices/K2.csv",row.names = FALSE)
#write.csv(K3, "data/Drell-Yan/response_matrices/K3.csv",row.names = FALSE)
#write.csv(K4, "data/Drell-Yan/response_matrices/K4.csv",row.names = FALSE)
#write.csv(Knaive, "data/Drell-Yan/response_matrices/Knaive.csv", row.names = FALSE)
#write.csv(Knaive_full, "data/Drell-Yan/response_matrices/Knaive_full.csv",row.names = FALSE)

#write.csv(x_expr, "data/Drell-Yan/true_histograms/x.csv")
#write.csv(y_expr, "data/Drell-Yan/smeared_histograms/y.csv")


## Visualize the estimated response matrices -----------------------------------
limits = c(-0.05,0.7)
p1 <- visualize(K1, xbins, ybins, limits = limits, title = "location-scale")
p2 <- visualize(K2, xbins, ybins, limits = limits, title = "kernel cde")
p3 <- visualize(K3, xbins, ybins, limits = limits, title = "local kernel cde")
p4 <- visualize(K4, xbins, ybins, limits = limits, title = "local linear cde")
pnaive <- visualize(Knaive, xbins, ybins, limits = limits, title = "histogram estimate")
pnaive_full <- visualize(Knaive_full, xbins, ybins, limits = limits, title = "histogram estimate (full sample)")

p1