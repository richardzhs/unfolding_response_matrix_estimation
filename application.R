## -------------------------------------------------------
## Application to jet transverse momentum spectrum in Drell-Yan events from the
## LHC data
##
## This R file contains the codes for applying the cde methods on estimating the
## response matrix to the real data from the LHC
##

source("cde/bandwidth/bandwidth.selection.R")
source("cde/bandwidth/cde.bandwidths.R")
source("cde/cde.est.R")
source("utils.R")

## Helper functions -----------------------------------------------------------
#' Whether the values in the vector is NA
is_NA <- function(x) {
  x == "NA0" | x == "NA1" | x == "NA2" | x == "NA3" | x == "NA4"
}

#' Make the NAs in the dataset to be the same.
same_NA <- function(dat) {
  dat_NA = dat
  dat_NA$PT_x[is_NA(dat_NA$PT_x)] = NA
  dat_NA$PT_y[is_NA(dat_NA$PT_y)] = NA
  dat_NA$PT_x = as.numeric(as.character(dat_NA$PT_x))
  dat_NA$PT_y = as.numeric(as.character(dat_NA$PT_y))
  return(dat_NA)
}

#' Clean the NA in the dataset
clean_NA <- function(dat) {
  dat_clean = dat[!(is_NA(dat$PT_x) | is_NA(dat$PT_y)),]
  dat_clean[,c("Event","PT_x","PT_y")] = apply(dat_clean[,c("Event","PT_x","PT_y")], MARGIN = 2, FUN = as.character)
  dat_clean[,c("Event","PT_x","PT_y")] = data.frame(apply(dat_clean[,c("Event","PT_x","PT_y")], MARGIN = 2, FUN = as.numeric))
  return(dat_clean)
}


## Read Data ------------------------------------------------------------------
#jet <- read.delim("jet.txt", sep = "")
# jet data after removing NAs (i.e. missing events and fake events)
#jet_clean = clean_NA(jet)
jet_clean = read.csv("jet_clean.csv")
# for reproducibility
set.seed(100)
mc_sample <- sample(1:nrow(jet_clean),100000,replace = FALSE)
dat_mc <- data.frame(x = jet_clean$x[mc_sample], y = jet_clean$y[mc_sample])
dat_expr <- data.frame(x = jet_clean$x[-mc_sample], y = jet_clean$y[-mc_sample])

xmargin = seq(20, 400, length.out = 1000)
short.idx = seq(1,1000,by=10)
xmargin.short = xmargin[short.idx]
ymargin = seq(20, 400, length.out = 1000)

xbins = seq(20, 400, length.out = 41)
ybins = seq(20, 400, length.out = 41)

## Estimate conditional densities and response matrices ------------------------
# CDE
cde.out1 = cde_ls(dat_mc$x, dat_mc$y, xmargin, ymargin, rescale = TRUE)
cde.out2 = hdrcde::cde(x = dat_mc$x, y = dat_mc$y, x.margin = xmargin, y.margin = ymargin, rescale = TRUE)
cde.out3 = cde_local_x(dat_mc$x, dat_mc$y, xmargin, ymargin, d1 = 1, d2 = 1/70,
                       rescale = TRUE, parallel = FALSE)
cde.out4 = cde(x = dat_mc$x, y = dat_mc$y, x.margin = xmargin, y.margin = ymargin, 
                       rescale = TRUE, deg = 1, maxk = 1000)

# estimated particle-level density ansatz
f1 <- density(x = dat_mc$x, n = length(xmargin), from = xmargin[1], to = xmargin[length(xmargin)])$y
f1.short <- density(x = dat$x, n = length(xmargin.short), from = xmargin.short[1], 
                    to = xmargin.short[length(xmargin.short)])$y

# estimated response matrices
K1 <- estimate_K(k = t(cde.out1$z), x = cde.out1$x, y = cde.out1$y, f = f1, 
                 xbins = xbins, ybins = ybins)
K2 <- estimate_K(k = t(cde.out2$z), x = cde.out2$x, y = cde.out2$y, f = f1, 
                 xbins = xbins, ybins = ybins)
K3 <- estimate_K(k = t(cde.out3$z), x = cde.out3$x, y = cde.out3$y, f =  f1, 
                 xbins = xbins, ybins = ybins)
K4 <- estimate_K(k = t(cde.out4$z), x = cde.out4$x, y = cde.out4$y, f =  f1, 
                 xbins = xbins, ybins = ybins)
Knaive <- estimate_K_naive(dat_mc$x[dat_mc$x<xbins[length(xbins)]&dat_mc$y<ybins[length(ybins)]], 
                           dat_mc$y[dat_mc$x<xbins[length(xbins)]&dat_mc$y<ybins[length(ybins)]], xbins, ybins)

# experimental histograms
x_expr <- hist(dat_expr$x[dat_expr$x<xbins[length(xbins)]&dat_expr$y<ybins[length(ybins)]], breaks = xbins)$count
y_expr <- hist(dat_expr$y[dat_expr$x<xbins[length(xbins)]&dat_expr$y<ybins[length(ybins)]], breaks = xbins)$count

## Save the results ------------------------------------------------------------
write.csv(K1, strcat(c("../unfolding-main/data/Drell-Yan/response_matrices/K1", ".csv")),row.names = FALSE)
write.csv(K2, strcat(c("../unfolding-main/data/Drell-Yan/response_matrices/K2", ".csv")),row.names = FALSE)
write.csv(K3, strcat(c("../unfolding-main/data/Drell-Yan/response_matrices/K3", ".csv")),row.names = FALSE)
write.csv(K4, strcat(c("../unfolding-main/data/Drell-Yan/response_matrices/K4", ".csv")),row.names = FALSE)
write.csv(Knaive, strcat(c("../unfolding-main/data/Drell-Yan/response_matrices/Knaive", ".csv")),row.names = FALSE)

write.csv(x_expr, strcat(c("../unfolding-main/data/Drell-Yan/true_histograms/x", ".csv")))
write.csv(y_expr, strcat(c("../unfolding-main/data/Drell-Yan/smeared_histograms/y", ".csv")))

## Visualize the estimated response matrices -----------------------------------
# read the matrices
K1 <- as.matrix(read.csv("../unfolding-main/data/Drell-Yan/response_matrices/K1.csv"))
K2 <- as.matrix(read.csv("../unfolding-main/data/Drell-Yan/response_matrices/K2.csv"))
K3 <- as.matrix(read.csv("../unfolding-main/data/Drell-Yan/response_matrices/K3.csv"))
K4 <- as.matrix(read.csv("../unfolding-main/data/Drell-Yan/response_matrices/K4.csv"))
Knaive <- as.matrix(read.csv("../unfolding-main/data/Drell-Yan/response_matrices/Knaive.csv"))


limits = c(-0.05,0.7)
p1 <- visualize(K1, xbins, ybins, limits = limits, title = "location-scale")
p2 <- visualize(K2, xbins, ybins, limits = limits, title = "kernel cde")
p3 <- visualize(K3, xbins, ybins, limits = limits, title = "local kernel cde")
p4 <- visualize(K4, xbins, ybins, limits = limits, title = "local linear cde")
pnaive <- visualize(Knaive, xbins, ybins, limits = limits, title = "histogram estimate")

p1
p2
p3
p4
pnaive


grid.arrange(pnaive, p2, p4, p3, p1, ncol=2,
             top = textGrob("Response Matrix Estimation for the simuated Drell-Yan Events (K=40x40, n=100000)",
                            gp=gpar(fontsize=20,font=3)))



## Sub-sample data and compute the response matrices/MSE ------------------------
# use the histogram estimate using all the data as the true matrix
Ktrue <- estimate_K_naive(dat$x[dat$x<300&dat$y<300], 
                          dat$y[dat$x<300&dat$y<300], xbins, ybins)
write.csv(Ktrue, strcat(c("unfolding_osb_po_uq-main/data/jet_pt/config79/response_matrices/Ktrue", ".csv")),row.names = FALSE)

registerDoParallel(50)
start.time = Sys.time()
foreach (i=1:100) %dopar% {
  idx <- sample(1:nrow(dat), size = 100000, replace = FALSE)
  sub_dat <- dat[idx,]
  #cde.out1 = cde_ls(sub_dat$x, sub_dat$y, xmargin, ymargin, rescale = TRUE)
  #cde.out2 = hdrcde::cde(x = sub_dat$x, y = sub_dat$y, x.margin = xmargin, y.margin = ymargin, rescale = TRUE)
  #cde.out3 = cde_local_x(sub_dat$x, sub_dat$y, xmargin, ymargin, d1 = 1, d2 = 1/70,
  #                       rescale = TRUE, parallel = FALSE)
  cde.out4 = hdrcde::cde(x = sub_dat$x, y = sub_dat$y, x.margin = xmargin, y.margin = ymargin, 
                         rescale = TRUE, deg = 1, maxk = 1000)
  
  # estimated particle-level density ansatz
  f1 <- density(x = sub_dat$x, n = length(xmargin), from = xmargin[1], 
                to = xmargin[length(xmargin)])$y
  
  # estimated response matrices
  #K1 <- estimate_K(k = t(cde.out1$z), x = cde.out1$x, y = cde.out1$y, f = f1, 
  #                 xbins = xbins, ybins = ybins)
  #K2 <- estimate_K(k = t(cde.out2$z), x = cde.out2$x, y = cde.out2$y, f = f1, 
  #                 xbins = xbins, ybins = ybins)
  #K3 <- estimate_K(k = t(cde.out3$z), x = cde.out3$x, y = cde.out3$y, f =  f1, 
  #                 xbins = xbins, ybins = ybins)
  K4 <- estimate_K(k = t(cde.out4$z), x = cde.out4$x, y = cde.out4$y, f =  f1, 
                   xbins = xbins, ybins = ybins)
  #Knaive <- estimate_K_naive(sub_dat$x[sub_dat$x<300&sub_dat$y<300], 
  #                           sub_dat$y[sub_dat$x<300&sub_dat$y<300], xbins, ybins)
  
  #write.csv(K1, pracma::strcat(c("unfolding_osb_po_uq-main/data/jet_pt/config79/response_matrices/K1-", i, ".csv")),row.names = FALSE)
  #write.csv(K2, pracma::strcat(c("unfolding_osb_po_uq-main/data/jet_pt/config79/response_matrices/K2-", i, ".csv")),row.names = FALSE)
  #write.csv(K3, pracma::strcat(c("unfolding_osb_po_uq-main/data/jet_pt/config79/response_matrices/K3-", i, ".csv")),row.names = FALSE)
  write.csv(K4, strcat(c("resp_mat/K4-", i, ".csv")),row.names = FALSE)
  #write.csv(Knaive, pracma::strcat(c("unfolding_osb_po_uq-main/data/jet_pt/config79/response_matrices/Knaive-", i, ".csv")),row.names = FALSE)
  cat(i, " ")
}
end.time = Sys.time()
cat("spent:", end.time-start.time)


# calculate the MAE
K <- read.csv("unfolding_osb_po_uq-main/data/jet_pt/config79/response_matrices/Ktrue.csv")
mat_mae <- array(0, c(nrow(K), ncol(K), 4))
niter <- 100
for (i in 1:niter) {
  K1 <- read.csv(strcat(c("unfolding_osb_po_uq-main/data/jet_pt/config79/response_matrices/K1-", i, ".csv")))
  K2 <- read.csv(strcat(c("unfolding_osb_po_uq-main/data/jet_pt/config79/response_matrices/K2-", i, ".csv")))
  K3 <- read.csv(strcat(c("unfolding_osb_po_uq-main/data/jet_pt/config79/response_matrices/K3-", i, ".csv")))
  #K4 <- read.csv(strcat(c("unfolding_osb_po_uq-main/data/jet_pt/config79/response_matrices/K4-", i, ".csv")))
  Knaive <- read.csv(strcat(c("unfolding_osb_po_uq-main/data/jet_pt/config79/response_matrices/Knaive-", i, ".csv")))
  mat_mae[,,1] <- mat_mae[,,1] + 1/niter * as.matrix(abs(K1-K))
  mat_mae[,,2] <- mat_mae[,,2] + 1/niter * as.matrix(abs(K2-K))
  mat_mae[,,3] <- mat_mae[,,3] + 1/niter * as.matrix(abs(K3-K))
  #mat_mae[,,4] <- mat_mae[,,4] + 1/niter * as.matrix(abs(K4-K))
  mat_mae[,,4] <- mat_mae[,,4] + 1/niter * as.matrix(abs(Knaive-K))
  cat(i, " ")
}

limits = c(0, 0.25)
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
pnaive <- visualize(mat_mae[,,4],xbins=xbins,ybins=ybins, limits=limits,
                title = "histogram estimate", legend_name = "mae", low_col = low, high_col = high)


grid.arrange(p1, p2, p3, pnaive, ncol=2,
             top = textGrob("MAE For Response Matrix Estimation (Simulated Drell-Yan Events, K=40x40, n=100000)",
                            gp=gpar(fontsize=20,font=3)))
