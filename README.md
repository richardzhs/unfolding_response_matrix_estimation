# Response Matrix Estimation in Unfolding Elementary Particle Spectra
This repository contains code for estimating response matrices and performing the corresponding unfolding tasks in particle physics.

## Overview
The repository includes the following scripts:

* `simulation_jets.R`: Estimates the response matrix for the inclusive jet transverse momentum spectrum.

* `Drell-Yan.R`: Estimates the response matrix for Drell–Yan events.

* `demo.ipynb`: Demonstrates the unfolding procedure using the estimated response matrices. This script reproduces the main results presented in the paper.

* `utils.R`: Utility functions used for estimating the response matrices.

* `cde.bandwidths.R`: Wrapper for bandwidth selection functions from [hdrcde] (https://github.com/robjhyndman/hdrcde), applied to kernel and local-linear conditional density estimation.

* `data.zip`: Contains simulated data as well as real data for Drell-Yan events.