# SFDA-Horseshoe-Wavelets

[![R](https://img.shields.io/badge/R-%23276DC3.svg?style=flat&logo=r&logoColor=white)](https://www.r-project.org/)
[![Stan](https://img.shields.io/badge/Stan-%23B22B31.svg?style=flat&logo=stan&logoColor=white)](https://mc-stan.org/)

Official repository containing the code and data associated with the manuscript: **"Wavelet-Based Bayesian Hierarchical Modeling with Regularized Horseshoe Priors for Spatially Correlated Functional Data"**.

## Overview

This repository provides the implementation of a parsimonious Wavelet-Based Bayesian Hierarchical Model specifically tailored for **Spatial Functional Data Analysis (Spatial FDA)**. 

To overcome the severe parameter inflation that typically arises when modeling wavelet coefficients with spatial dependence, our methodology anchors a single, shared Matérn spatial correlation matrix at each resolution level. Adaptive sparsity is enforced through a spatially informed **Regularized Horseshoe (RHS)** prior. To ensure computational efficiency and avoid challenging posterior geometries when using the No-U-Turn Sampler (NUTS), the model is implemented using a strictly Non-Centered Parameterization (NCP).

## Repository Structure

The repository is organized into two main sections, reflecting the manuscript's core analyses: the real-world environmental application and the Monte Carlo simulation studies.

```text
SFDA-Horseshoe-Wavelets/
├── Application_PM10/              # Real-world application (Air pollution in Mexico City)
│   ├── Fitting_Code/              # Inference and reconstruction scripts
│   │   ├── Model_horseshoe_wavelets.R  # Fits the proposed spatial model
│   │   ├── Model_OKFD.R                # Fits the baseline model (Ordinary Kriging for Functional Data)
│   │   └── Curve_Reconstruction.R      # Reconstructs the signal in the original domain
│   └── Main_Data/                 # Georeferenced databases and time series
│       ├── PM10_2025_day.RData         # Imputed PM10 concentration data
│       ├── coord.RData                 # Geographic coordinates of the stations
│       └── Results/                    # Model outputs and results
│
└── Simulation_Studies/            # Validation via simulation studies (Noise Regimes)
    ├── Main_Data/                 # Synthetic data generation and storage
    │   ├── Data_Generation_Process.R   # Script to simulate the B-spline spatial process
    │   ├── SNR/                        # Data categorized by Signal-to-Noise Ratio (SNR)
    │   │   ├── High/                   # Scenario 1: High SNR (Low noise)
    │   │   ├── Medium/                 # Scenario 2: Moderate SNR
    │   │   └── Low/                    # Scenario 3: Low SNR (High noise)
    │   └── Ytrue.RData                 # True latent curves
    └── MCMC/                      # Directories for Stan posterior samples
        ├── Simulation_1/
        ├── Simulation_2/
        └── Simulation_3/

## Prerequisites and Dependencies
The model implementation is conducted in R, with HMC sampling performed via Stan. Please ensure your environment has the following packages installed:
