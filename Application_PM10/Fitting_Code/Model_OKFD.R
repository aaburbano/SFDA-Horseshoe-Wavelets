rm(list=ls(all=TRUE))
# Load the libraries
library("fdagstat")
library("fda")
library("ggplot2")
library("dplyr")
library("tidyr")
library("sf")
library("gstat")

#-------------------------------------------------------------------------------
# Set directories and load data
#-------------------------------------------------------------------------------
load("../Main_Data/PM10_2025_day.RData")
load("../Main_Data/coordinates_utm_full_km.RData")
load("../Main_Data/coordinates_utm_obs_km.RData")
#-------------------------------------------------------------------------------

# Pre-processing: Apply log transformation
Sample = PM10_2025_day %>%
         select(-Day) %>%
         as.matrix()

Sample = apply(Sample, c(1, 2), function(x) ifelse(is.na(x), NA, log(x)))
Sample = as.data.frame(Sample)

# Select observed stations for modeling (excluding prediction sites)
PM102025_log = Sample %>%
                dplyr::select(-ACO, -IZT, -TAH)

PM102025_log = as.matrix(PM102025_log)

# Define smoothing parameters
nbasis = 72           # Number of basis functions
argvals = 1:128       # Time points (days)
rangeval = range(argvals)
lambda = 0            # Smoothing penalty (tuning parameter)

# Define prediction coordinates (UTM in km)
coordinates_pred_utm_km = coordinates_utm_full_km[18:20,]

# Setup Basis and Functional Parameter objects
bspline_basis = create.bspline.basis(rangeval = rangeval, nbasis = nbasis, norder = 4)
fdPar_obj = fdPar(fdobj = bspline_basis, Lfdobj = 2, lambda = lambda)

## Perform Functional Smoothing
smoothed_fd = smooth.basis(argvals = argvals, y = PM102025_log, fdParobj = fdPar_obj)
smoothed_values = eval.fd(argvals, smoothed_fd$fd)

# Visualize smoothed training curves
matplot(argvals, smoothed_values, type = 'l', lty = 1,
        xlab = "Day of the Year", ylab = "Smoothed Value",
        main = "Smoothed Curves (17 Observed Stations)",
        col = rainbow(17))

#-------------------------------------------------------------------------------
# Ordinary Kriging for Functional Data (OKFD) using fdagstat
#-------------------------------------------------------------------------------

# 1. Define the fstat object
g = fstat(g = NULL,  vName = "PM10", 
          Coordinates = coordinates_utm_obs_km, 
          Functions = data.frame(smoothed_values))

# 2. Define Drift/Trend for Ordinary Kriging ("~1" assumes a constant unknown mean)
g = estimateDrift("~1", g, Intercept = TRUE)

# 3. Calculate the maximum observed distance for variogram scaling
max_dist = max(dist(coordinates_utm_obs_km))

# 4. Compute the empirical trace-variogram
# "~Easting+Northing" allows the package to calculate spatial distances
g = fvariogram("~Easting+Northing", g, Nlags = 15, LagMax = max_dist, ArgStep = 1, comments = TRUE)

# 5. Plot the empirical variogram
plotVariogram(g)

# 6. Theoretical variogram model fitting
# Set initial values to assist solver convergence
modelo_inicial = vgm(psill = 1, model = "Sph", range = max_dist / 2, nugget = 0.1)
g = fitVariograms(g, modelo_inicial, fitRanges = TRUE, forceNugget = TRUE)

# 7. Add covariance structure based on fitted model
g = addCovariance(g)

# 8. Generate OKFD spatial predictions for new locations
forecasts_okf = predictFstat(g, .newCoordinates = coordinates_pred_utm_km, .what = "PM10", .type = "OK")

# Display final prediction results
forecasts_okf$Forecast

# Organize predicted curves into a data frame
Ypred_OKF = data.frame(ok1 = forecasts_okf$Forecast[,1], 
                       ok2 = forecasts_okf$Forecast[,2], 
                       ok3 = forecasts_okf$Forecast[,3])

# Visualize the 3 predicted curves
matplot(1:128, Ypred_OKF, type = 'l', lty = 1,
        xlab = "Day of the Year", ylab = "Smoothed Value",
        main = "Predicted Curves for Ungauged Sites (OKFD)")
