# Packages -----------------------------------------------------------------

library(geoR)         # Geostatistical analysis
library(geobr)        # Obtaining IBGE shapefiles
library(dplyr)        # Data manipulation
library(MASS)         # Boxcox function (Box-Cox transformation)
library(sf)           # Spatial data manipulation
library(ggplot2)      # Plot creation

# Data Reading ------------------------------------------------------------

# Load precipitation data from a CSV file
df_precipitacao <- read.csv("Data/focos_qmd_inpe_2025-01-01_2025-01-31_36.470290.csv")

# Convert the data to the 'geodata' format of the geoR package,
# specifying the coordinate columns (columns 11 and 10) and the data column (column 8)
tb1 <- as.geodata(df_precipitacao, coords.col = c(11,10), data.col = 8)
summary(tb1)        # Display a summary of the data
plot(tb1, low = TRUE)  # Plot the data with the 'low' option to visualize low-density points

# Adjust data to avoid issues with logarithm of zero
tb1$data <- tb1$data + 0.01  
boxcox(tb1)         # Apply Box-Cox transformation (to check the best transformation parameter)

# Logarithmic transformation of precipitation data
df_precipitacao$Precipitacao <- log10(df_precipitacao$Precipitacao + 0.01)
# Reconvert the transformed data to the 'geodata' object
tb1 <- as.geodata(df_precipitacao, coords.col = c(11,10), data.col = 8)
summary(tb1)        # Display summary after transformation

# Plot the data again, with different trends
plot(tb1, low = TRUE)
plot(tb1, low = TRUE, trend = "1st")  # Consider linear trend (1st order)
plot(tb1, low = TRUE, trend = "2nd")  # Consider quadratic trend (2nd order)


# Definition of borders for the prediction area ----------------------------

# Read the shapefile for the state of Minas Gerais (year 2020)
bor <- read_state(year = 2020) |> filter(name_state == "Minas Gerais")
# Transform the polygon into a data frame with coordinates (longitude and latitude)
bor_coords <- st_coordinates(bor)       # Extract the polygon coordinates
bor_df <- as.data.frame(bor_coords[, 1:2]) # Select only X (Longitude) and Y (Latitude)
colnames(bor_df) <- c("Longitude", "Latitude")

# Add the borders to the 'tb1' object
tb1$borders <- with(bor_df, cbind(Longitude, Latitude))


# Variogram Calculation ---------------------------------------------------

# Calculate the empirical variogram with 1st order trend and maximum distance of 5
grama <- variog(tb1, max.dist = 5, trend = "1st")
# Generate a variogram envelope (simulation) for uncertainty assessment with 100 simulations
env.var <- variog.mc.env(tb1, obj.v = grama, nsim = 100)
ef1 = eyefit(grama)  # Code commented for visual variogram fitting
summary(ef1)
# Plot the empirical variogram along with the simulation envelope
plot(grama, env = env.var)


# Variogram Model Fitting using Maximum Likelihood ------------------------

# Initial values and fixed nugget (sigma2, phi, and tau2)
# sigma2 = 1.31; phi = 1.6; tau2 = 0.16

# Create a list to store the fits of different models
tb1.ml0 <- list()
# Fit a model with 'matern' covariance
tb1.ml0$fit1 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "matern", trend = "1st")
# Fit a model with 'exponential' covariance
tb1.ml0$fit2 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "exponential", trend = "1st")
# Fit a model with 'gaussian' covariance
tb1.ml0$fit3 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "gaussian", trend = "1st")
# Fit a model with 'spherical' covariance
tb1.ml0$fit4 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "spherical", trend = "1st")
# Fit a model with 'circular' covariance
tb1.ml0$fit5 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "circular", trend = "1st")
# Fit a model with 'cubic' covariance
tb1.ml0$fit6 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "cubic", trend = "1st")
# Fit a model with 'wave' covariance
tb1.ml0$fit7 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "wave", trend = "1st")
#tb1.ml0$fit8 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "power", trend = "1st") # Code commented
# Fit a model with 'powered.exponential' covariance
tb1.ml0$fit9 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "powered.exponential", trend = "1st")
# Fit a model with 'cauchy' covariance
tb1.ml0$fit10 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "cauchy", trend = "1st")
#tb1.ml0$fit11 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "gencauchy", trend = "1st") # Code commented
# Fit a model with 'gneiting' covariance
tb1.ml0$fit12 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "gneiting", trend = "1st")
#tb1.ml0$fit13 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "gneiting.matern", trend = "1st") # Code commented
# Fit a model with 'pure.nugget' covariance (only nugget, no spatial structure)
tb1.ml0$fit14 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "pure.nugget", trend = "1st")

# Compare models using the Akaike Information Criterion (AIC)
sapply(tb1.ml0, AIC)

# Display the summary of the model with 'cauchy' covariance
cauchymodel <- summary(likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "cauchy", trend = "1st"))
cauchymodel

# Plot the data points divided into quintiles with legends
points(tb1, pt.div = "quintile", xlab = "Coord X", 
       ylab = "Coord Y", x.leg = -15.8, y.leg = -48, trend = "1st")

# Add legend

# Create a prediction grid within the defined borders ---------------------

gr <- pred_grid(tb1$borders, by = .1)   # Create a grid with spacing of 0.1
gr0 <- gr[.geoR_inout(gr, tb1$borders), ] # Select points that are within the borders


# Kriging ----------------------------------------------------------------

# Define kriging control using the previously fitted 'cauchy' model
KC <- krige.control(type = "SK", obj.model = tb1.ml0$fit10)
# Perform simple kriging on the grid points
tb1.kc <- krige.conv(tb1, loc = gr0, krige = KC)


# Mapping Kriging Results ------------------------------------------------

# Image of the interpolated value (using the 'tb1.kc' object)
image(tb1.kc, gr,
      tb1$borders,
      col = terrain.colors(200),
      x.leg = c(-45, -40.8), y.leg = c(-24, -23.5), cex = 0.7, trend = "1st")

# Image of kriging variance (uncertainty)
image(tb1.kc, loc = gr,
      bor = tb1$borders, col = terrain.colors(21),
      val = sqrt(tb1.kc$krige.var),
      x.leg = c(-50.5, -40.8), y.leg = c(-24, -23.5), cex = 0.6)


# Conditional Kriging Output ---------------------------------------------

# Define output controls with threshold and quantiles for posterior analysis
OC <- output.control(thres = median(tb1$data),
                     quantile = c(0.1, 0.9))

# Perform kriging with output control
tb1.kc <- krige.conv(tb1, loc = gr0, krige = KC, out = OC)

# Image of conditional probability (1 - probability)
image(tb1.kc, loc = gr, bor = tb1$borders, col = terrain.colors(21), 
      val = 1 - tb1.kc$prob, x.leg = c(-50.5, -40.8), y.leg = c(-24, -23.5), cex = 0.7)

# Image of the 90% quantile of kriging
image(tb1.kc, loc = gr, bor = tb1$borders, col = terrain.colors(21), 
      val = tb1.kc$quantile$q90, x.leg = c(-50.5, -40.8), y.leg = c(-24, -23.5), cex = 0.7)

# Image of the 10% quantile of kriging
image(tb1.kc, loc = gr, bor = tb1$borders, col = terrain.colors(21), 
      val = tb1.kc$quantile$q10, x.leg = c(-50.5, -40.8), y.leg = c(-24, -23.5), cex = 0.7)


# Simulations  ------------------------------------------------------------

# Calculate the probability that simulations exceed the median of the data
p091 <- apply(tb1.kc$simul, 2, function(x) mean(x > median(tb1$data)))
# Plot a histogram of the obtained probabilities with overlaid density
hist(p091, prob = TRUE)
lines(density(p091))