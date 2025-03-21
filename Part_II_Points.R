# Packages -----------------------------------------------------------------

# Install only the packages that are not already installed
install.packages(setdiff(c("geobr", "dplyr", "sp", "sf", "spatstat", "ggspatial",
                           "ggrepel"), 
                         rownames(installed.packages())), dependencies = TRUE)

library(geobr)        # Obtaining IBGE shapefiles
library(dplyr)        # Data manipulation
library(sp)           # Spatial data analysis
library(sf)           # Spatial data manipulation
library(spatstat)     # Point pattern analysis
library(ggspatial)    # ggplot2 extension for spatial data
library(ggrepel)      # Avoids label overlap in ggplot2 plots

# Data Reading ------------------------------------------------------------

# Load wildfire data from a CSV file
df_queimadas <- read.csv("Data/focos_qmd_inpe_2025-01-01_2025-01-31_36.470290.csv")

# Display a statistical summary of the data
summary(df_queimadas)

# Display the total number of wildfire foci in the dataset
nrow(df_queimadas)

# State Mesh --------------------------------------------------------------

# Get the state mesh for Minas Gerais (year 2020)
MGbord <- read_state(year = 2020) |> filter(name_state == "Minas Gerais")

# Transform to the UTM coordinate system (Zone 22S)
MGbord_projetada <- st_transform(MGbord, crs = 32722)

# Convert coordinates to kilometers for better interpretation
MGbord_projetada_KM <- st_geometry(MGbord_projetada) / 1000

# Distribution ------------------------------------------------------------

# Select only the Longitude and Latitude columns of the wildfire foci
xy_vetor <- df_queimadas[,c("Longitude","Latitude")]
head(xy_vetor)  # Display the first few rows of coordinate data

# Define the projection system for geographic coordinates (WGS84)
proj <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# Create a SpatialPointsDataFrame object with the wildfire points
sp.points <- SpatialPointsDataFrame(coords=xy_vetor, data=df_queimadas, proj4string=proj)

# Get the bounding box of the spatial points
bbox(sp.points)

# Convert the state mesh to a spatial observation window
window <- as.owin(MGbord_projetada_KM)

# Convert the spatial points to sf (Simple Features) format
sp.points_sf <- st_as_sf(sp.points)

# Transform the points to the UTM coordinate system (Zone 22S)
sp.points_projetado <- st_transform(sp.points_sf, crs = 32722)

# Convert the point coordinates to kilometers
sp.points_projetado <- st_geometry(sp.points_projetado) / 1000

# Extract the coordinates of the projected points
coords <- st_coordinates(sp.points_projetado)

# Create a Planar Point Pattern (PPP) object
MG.ppp <- ppp(x=coords[,1], y=coords[,2], window=window)

# Plot the points on the region map
plot(MG.ppp, pch=16, cex=0.5, main=" ")


# Kolmogorov-Smirnov Intensity --------------------------------------------

# Estimating the average intensity of points (lambda)
lamb <- summary(MG.ppp)$intensity
lamb  # Display the average intensity

# Kolmogorov-Smirnov test to evaluate if the points follow a homogeneous Poisson process
KS <- cdf.test(MG.ppp, "x", test="ks")

# Plot the result of the Kolmogorov-Smirnov test
plot(KS, main= "")

# Get and display the p-value of the test
pval <- KS$p.value
pval
KS


# Intensity Maps ----------------------------------------------------------

# Estimate the density of points (smoothing of spatial intensity)
density <- density(MG.ppp)

# Plot the intensity map of wildfire foci
plot(density, xlab="Distance", col=terrain.colors(200),
     ylab="Distance (kilometers)", use.marks=TRUE,
     main="")

# Add the points to the intensity map
plot(MG.ppp, add=T, cex=0.1, use.marks=TRUE)

# Add contours to the intensity map
contour(density, main="Contour Map", add=TRUE)


# F Function --------------------------------------------------------------

# Define a distance vector for analysis of F, G, and K functions (up to 50 km)
r <- seq(0, 50, length=10000)

# Calculate envelopes for the F function
F <- envelope(MG.ppp, Fest, nsim = 500, r=r, correction="border")

# Plot the observed F function and theoretical envelopes
plot(F$r, F$theo, main = "", type="n",
     xlab="Distance (Kilometers)", ylab="F(x)")

# Add lines representing the different values of the F function
lines(F$r, F$theo, lty=3, col="blue")  # Theoretical F
lines(F$r, F$hi, lty=2, col="red")    # Upper envelope
lines(F$r, F$lo, lty=2, col="red")    # Lower envelope
lines(F$r, F$obs, lty=1, col="black") # Observed F

# Add a legend to the plot
legend("bottomright", legend=c("Theoretical F", "Upper F", "Lower F", "Observed F"),
       col=c("blue", "red", "red", "black"),
       lty=c(3, 2, 2, 1), cex = 0.7)
