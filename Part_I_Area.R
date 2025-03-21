# Packages -----------------------------------------------------------------

# Install only the packages that are not already installed
install.packages(setdiff(c("tidyverse", "dplyr", "geobr", "sf", "readxl", 
                           "stringr", "ggspatial", "spdep", "sp"), 
                         rownames(installed.packages())), dependencies = TRUE)

# Load the packages
library(tidyverse) # Data manipulation
library(dplyr) # Data manipulation
library(geobr) # IBGE shapefiles
library(sf) # Spatial data manipulation
library(readxl) # Reading Excel data
library(stringr) # Qualitative data manipulation
library(ggspatial) # ggplot2 extension for spatial data
library(spdep) # Spatial data manipulation
library(sp) # Spatial data manipulation


# Polygons ----------------------------------------------------------------

# Get the IBGE mesh / shapefile (Brazilian states)
shp <- geobr::read_municipality(year = 2022) # latest available year: 2022

# Graphical visualization - Brazil
shp |> ggplot() + geom_sf()

# Get the mesh for the state of Minas Gerais
shp_MG <- shp |> filter(name_state=="Minas Gerais")

# Graphical visualization - Minas Gerais
shp_MG |> ggplot() + geom_sf()

# Save the mesh in the project folder
shp_MG |> st_write("Malhas/MG_Estado_2022.shp", delete_layer = TRUE)

# Additional: Neighborhoods
read_neighborhood(year = 2022) |> 
  filter(name_muni=="Belo Horizonte") |> 
  ggplot() + geom_sf()

# It is also possible to obtain the meshes from the IBGE website:
# https://www.ibge.gov.br/geociencias/organizacao-do-territorio/malhas-territoriais/15774-malhas.html


# Data Reading ------------------------------------------------------------

df_dengue_MG <- read_excel("Data/Dados_Dengue_Pop.xlsx")
head(df_dengue_MG) # Check the structure of the dataset
summary(df_dengue_MG) # Descriptive summary of the data

# Preparation of Spatial Data ---------------------------------------------

ind <- match(shp_MG$name_muni, df_dengue_MG$Município); sum(is.na(ind))
shp_MG$name_muni[is.na(ind)] # Check the divergent municipalities

# Correction in the df_dengue_MG dataset
df_dengue_MG$Município <- str_replace(df_dengue_MG$Município, "Barão do Monte Alto", "Barão de Monte Alto")
df_dengue_MG$Município <- str_replace(df_dengue_MG$Município, "Olhos-d'Água", "Olhos-D'água")
df_dengue_MG$Município <- str_replace(df_dengue_MG$Município, "Pingo-d'Água", "Pingo-D'água")
df_dengue_MG$Município <- str_replace(df_dengue_MG$Município, "São João del Rei", "São João Del Rei")

ind <- match(shp_MG$name_muni, df_dengue_MG$Município); sum(is.na(ind))

# Reordering and adjusting the dataset
df_dengue_ind <- df_dengue_MG[ind,]
row.names(df_dengue_ind) <- shp_MG$name_muni

# Merge the dataframes
ShapeMG <- left_join(shp_MG, df_dengue_ind, by = c("name_muni" = "Município"))
head(ShapeMG)


# Quantile Map ------------------------------------------------------------

# Define the intervals
quantile(ShapeMG$Casos_mil_hab, seq(0.2,1,by=0.2)) # Calculate the quantiles
intervals <- c(-1, 24.424, 43.564, 70.510, 114.416, 319.930) # set the intervals (subtract 1 from the first and add 1 to the last)
legends <- c("Between 0 and 24.424", "Between 24.424 and 43.564", "Between 43.564 and 70.510", 
             "Between 70.510 and 114.416","Between 114.416 and 318.930")
colors <- c("Between 0 and 24.424"="white", "Between 24.424 and 43.564"="#FFC600", "Between 43.564 and 70.510"="#FF8D00", 
            "Between 70.510 and 114.416"="#FF3800", "Between 114.416 and 318.930"="darkred")

# Add the class column to ShapeMG
ShapeMG$quantis <- cut(ShapeMG$Casos_mil_hab, breaks = intervals, labels = legends)

# Map
ggplot(data = ShapeMG) +
  geom_sf(aes(fill = quantis), color = "darkgray", size = 0.2) +
  scale_fill_manual(values = colors, name = "Rate 1:1000") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = "white", color = "black"),
    axis.title = element_blank()  # Remove axis titles
  ) +
  guides(fill = guide_legend(override.aes = list(size = 1))) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.5, "cm"), pad_y = unit(0.5, "cm"),
                         style = north_arrow_fancy_orienteering)


# Neighborhood ------------------------------------------------------------

ShapeMG.nb <- poly2nb(ShapeMG, queen = T)
neighborhood <- nb2listw(ShapeMG.nb, style="C"); neighborhood # Symmetric weight matrix

plot(ShapeMG$geom, col = "lightgray")
plot(ShapeMG.nb, ShapeMG$geom, add = TRUE, col = "blue", lwd = 0.5)


# Global Moran's Index ----------------------------------------------------

Mglobal <- moran.test(ShapeMG$Casos_mil_hab, listw=neighborhood)
Mglobal


# Local Moran's Index -----------------------------------------------------

ShapeMG.mloc <- localmoran(ShapeMG$Casos_mil_hab, listw=neighborhood,
                           zero.policy=T, 
                           alternative = "two.sided")

head(ShapeMG.mloc)


# LISA Cluster ------------------------------------------------------------

# Calculate deviations
Sd_1 <- ShapeMG$Casos_mil_hab - mean(ShapeMG$Casos_mil_hab)
mI_1 <- ShapeMG.mloc[, 1]
# Determine the quadrants
quadrant <- vector(mode = "numeric", length = nrow(ShapeMG))
quadrant[Sd_1 >= 0 & mI_1 >= 0] <- 1
quadrant[Sd_1 <= 0 & mI_1 >= 0] <- 2
quadrant[Sd_1 >= 0 & mI_1 <= 0] <- 3
quadrant[Sd_1 <= 0 & mI_1 <= 0] <- 4
# Significance
signif <- 0.05
quadrant[ShapeMG.mloc[, 5] > signif] <- 5
# Add quadrants to the dataframe
ShapeMG$quadrant <- factor(quadrant, levels = 1:5, 
                           labels = c("high-high", "low-low", "high-low", 
                                      "low-high", "N.Sig"))
# Define the colors
colors <- c("high-high" = "red", "low-low" = "blue", "high-low" = "lightpink", 
            "low-high" = "skyblue2", "N.Sig" = "white")
# Create the plot with ggplot2
ggplot(data = ShapeMG) +
  geom_sf(aes(fill = quadrant), color = "gray", size = 0.2) +
  scale_fill_manual(values = colors, name = "Moran's Quadrants") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = "white", color = "black"),
    axis.title = element_blank()  # Remove axis titles
  ) +
  # Add the legend layer manually
  guides(fill = guide_legend(override.aes = list(size = 1))) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.5, "cm"), pad_y = unit(0.5, "cm"),
                         style = north_arrow_fancy_orienteering)

# LISA Significance -------------------------------------------------------

# Define the LISA significance intervals and colors
intervals <- c(0, 0.001, 0.01, 0.05, 1)
legends <- c("0.1%", "1.0%", "5.0%", "N.Sig")
colors <- c("0.1%"="red", "1.0%"="blue", "5.0%"="skyblue", "N.Sig"="white")

# Add the significance class column to ShapeMG
ShapeMG$signif_class <- cut(ShapeMG.mloc[,5], breaks = intervals, labels = legends, right = TRUE)

# Create the plot with ggplot2
ggplot(data = ShapeMG) +
  geom_sf(aes(fill = signif_class), color = "darkgray", size = 0.2) +
  scale_fill_manual(values = colors, name = "Significance") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = "white", color = "black"),
    axis.title = element_blank()  # Remove axis titles
  ) +
  guides(fill = guide_legend(override.aes = list(size = 1))) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.5, "cm"), pad_y = unit(0.5, "cm"),
                         style = north_arrow_fancy_orienteering)
