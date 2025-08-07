library(ggplot2)
library(gridExtra)
library(MASS)
library(dplyr)
library(pbapply)
library(progress)
library(progressr)
library(sf)
library(writexl)
library(readxl)
library(terra)
library(raster)
library(tmap)
library(tidyverse)   # For data manipulation
library(exactextractr) # For efficient zonal stats (optional)
library(mapview)
library(fitdistrplus)
library(parallel)
library(arcgisbinding)


############################################################################################
### Calculating runoff from original vector data using Curve Number (CN) method

# Load required library
library(sf)

# Define the path to your local geodatabase (modify as needed)
geodatabase_path <- "path/to/your/Geodatabase.gdb"

# List all layers in the geodatabase to verify the layer names
layers <- st_layers(geodatabase_path)
print(layers)

# Read the CN_Socio vector data from the geodatabase
CN <- st_read(geodatabase_path, layer = "CN_Socio")

# Inspect the data
class(CN)
names(CN)
head(CN)
summary(CN)

# Set the geometry column explicitly
st_geometry(CN) <- "Shape"

# Initialize new columns for precipitation (P) and runoff (Qmm)
CN$P <- NA
CN$Qmm <- NA

# Sample precipitation values from a Monte Carlo simulation (user must define `simulated_precipitation`)
set.seed(123)  # For reproducibility
num_polygons <- nrow(CN)
CN$P <- sample(simulated_precipitation, num_polygons, replace = TRUE)

# Compare sampled precipitation to the simulation
summary(CN$P)
summary(simulated_precipitation)

# Calculate runoff (Qmm) using the SCS-CN method
# Assumes CN$S column already exists and is properly populated
CN$Qmm <- ifelse(
  CN$P < 0.2 * CN$S,
  0,
  (CN$P - 0.2 * CN$S)^2 / (CN$P + 0.8 * CN$S)
)

# Count and print the number of polygons with runoff > 0
num_polygons_with_runoff <- sum(CN$Qmm > 0)
cat("Number of polygons with Qmm > 0:", num_polygons_with_runoff, "\n")

# Summary of runoff results
summary(CN$Qmm)

# Save the updated vector data back to the geodatabase
# Replace the path below with your preferred output location
st_write(CN, dsn = "path/to/your/output/Geodatabase.gdb", layer = "CN_Socio", append = FALSE)



#############################################################################################
### Estimating Runoff from Curve Number (CN) Raster and Simulated Precipitation
### Based on the SCS-CN method using gamma-distributed precipitation

library(terra)

# Set working directory containing CN raster (optional: adjust for reproducibility)
# setwd("path_to_your_raster_directory")  

# Load the CN raster
cn <- rast("CN.tif")  

# Clamp CN values to a realistic range (30–100)
cn_clean <- clamp(cn, lower = 30, upper = 100)
cn_clean[is.na(cn_clean)] <- 70  # Fill NA values with average CN

# Create a precipitation raster using gamma-distributed simulations
# (Assumes 'simulated_precipitation' is already loaded in the environment)
p <- cn_clean  # Use CN raster as a template for dimensions
set.seed(123)  # For reproducibility
values(p) <- sample(simulated_precipitation, size = ncell(p), replace = TRUE)

# Calculate potential maximum retention (S) in mm
S <- (25400 / cn_clean) - 254

# Compute direct runoff (Qmm) using the SCS-CN formula
Qmm <- ifel(p > 0.2 * S, (p - 0.2*S)^2 / (p + 0.8*S), 0)

# Summarize runoff statistics
global(Qmm, c("min", "max", "mean"), na.rm = TRUE)

# Visualize runoff depth
plot(Qmm, main = "Runoff Depth (mm) from Simulated Precipitation")

# Save runoff raster (GeoTIFF format)
writeRaster(
  Qmm, 
  filename = "Qmm.tif",
  overwrite = TRUE
)



# =============================================
# PART 1: LIBRARIES AND DATA LOADING
# =============================================

library(sf)
library(terra)
library(raster)  # For rasterize()

# Define paths
gdb_path <- "path_to_your_geodatabase.gdb"
ksat_path <- "path_to_Ksat.tif"
output_dir <- "path_to_output_folder/"

# Load spatial layers
Study <- st_read(dsn = gdb_path, layer = "Study")
Indy_roads <- st_read(dsn = gdb_path, layer = "Indy_roads")

# Reproject to UTM Zone 16N
utm_crs <- 32616
Study_proj <- st_transform(Study, utm_crs)
Indy_roads_proj <- st_transform(Indy_roads, utm_crs)

# =============================================
# PART 2: CALCULATE DRAINAGE CAPACITY δ_i
# =============================================

# Set raster resolution (meters)
resolution <- 9.43
template <- raster(extent(Study_proj), res = resolution, crs = st_crs(Study_proj)$wkt)

# Add road segment lengths in meters
Indy_roads_proj$length_m <- as.numeric(st_length(Indy_roads_proj))

# Rasterize total road length per cell
r_length <- rasterize(as(Indy_roads_proj, "Spatial"), template, field = "length_m", fun = sum, background = 0)

# Compute road density r_i (m of road / m²)
cell_area <- resolution^2
r_density <- r_length / cell_area

# Calculate mean road density (r)
mean_road_density <- cellStats(r_density, stat = "mean", na.rm = TRUE)

# Load saturated hydraulic conductivity (K_sat, mm/hr)
K_sat <- rast(ksat_path)

# Resample K_sat if not aligned
if (!compareGeom(r_density, K_sat, stopOnError = FALSE)) {
  K_sat <- resample(K_sat, r_density, method = "bilinear")
}

# Set reduction coefficient alpha (0 < α < 1)
alpha <- 0.7  # adjust as needed

# Calculate normalized road density
r_density_norm <- r_density / mean_road_density

# Compute drainage capacity δ_i
delta_i <- alpha * K_sat * r_density_norm

# Save drainage capacity δ_i raster
writeRaster(delta_i, paste0(output_dir, "Drainage_Capacity_delta_i.tif"), overwrite = TRUE)

# Optional plot
plot(delta_i, main = "Drainage Capacity δᵢ")

# =============================================
# PART 3: PONDING CALCULATION
# =============================================

# Load runoff raster Qmm (replace path accordingly)
Qmm <- rast("path_to_Qmm.tif")

# Crop Qmm to δ_i extent
Qmm_cropped <- crop(Qmm, ext(delta_i))

# Match resolution and extent
delta_i_resampled <- resample(delta_i, Qmm_cropped, method = "bilinear")

# Calculate ponding depth H_i = Qmm - δ_i
H_i <- Qmm_cropped - delta_i_resampled
H_i <- clamp(H_i, lower = 0)

# Save ponding output
writeRaster(H_i, paste0(output_dir, "Ponding_Depth_Fi.tif"), overwrite = TRUE)

# =============================================
# PART 4: VALIDATION OUTPUTS
# =============================================

cat("\nRunoff (Qmm) Statistics [mm]:\n")
print(global(Qmm, c("min", "max", "mean", "sd"), na.rm = TRUE))

cat("\nDrainage Capacity (δᵢ) Statistics [mm/hr]:\n")
print(global(delta_i, c("min", "max", "mean", "sd"), na.rm = TRUE))

cat("\nPonding Depth (Fᵢ) Statistics [mm]:\n")
print(global(H_i, c("min", "max", "mean", "sd"), na.rm = TRUE))

# Visual checks
plot(Qmm, main = "Runoff Qmm")
plot(delta_i, main = "Drainage Capacity δᵢ")
plot(H_i, main = "Ponding Depth Fᵢ")






################################################
### PF - Flood Probability Simulation Script
################################################

# Load necessary libraries
library(terra)
library(MASS)

# -------------------------------------------
# STEP 1: Use fitted gamma parameters from real precipitation
# (Assuming prcp_data was previously loaded and gamma_fit computed)
# Example from MARTINSVILLE or selected station:
prcp_data <- prcp_data[prcp_data > 0]  # Ensure positive-only values
gamma_fit <- fitdistr(prcp_data, "gamma")
alpha <- gamma_fit$estimate["shape"]
rate <- gamma_fit$estimate["rate"]
beta <- 1 / rate
mean_prcp <- mean(prcp_data)

# -------------------------------------------
# STEP 2: Load ponding depth raster
# -------------------------------------------
surf_ponding <- rast("surf_ponding.tif")

# Validate raster content
if (!hasValues(surf_ponding)) {
  stop("Error: 'surf_ponding.tif' has no values.")
}

# Quick checks
cat("Max ponding:", global(surf_ponding, max, na.rm = TRUE), "\n")
cat("Min ponding:", global(surf_ponding, min, na.rm = TRUE), "\n")

# -------------------------------------------
# STEP 3: Monte Carlo Simulation
# -------------------------------------------
thresholds <- c(0, 25, 50, 75, 100)  # Very-low to very-high (in mm)
labels <- c("very_low", "low", "moderate", "high", "very_high")

# Initialize stack to store flood counts for each class
n_sim <- 10000
set.seed(123)

PF_stack <- rast(H_i)
names(PF_stack) <- labels
values(PF_stack) <- 0

for (i in 1:n_simulations) {
  rand_prcp <- rgamma(1, shape = alpha, scale = beta)
  scale_factor <- rand_prcp / mean_prcp
  
  # Compute flood condition
  scaled_ponding <- surf_ponding * scale_factor
  flood_events <- scaled_ponding > threshold
  
  # Accumulate flood events
  flood_count <- flood_count + flood_events
  
  cat("Simulation", i, "– Max count so far:", global(flood_count, max, na.rm = TRUE), "\n")
}

# -------------------------------------------
# STEP 4: Compute PF = flood_count / n_simulations
# -------------------------------------------
PF <- flood_count / n_simulations

# Save output
writeRaster(PF, "PF.tif", format = "GTiff", overwrite = TRUE)

# -------------------------------------------
# STEP 5: Visualization
# -------------------------------------------
plot(PF,
     main = "Flood Probability (PF)",
     col = hcl.colors(100, "Blues", rev = TRUE),
     zlim = c(0, 1))

# Histogram of ponding
hist(values(surf_ponding), main = "Ponding Depth Distribution")

# Optional: Preview scale factors
scale_samples <- rgamma(100, shape = alpha, scale = beta) / mean_prcp
summary(scale_samples)

# Optional: Single test iteration
test_scale <- rgamma(1, shape = alpha, scale = beta) / mean_prcp
test_result <- surf_ponding * test_scale > threshold
plot(test_result, main = "Example Flood Event")




### Memory-Efficient Flood Probability Simulation
# ----------------------------------------------
# STEP 1: Use fitted gamma parameters from real precipitation
# (Assuming prcp_data was previously loaded and gamma_fit computed)
# Example from MARTINSVILLE or selected station:
prcp_data <- prcp_data[prcp_data > 0]  # Ensure positive-only values
gamma_fit <- fitdistr(prcp_data, "gamma")
alpha <- gamma_fit$estimate["shape"]
rate <- gamma_fit$estimate["rate"]
beta <- 1 / rate
mean_prcp <- mean(prcp_data)

# ----------------------------------------------
# STEP 2: Load ponding raster (H_i)
# ----------------------------------------------
surf_ponding <- raster("rasters/surf_ponding.tif")  # Replace with actual filename

if (!hasValues(surf_ponding)) {
  stop("Ponding raster has no values. Check file path or file integrity.")
}

print(paste("Max ponding:", maxValue(surf_ponding)))
print(paste("Min ponding:", minValue(surf_ponding)))

# ----------------------------------------------
# STEP 3: Monte Carlo simulation of flood probability (PF)
# ----------------------------------------------
thresholds_mm <- c(0, 25, 50, 75, 100)  # In mm: 0–2.5, 2.5–5, 5–7.5, 7.5–10, 10+
labels <- c("very_low", "low", "moderate", "high", "very_high")
n_classes <- length(labels)

output_files <- paste0("outputs/flood_count_", labels, ".tif")
flood_counts <- vector("list", n_classes)

for (k in 1:n_classes) {
  flood_counts[[k]] <- writeStart(
    raster(surf_ponding),
    filename = output_files[k],
    format = "GTiff",
    overwrite = TRUE
  )
}

n_simulations <- 10000
set.seed(123)
random_precips <- rgamma(n_simulations, shape = alpha, scale = beta)

bs <- blockSize(surf_ponding, minblocks = 20)
cat("Processing", bs$n, "blocks\n")

for (j in 1:bs$n) {
  row_start <- bs$row[j]
  n_rows <- bs$nrows[j]
  
  block_values <- getValues(surf_ponding, row = row_start, nrows = n_rows)
  block_counts <- matrix(0L, nrow = length(block_values), ncol = n_classes)
  
  for (i in 1:n_simulations) {
    scale_factor <- random_precips[i] / mean_prcp
    scaled_values <- block_values * scale_factor
    
    for (k in 1:n_classes) {
      lower <- thresholds_mm[k]
      upper <- ifelse(k < n_classes, thresholds_mm[k + 1], Inf)
      block_counts[, k] <- block_counts[, k] + (scaled_values > lower & scaled_values <= upper)
    }
  }
  
  for (k in 1:n_classes) {
    flood_counts[[k]] <- writeValues(flood_counts[[k]], block_counts[, k], row_start)
  }
  
  if (j %% 5 == 0) {
    cat("Completed block", j, "of", bs$n, "\n")
    gc()
  }
}

# ----------------------------------------------
# STEP 4: Finalize all outputs
# ----------------------------------------------
for (k in 1:n_classes) {
  flood_counts[[k]] <- writeStop(flood_counts[[k]])
  
  # Convert to probability and overwrite
  PF_raster <- flood_counts[[k]] / n_simulations
  writeRaster(PF_raster,
              filename = paste0("outputs/PF_", labels[k], ".tif"),
              format = "GTiff",
              overwrite = TRUE)
}

cat("All PF rasters saved.\n")



# --------------------------------------
# Final Flood Probability (PF) Calculation
# --------------------------------------

# Compute Flood Probability
PF <- flood_count / n_simulations

# Verify distribution
hist(PF, main = "Flood Probability Distribution", xlab = "Flood Probability")

# Save result
writeRaster(PF, "PF.tif", overwrite = TRUE)




