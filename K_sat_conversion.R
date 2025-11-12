
library(terra)

# --- Define Input and Output Paths ---
# (Replace with your actual file paths)
ksat_input  <- "data/SSURGO_Ksat_um_per_s.tif"       
ksat_output <- "outputs/SSURGO_Ksat_mm_per_day.tif"  

# --- Load Ksat Raster ---
Ksat_um_s <- rast(ksat_input)
cat("Original Ksat units: micrometers per second (µm/s)\n")
cat("Raster summary before conversion:\n")
print(global(Ksat_um_s, c("min", "max", "mean"), na.rm = TRUE))

# --- Convert Units ---
# 1 µm/s = 0.001 mm/s
# 1 day = 86400 seconds
# Therefore, 1 µm/s = 86.4 mm/day
Ksat_mm_day <- Ksat_um_s * 86.4

#Define Δt (in days)
# Example: Δt = 1 for daily rainfall, 0.25 for 6-hour storm, etc.
delta_t <- 1  # <-- Adjust this as needed

# Compute event-based drainage depth (mm)
Ksat_mm_event <- Ksat_mm_day * delta_t

#Save Output ---
writeRaster(Ksat_mm_event, ksat_output, overwrite = TRUE)
cat("\nConversion complete.\n")
cat("Saved Ksat (drainage depth, mm over event duration Δt =", delta_t, "days):\n", ksat_output, "\n")


