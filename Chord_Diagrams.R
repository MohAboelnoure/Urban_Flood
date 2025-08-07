# -----------------------------------------------------------
# Chord Diagram of Flood Exposure Pathways (FEP) - Historical to RCP4.5
# -----------------------------------------------------------

# Load required libraries
library(circlize)
library(dplyr)
library(tidyr)

# -----------------------------
# Load and prepare the data
# -----------------------------

# Read the input CSV file
# Replace with your own relative or absolute file path
df_raw <- read.csv("data/FEP_Chord_Hist_rcp45.csv")

# Assign raw data to working dataframe
df <- df_raw

# Define the class order
levels_order <- c("Very low", "Low", "Medium", "High", "Very high")
df$CLASS_FROM <- factor(df$CLASS_FROM, levels = levels_order)
df$CLASS_TO   <- factor(df$CLASS_TO, levels = levels_order)

# Prepare matrix for chord diagram
chord_matrix <- df %>%
  filter(CLASS_FROM != "Same", CLASS_TO != "Same") %>%
  select(CLASS_FROM, CLASS_TO, AREA) %>%
  pivot_wider(names_from = CLASS_TO, values_from = AREA, values_fill = 0)

chord_matrix <- as.data.frame(chord_matrix)
rownames(chord_matrix) <- chord_matrix$CLASS_FROM
chord_matrix <- chord_matrix[, -1]
chord_matrix <- as.matrix(chord_matrix)

# Clear any existing plots
circos.clear()

# Define color mapping for each category
grid_colors <- c(
  "Very low"  = "#006400",
  "Low"       = "#9ACD32",
  "Medium"    = "#FFFF00",
  "High"      = "#FFA500",
  "Very high" = "#FF0000"
)

# Calculate total area per source class for axis scaling
area_totals <- df %>%
  filter(CLASS_FROM != "Same", CLASS_TO != "Same") %>%
  group_by(CLASS_FROM) %>%
  summarise(total_area = sum(AREA, na.rm = TRUE)) %>%
  arrange(factor(CLASS_FROM, levels = levels_order))

# -----------------------------
# Plot the chord diagram
# -----------------------------

chordDiagram(
  chord_matrix,
  grid.col = grid_colors,
  transparency = 0.2,
  directional = 1,
  direction.type = c("arrows", "diffHeight"),
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.05)
)

# Customize the outer ring
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  
  sector_name <- CELL_META$sector.index
  sector_xlim <- CELL_META$xlim
  
  # Match sector name to area total
  total_area <- area_totals$total_area[area_totals$CLASS_FROM == sector_name]
  
  # Add axis and ticks
  circos.axis(
    h = "top",
    major.at = seq(0, total_area, length.out = 4),
    labels = round(seq(0, total_area, length.out = 4), 0),
    labels.cex = 0.7,
    major.tick.length = 0.5,
    sector.index = sector_name,
    track.index = 2
  )
  
  # Add sector label
  circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(5),
              sector_name, facing = "bending.outside",
              niceFacing = TRUE, adj = c(0.5, 0),
              cex = 1.5, font = 2)
  
}, bg.border = NA)

# -----------------------------
# Save the plot as PNG
# -----------------------------

# Store the plot
saved_plot <- recordPlot()

# Create output directory if it doesn't exist
if (!dir.exists("outputs")) dir.create("outputs")

# Save the plot to PNG
png("outputs/FEP_Chord_Hist_rcp45.png", width = 10, height = 10, units = "in", res = 600)
replayPlot(saved_plot)
dev.off()
