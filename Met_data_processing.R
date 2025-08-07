# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(MASS)

# Optional: Set working directory to current project folder
# setwd("your/working/directory")  # Uncomment and modify if needed

# Define relative file paths for each station
data_dir <- "Met_data/data"  # Make sure this folder is in your project directory

station_files <- list(
  MARTINSVILLE = "data_39.40625_-86.46875",
  FRANKLIN = "data_39.46875_-86.03125",
  SHELBYVILLE = "data_39.53125_-85.78125",
  AIRPORT = "data_39.71875_-86.28125",
  DANVILLE = "data_39.71875_-86.53125",
  GREENFIELD = "data_39.78125_-85.78125",
  CASTLETON = "data_39.90625_-86.03125",
  JAMESTOWN = "data_39.90625_-86.59375",
  ANDERSON = "data_40.09375_-85.65625",
  TIPTON = "data_40.21875_-86.09375"
)

# Load datasets
stations <- list()
for (name in names(station_files)) {
  filepath <- file.path(data_dir, station_files[[name]])
  if (file.exists(filepath)) {
    stations[[name]] <- read.table(filepath, quote = "\"", comment.char = "", header = TRUE)
  } else {
    warning(paste("File not found:", filepath))
  }
}

# Summary statistics
for (name in names(stations)) {
  cat("\n### Summary for", name, "###\n")
  print(summary(stations[[name]]))
}

# Annual mean precipitation
annual_means <- list()
for (name in names(stations)) {
  data <- stations[[name]]
  if (!("time" %in% names(data)) || !("prcp" %in% names(data))) next
  data$time <- as.Date(data$time)
  data$year <- format(data$time, "%Y")
  annual_totals <- aggregate(prcp ~ year, data, sum)
  annual_means[[name]] <- mean(annual_totals$prcp)
}
cat("\n### Annual Means ###\n")
print(annual_means)

# Fit Gamma distribution to one example station (e.g., MARTINSVILLE)
if ("MARTINSVILLE" %in% names(stations)) {
  prcp_data <- stations$MARTINSVILLE$prcp
  prcp_data <- prcp_data[prcp_data > 0]
  gamma_fit <- fitdistr(prcp_data, "gamma")
  alpha <- gamma_fit$estimate["shape"]
  rate <- gamma_fit$estimate["rate"]
  beta <- 1 / rate
  
  cat("Gamma fit for MARTINSVILLE:\n")
  cat("Shape (alpha):", alpha, "\n")
  cat("Scale (beta):", beta, "\n")
  
  hist(prcp_data, breaks = 30, prob = TRUE,
       main = "Fitted Gamma - MARTINSVILLE", col = "gray", border = "black")
  curve(dgamma(x, shape = alpha, rate = rate), col = "blue", lwd = 2, add = TRUE)
}

# Fit and simulate for all stations
for (name in names(stations)) {
  data <- stations[[name]]
  prcp_data <- data$prcp[data$prcp > 0]
  if (length(prcp_data) == 0) next
  gamma_fit <- fitdistr(prcp_data, "gamma")
  alpha <- gamma_fit$estimate["shape"]
  rate <- gamma_fit$estimate["rate"]
  beta <- 1 / rate
  simulated <- rgamma(10000, shape = alpha, scale = beta)
  
  cat("\n", name, "Simulated Mean:", mean(simulated), "\n")
  
  # Histogram
  hist(simulated, breaks = 30, prob = TRUE,
       main = paste("Simulated Precipitation -", name),
       col = "gray", xlab = "Precipitation")
  curve(dgamma(x, shape = alpha, scale = beta), col = "blue", lwd = 2, add = TRUE)
}

# Panel plot (all stations)
par(mfrow = c(2, 5), mar = c(4, 4, 2, 1))
for (name in names(stations)) {
  data <- stations[[name]]
  prcp_data <- data$prcp[data$prcp > 0]
  if (length(prcp_data) == 0) next
  gamma_fit <- fitdistr(prcp_data, "gamma")
  alpha <- gamma_fit$estimate["shape"]
  rate <- gamma_fit$estimate["rate"]
  beta <- 1 / rate
  simulated <- rgamma(10000, shape = alpha, scale = beta)
  
  hist(simulated, breaks = 30, prob = TRUE, main = name,
       xlab = "Precipitation", col = "gray", border = "black")
  curve(dgamma(x, shape = alpha, scale = beta), col = "blue", lwd = 2, add = TRUE)
}
par(mfrow = c(1, 1))

# CDF plots (all stations)
cdf_data_list <- list()
for (name in names(stations)) {
  prcp_data <- stations[[name]]$prcp
  prcp_data <- prcp_data[prcp_data > 0]
  if (length(prcp_data) == 0) next
  cdf <- ecdf(prcp_data)
  cdf_data <- data.frame(
    Precipitation = sort(prcp_data),
    CDF = cdf(sort(prcp_data)),
    Station = name
  )
  cdf_data_list[[name]] <- cdf_data
}
all_cdf_data <- do.call(rbind, cdf_data_list)

ggplot(all_cdf_data, aes(x = Precipitation, y = CDF, color = Station)) +
  geom_step(size = 0.7) +
  labs(title = "CDFs of Daily Precipitation", x = "Precipitation", y = "Cumulative Probability") +
  theme_minimal() +
  theme(legend.position = "bottom")


