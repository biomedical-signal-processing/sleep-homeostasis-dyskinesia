# -------------------------------------------------------------
#  GGIR Pipeline for Actigraphy Data with Sleep Log Integration
# -------------------------------------------------------------
# This script performs a complete GGIR analysis on accelerometer data.
# It includes all five processing stages and integrates a sleeplog
# to guide sleep detection and calculate Time in Bed (TiB) and Sleep Latency.
# -------------------------------------------------------------
 
#  Load required libraries
library(GGIR)
library(readxl)
library(dplyr)
 
#  Define data and output paths
datadir <- "C:/Users/Galati/Desktop/prova/file"       # Path to raw accelerometer data
outputdir <- "C:/Users/Galati/Desktop/prova/out65"    # Path to store GGIR output
 
#  Path to the sleep log in CSV format
sleeplog_csv <- "C:/Users/Galati/Desktop/prova/Sleep_log_7d/check_clean.csv"
 
#  Read the sleep log
sleeplog_df <- read.csv(sleeplog_csv, stringsAsFactors = FALSE)
 
#  Inspect the sleep log structure
print(head(sleeplog_df))   # Display first few rows
str(sleeplog_df)           # Display column structure and types
 
#  Run GGIR with all stages and sleep log integration
GGIR(
  datadir = datadir,
  outputdir = outputdir,
 
  #  Enable all analysis stages (1 to 5)
  mode = c(1, 2, 3, 4, 5),  
  overwrite = TRUE,
 
  #  Generate summary reports for specific stages
  do.report = c(2, 4, 5),
 
  #  Enable data visualizations and Cosinor model for circadian rhythm analysis
  do.visual = TRUE,
  cosinor = TRUE,
 
  #  Sleep log parameters
  loglocation = sleeplog_csv,  # Path to sleep log file
  coln1 = 2,                   # Column number for sleep onset times
  colid = 1,                   # Column number for subject ID
 
  #  Specify the sleep window type to compute Time in Bed (TiB) and Sleep Latency
  sleepwindowType = "TimeInBed"
)
 
#  Final check for sleep log inclusion
print(" Sleeplog successfully processed and integrated into GGIR!")
print(head(sleeplog_df))
 
