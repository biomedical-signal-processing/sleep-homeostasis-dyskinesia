## Overview

This folder contains a R script processing actigraphy data from GeneActiv inertial sensors using the [GGIR package](https://cran.r-project.org/web/packages/GGIR/vignettes/GGIR.html). 

The pipeline is the following:
1. Processes raw accelerometer data
2. Integrates sleep logs for improved sleep detection
3. Calculates key estimated(e) sleep-wake metrics:
   - Sleep Efficiency (eSE)
   - Total Sleep Time (eTST)
   - Wake After Sleep Onset (eWASO)
   - Time in Bed (eTIB)
   - Sleep Latency (eSL)

## Requirements

### Core Software

- R (≥ 4.0.0)
- Required packages:
```
install.packages(c(“GGIR”, “readxl”, “dplyr”))
```

### Input Data

| Data Type        | Format             | Description                    |
|------------------|--------------------|--------------------------------|
| Actigraphy Data  | `.bin` (GeneActiv) | Raw accelerometer data         |
| Sleep Log        | `.csv`             | Participant sleep diaries      |

## Configuration

### Key Parameters

```
# File paths
datadir <- “C:/Users/Galati/Desktop/prova/file”       # Raw accelerometer data
outputdir <- “C:/Users/Galati/Desktop/prova/out65”    # Output directory
sleeplog_csv <- “C:/Users/Galati/Desktop/prova/Sleep_log_7d/check_clean.csv”  # Sleep log

# GGIR parameters
mode = c(1, 2, 3, 4, 5)   # Enable all processing stages
sleepwindowType = “TimeInBed”  # Compute TiB and Sleep Latency
colid = 1                   # Subject ID column in sleep log
coln1 = 2                   # Sleep onset time column
```

## Workflow

### 1. Data Preparation

- **Actigraphy Data**: Place raw `.bin` files in `datadir`
- **Sleep Log**: CSV file with columns:

  | Subject_ID | Sleep_Onset | ... |
  |------------|-------------|-----|
  | P001       | 22:30       | ... |
  | P002       | 23:15       | ... |

### 2. GGIR Processing Stages

| Stage | Function                | Output                        |
|-------|-------------------------|-------------------------------|
| 1     | Raw data processing     | `.csv` time series            |
| 2     | Activity monitoring     | Daily activity summaries      |
| 3     | Sleep detection         | Sleep period identification   |
| 4     | Sleep analysis          | eSE, eTST, eWASO              |
| 5     | 24-hour rhythms         | Circadian patterns            |

### 3. Sleep Log Integration

GGIR uses sleep logs to:
- Improve sleep window detection
- Calculate Time in Bed
- Calculate Sleep Latency
- Validate automated sleep detection

## Output Files

GGIR generates in `outputdir`:
```
/outputdir
├── results/
│   ├── part2_summary.csv       # Activity summaries
│   ├── part3_summary.csv       # Sleep detection
│   ├── part4_summary.csv       # Sleep metrics
│   └── part5_summary.csv       # Circadian rhythms
└── results/QC/
    └── visualreport.pdf        # Data quality report
```

### Key Sleep Metrics

| Metric         | Description                      | Calculation                   |
|----------------|----------------------------------|-------------------------------|
| eSE            | Estimated Sleep Efficiency       | (TST / TiB) × 100             |
| eTST           | Estimated Total Sleep Time       | Total sleep minutes           |
| eWASO          | Estimated Wake After Sleep Onset | Wake minutes after sleep onset|
| eTiB           | Estimated Time in Bed            | Sleep log: Bed exit - entry   |
| eSleep Latency | Estimated Time to fall asleep    | Sleep onset - Bed entry       |

## Customization

### Sleep Detection Parameters

```
# Sleep algorithm tuning
def.noc.sleep = c(1, 1, 0),  # Sleep detection thresholds
HASPT.algo = “HorAngle”,     # Sleep period algorithm

# Advanced parameters
anglethreshold = 5,          # Wrist angle threshold
timethreshold = 5,           # Time threshold (minutes)
```

### Activity Analysis

```
# Activity intensity thresholds
mvpathreshold = c(100, 400),   # MVPA thresholds (mg)
bout.metric = 4,               # Bout detection metric
```

## Contact

For technical support:
[Luigi Fiorillo](mailto:luigi.fiorillo@eoc.ch),
[Giovanni Lombardi](mailto:giovanni.lombardi@eoc.ch)
