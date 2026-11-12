## Overview

This folder contains a MATLAB pipeline to process 64-channel ANT Neuro EEG data through:
1. `preprocessing.m`: Converts raw `.cnt` files to cleaned `.set` files
2. `feature_extraction.m`: Calculates theta-band (4-8 Hz) power spectral density (PSD) from `.set` files


## Requirements

### Core Software

- MATLAB (R2020a or newer)
- EEGLAB 2024.2 ([Installation Guide](https://eeglab.org/tutorials/Installation))

### Toolboxes

- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox

### Dependencies

| Tool | Purpose | Source |
|------|---------|--------|
| FieldTrip | Head modeling | [fieldtriptoolbox.org](https://www.fieldtriptoolbox.org/) |
| ICLabel | Component classification | [GitHub](https://github.com/sccn/ICLabel) |
| ANTeepimport | ANT Neuro data import | Included |
| csc-eeg-tools | Artifact removal | Included |

## Data Organization

Organize input data as follows:

```
ThetaProject/
└── EEGDATA/
    └── PXXX/             # Subject ID (e.g., P065)
        └── BSL/          # Session
            └── original/
                ├── *s1_oc.cnt   # Condition file
                └── …            # Additional conditions
```

## Preprocessing Workflow (`preprocessing.m`)

### Key Steps

1. **Data Import**  
   - Reads ANT Neuro `.cnt` files using FieldTrip
2. **Filtering**  
   - 0.5 Hz high-pass → 50 Hz notch → 100 Hz low-pass
3. **Artifact Handling**  
   - Manual bad channel identification
   - Artifact Subspace Reconstruction (ASR)
4. **ICA Processing**  
   - RUNICA decomposition
   - ICLabel classification (muscle/eye/heart components)
   - Interactive component rejection
5. **Output**  
   Saves `.set` files in group folders (CTL/EPD/ADV/DYS)

### Configuration

```
ID = ‘P065’;             % Subject ID
condition = ‘*s1_oc’;    % File pattern
group = ‘EPD’;           % Group classification
dataPath = ‘C:...\data_processed_BSL’; % Output directory
```

### Execution
Run script section-by-section (includes pauses for manual verification)

## Feature Extraction (`feature_extraction.m`)

### Key Parameters

```
window = 512;       % Welch window
noverlap = 256;     % 50% overlap
nfft = 512;         % FFT points
notchBand = 49.5, 50.5; % Notch filter
```

### Key Steps

1. **Channel Processing**  
   - Automatic removal of EOG channel (32)
2. **PSD Calculation**  
   - Welch's method per channel
   - Frequency-domain notch filtering
3. **Output**  
   Saves `psdIndividual.mat` containing:
   - Group/session-wise theta PSD
   - Frequency vector (4-8 Hz band)
   - 
### Configuration

```
dataPath = ‘C:...\data_processed_BSL\BRAIN_…’; % Input .set location
outputPath = ‘C:...\dati_PSD’;                   % Output directory
groups = {‘CTL’, ‘EPD’, ‘ADV’, ‘DYS’};            % Groups to process
```

## Output Structure

```
data_processed_BSL/
├── CTL/                             # Group folders
│   ├── 065_s1_oc_preprocessedBSL.set
│   └── … 
├── EPD/
└── dati_PSD/                        # Feature output
    ├── CTL/
    │   ├── morning/                 # Session data
    │   └── evening/
    ├── EPD/
    └── psdIndividual.mat           # Consolidated features
```

## Contact

For technical support:
[Luigi Fiorillo](mailto:luigi.fiorillo@eoc.ch),
[Giovanni Lombardi](mailto:giovanni.lombardi@eoc.ch)
