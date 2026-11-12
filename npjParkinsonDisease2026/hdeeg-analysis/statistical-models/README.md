## Overview

This folder contains three MATLAB scripts for analyzing EEG theta-band (4-8 Hz) power:
1. **Qualitative Topoplots (`qualitative_topoplots.m`)**  
   Visualizes spatial distributions of theta power across groups and sessions.
2. **Morning Analysis (`morning_theta.m`)**  
   Compares group differences (e.g., CTL vs. EPD/ADV/DYS) in morning sessions using linear models with clinical covariates.
2. **Diurnal Analysis (`diurnal_theta.m`)**  
   Examines `Group × Session` interactions using linear mixed-effects models with permutation-based TFCE correction.

## Requirements

### Core Software

- MATLAB (R2020a+)
- Statistics and Machine Learning Toolbox
- Parallel Computing Toolbox (recommended)

### Dependencies

| Tool | Purpose | Source |
|------|---------|--------|
| `ANT_64.ced` | Electrode locations | External input |
| `build_PSD_data_from_struct.m` | Data organizer | Included |
| `lmeEEG_TFCE_mod.m` | TFCE implementation | External dependency |

## Data Preparation

1. **Input**: `psdIndividual.mat` from feature extraction pipeline
2. **Clinical Data**: `clinical_data.xlsx` with columns: ID, Age, Sex, Disease_duration, LEDD

## Configuration

### Key Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `referenceGroup` | Baseline group | `'CTL'` |
| `referenceSession` | Baseline session | `'morning'` |
| `theta_band` | Frequency range | 4-8 Hz |
| `exclude` | Channels to exclude | `[13 19 32]` |
| `nPerms` | TFCE permutations | `2000` |

### Customization

- Change baseline group/session by modifying:
```
referenceGroup = ‘DYS’;  % in script diurnal_theta.m
referenceSession = ‘evening’;
```

- Adjust theta band:
```
theta_idx = find(freqs >= 4 & freqs <= 7);  % Modify in both scripts
```

## Qualitative Topoplots (`qualitative_topoplots.m`)

### Workflow

1. **Data Loading**:  
   - Loads group-level PSD data (morning and evening sessions)
2. **Theta Power Calculation**:  
   - Averages PSD within theta band (default: 4-8 Hz)
   - Computes log10 power and session differences
3. **Topoplot Generation**:  
   - Creates 4×3 panel figure (groups × sessions/differences)
   - Morning, evening, and evening-morning difference maps

### Outputs

- `report/theta_X-YHz_clean.png` (and .fig/.tiff) containing:
  - Rows: CTL, EPD, ADV, DYS groups
  - Columns: Morning, Evening, Evening-Morning difference

## Morning Analysis (`morning_theta.m`)

### Workflow

1. **Data Loading**:  
 - Loads PSD data and clinical covariates
2. **Model Fitting**:  
 Robust linear model per channel:  
 `Theta Power ~ Group + Age + Sex + Duration + LEDD`
3. **Effect Size**:  
 Computes Cohen's d for group differences
4. **TFCE Correction**:  
 Non-parametric permutation testing with TFCE

### Outputs

- `result_morning.mat` containing:
- t-values, p-values, Cohen's d per channel
- TFCE-corrected statistical maps
- Group comparison topoplots

## Diurnal Analysis (`diurnal_theta.m`)

### Workflow

1. **Linear Mixed-Effects Model**:  
 Robust linear model per channel:  
 `Theta Power ~ Group*Session + Age + Sex + LEDD + Duration + (1|Id)`


2. **Effects Estimated**:
- Group × Session interaction
- Main effect of Group
- Main effect of Session

3. **TFCE Correction**:  
2000 permutations with channel-wise TFCE

### Outputs
- `results_interaction.mat` containing coefficients for all model terms, t-values and effect sizes
- `tfce_results_interaction.mat` with TFCE outputs, i.e., 2×3 topoplot matrix (TFCE and -log10(p))

## Statistical Parameters

### Model Specifications

| Analysis Type | Model Family | Random Effects | Covariates |
|---------------|--------------|------------|---------|
| Morning | Robust Linear | None       | Age, Sex, LEDD, Duration |
| Diurnal | Linear Mixed-Effects | Subject ID | Age, Sex, LEDD, Duration |

### TFCE Parameters

```
E_H = 0.66 2;  % Height/Extent parameters
nPerms = 2000;    % Permutation count
```

## Contact

For technical support:
[Luigi Fiorillo](mailto:luigi.fiorillo@eoc.ch),
[Giovanni Lombardi](mailto:giovanni.lombardi@eoc.ch)
