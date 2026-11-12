# Theta and clinical correlations

[`partial_correlations.m`](partial_correlations.m) analyzes morning and evening-minus-morning theta power against levodopa equivalent daily dose (LEDD) in EPD, ADV, and DYS, and against UDysRS Part IV in DYS. Its Figure 5 block uses the prefrontal electrodes F3, F4, and Fz.

## Inputs and setup

- MATLAB with the Statistics and Machine Learning Toolbox.
- `psdIndividual.mat` from EEG feature extraction and `clinical_data.xlsx` with subject ID, group, age, sex, LEDD, and UDysRS fields.
- `build_PSD_data_from_struct.m` from `../hdeeg-analysis/statistical-models/` on the MATLAB path.
- Set the placeholder paths near the start of the script to the local input files. The input datasets are not included in this repository.

## Analysis and outputs

The script extracts theta power over 4–8 Hz, computes morning power and the evening-minus-morning log-power difference, and aligns these measures with clinical records. It regresses age and sex from each EEG and clinical measure, then correlates the residuals. It reports Pearson correlations and Fisher-z confidence intervals. The Figure 5 block applies Benjamini–Hochberg FDR correction to the six LEDD tests; the two UDysRS tests are exploratory.

Outputs include `ThetaMorning_vs_LEDD_results.xlsx`, `Fig5_All_Correlations.xlsx`, and PNG correlation figures in the MATLAB working directory.

This script does not include the paper's supplementary correlations between theta and actigraphy-derived sleep measures.
