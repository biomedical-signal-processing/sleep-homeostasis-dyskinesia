# 🧠 Sleep, Synaptic Homeostasis and Dyskinesia in Parkinson’s Disease

Welcome to this **umbrella repository** for code and resources related to the macro-project:

> **Disabling levodopa-induced dyskinesia in Parkinson's disease: the role of sleep**

The central aim is to understand how sleep interacts with levodopa-induced dyskinesia (LID) in Parkinson’s disease. 
During wake, synaptic strength tends to build up, while deep sleep and slow-wave activity support synaptic downscaling and network renormalization, processes that appear altered in patients with LID.
Investigating the association between the build-up process, the downscaling process, and LID would help 
establish a causal relationship between an alteration in this homeostatic process and the development of 
a disabling treatment-related complication.

This work is related to the following Swiss National Science Foundation project:  
<https://data.snf.ch/grants/grant/10001305>

---

## Repository Structure

This repository is organized into **subfolders**, each corresponding to a specific published study.

Each study folder contains:

- Analysis code (preprocessing, feature extraction, statistics, figures, etc..).  
- A **study-specific `README.md`** describing the analysis layout, inputs, and study results.

Use this top-level README as an entry point, then dive into the subfolder of the study you are interested in.

---

## Included Studies

### 1. Altered wakeful theta activity characterizes levodopa-induced dyskinesia in Parkinson’s disease

The present study aims to investigate **wake-theta activity** in resting states in advanced PD patients with and without dyskinesia. 
The central focus is to analyze and demonstrate that theta activity in advanced dyskinetic patients is already significantly elevated in the morning and remains stable throughout the day, unlike in controls and other PD groups, where theta activity is lower in the morning and progressively increases as wakefulness extends.

Fiorillo, L., Lombardi, G., La Porta, N., Arnaud, L., Veneruso, M., Castelnovo, A., Bertaina, I., Staedler, C., Kaelin-Lang, A., and Galati, S.

*Altered wakeful theta activity characterizes levodopa-induced dyskinesia in Parkinson’s disease.* npj Parkinson’s Disease (2026).  
<https://doi.org/10.1038/s41531-026-01320-z>


👉 Start here for this project:  
[`npjParkinsonDisease2026/README.md`](npjParkinsonDisease2026/README.md)

---

### 2. Dynamics of slow wave sleep in advanced Parkinson's disease: An entropy-based study

This study focuses on **slow-wave sleep dynamics** in advanced Parkinson’s disease, using entropy-based metrics to capture nocturnal alterations that may relate to dyskinesia beyond classical slow-wave power.
The goal is to determine whether nonlinear signatures could capture complementary or more robust alterations in sleep-related brain activity associated with LID. 


Lombardi, G., Fiorillo, L., Castelnovo, A., Bertaina, I., Staedler, C., Kaelin-Lang, A., and Galati, S.

*Dynamics of slow wave sleep in advanced Parkinson's disease: An entropy-based study.* Neurobiology of Disease (2026).  
<https://doi.org/10.1016/j.nbd.2026.107395>


👉 Start here for this project:  
[`NeurobiologyOfDisease2026/README.md`](NeurobiologyOfDisease2026/README.md)

---

## How to Use This Repository

- Choose a study above and open its README for the script layout and study context.
- Configure each script's input and output paths and obtain the required study data and external dependencies before running it.
- The two studies use different group labels and analysis pipelines; follow the conventions documented within each study folder.
- The study datasets are not included here. Both papers state that supporting data may be requested from the corresponding author; see the Data Access section in each study README.
