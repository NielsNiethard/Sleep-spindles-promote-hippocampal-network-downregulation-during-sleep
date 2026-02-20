# Sleep spindles promote hippocampal network downregulation during sleep

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18714655.svg)](https://doi.org/10.5281/zenodo.18714655) 

This repository contains the MATLAB code required to reproduce the analysis and figures for the publication: **"Sleep spindles promote hippocampal network downregulation during sleep"** by Niels Niethard et al.

## Overview
The code provided here processes calcium imaging data across different sleep states (NREM, REM, and Wake) to analyze the frequency and amplitude of calcium transients. Specifically, it includes algorithms to:
* Detect and define active vs. inactive cells during specific sleep events (Slow Oscillations and Spindles).
* Extract sleep state triplets (e.g., SWS-REM-SWS, SWS-Wake-SWS) and doublets to analyze cross-state network dynamics.
* Perform statistical evaluations (Linear Mixed-Effects Models, ANOVAs, and non-parametric tests) on state-dependent calcium activity.

## System Requirements
* **MATLAB**: The scripts were developed and tested using MATLAB R2024b.
* **Dependencies**: 
  * [Mass Univariate ERP Toolbox](https://github.com/dmgroppe/Mass_Univariate_ERP_Toolbox) (must be downloaded and added to your MATLAB path).
  * Standard MATLAB toolboxes (Statistics and Machine Learning Toolbox, Curve Fitting Toolbox).
