# FIES Spectroscopy: Analysis of Stellar Activity and Rotation Parameters
[![DOI](https://zenodo.org/badge/1029184911.svg)](https://doi.org/10.5281/zenodo.16622703)

Spectral analysis of FIES data for determination of stars' rotational velocity and chromospheric activity.

## Overview

This project focuses on the analysis of spectroscopic data from the FIber-fed Echelle Spectrograph (FIES) at the Nordic Optical Telescope (NOT) located in La Palma, Spain. The primary goal is to study young, magnetically active solar-type stars through high-resolution spectroscopy.

Refer to the [report](analysis_report.pdf) for a detailed analysis and interpretation of the results.

## Quickstart

### Installation

1. Clone the repository:
   ```
   git clone https://github.com/sechlol/fies-spectroscopy.git
   cd fies-spectroscopy
   ```

2. Install the required dependencies:
   ```
   pip install -r requirements.txt
   ```

3. Download the spectral data and run the analysis pipeline:
   ```
   python -m spectroscopy.pipeline
   ```

The `-m` flag tells Python to run the module as a script, ensuring proper package imports regardless of the current working directory.


### Project Objectives

1. **Correction of Spectra for Radial Velocity (vrad)**: Correct spectral distortions caused by stellar motion along the line of sight.
2. **Continuum Normalization**: Normalize spectra against a well-fitted continuum curve to isolate regions of interest.
3. **Calculating log R'HK and Comparing with Rotation**: Analyze the relationship between chromospheric activity and stellar rotation.
4. **Comparing Line Profiles of Chromospheric Lines**: Study variations in line profiles across a sample of young solar-type stars.
5. **Estimating v sin(i)**: Determine the projected rotational velocity of stars.

## Project Structure

```
fies-spectroscopy/
├── analysis/           # Analysis notebooks
├── data/               # Data directory
│   ├── fits/           # FITS files
│   └── out/            # Results
├── spectroscopy/       # Python modules
└── paths.py            # Path management
```


## Data Analysis


The project analyzes spectroscopic data from several young solar-type stars:
- V889Her
- V383Lac
- V453And
- BCet
- HNPeg
- EXCet
- V774Tau
- V834Tau

The raw data consists of FITS images processed using the FIESTool pipeline.

### Jupyter Notebooks
The project includes several Jupyter notebooks in the `analysis` directory:

1. `1_spectrum_visualization.ipynb`: Visualizes raw spectra and demonstrates plotting capabilities for examining spectral features across different targets.

2. `2_vrad_estimation.ipynb`: Estimates and corrects for radial velocity (${\text{vrad}}$) to account for Doppler shifting due to stellar motion and Earth's rotation.

3. `3_spectrum_normalization.ipynb`: Normalizes the continuum in spectral regions of interest, particularly around ${\text{CaII}}$ emission bands, using polynomial fitting.

4. `4_line_emission.ipynb`: Calculates chromospheric activity indices ($\log R'_{HK}$) and compares them with rotation parameters (Rossby number) to analyze magnetic activity.

5. `5_line_comparison.ipynb`: Compares spectral line profiles across different stars, focusing on chromospheric emission lines to identify patterns in magnetic activity.

6. `6_vsini_estimation.ipynb`: Estimates projected rotational velocity ($v \sin i$) by analyzing the broadening of absorption lines and comparing with literature values.

## Results

The analysis produces several outputs:
- Radial velocity corrections for each star
- Normalized spectra for analysis of specific spectral features
- Calculation of $R'_{HK}$ and Rossby number ($R_o$) for each star
- Estimation of $v \sin i$ values
- Comparative plots of chromospheric line features

All results are saved under the `data/out/` directory.

## How to Cite

If you use this code or data in your research, please cite it as:

```bibtex
@software{sechlol_fies_spectroscopy_2025,
  author       = {Cardin, Christian},
  title        = {{FIES Spectroscopy: Analysis of Stellar Activity and Rotation Parameters}},
  month        = jul,
  year         = 2025,
  publisher    = {Zenodo},
  version      = {v1.0.0},
  doi          = {10.5281/zenodo.16622703},
  url          = {https://doi.org/10.5281/zenodo.16622703}
}
```

This project is archived on Zenodo: https://doi.org/10.5281/zenodo.16622703