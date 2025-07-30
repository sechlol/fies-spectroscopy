# FIES-Spectroscopy

Spectral analysis of FIES data for determination of stars' rotational velocity and chromospheric activity.

## Overview

This project focuses on the analysis of spectroscopic data from the FIber-fed Echelle Spectrograph (FIES) at the Nordic Optical Telescope (NOT) located in La Palma, Spain. The primary goal is to study young, magnetically active solar-type stars through high-resolution spectroscopy.

### Project Objectives

1. **Correction of Spectra for Radial Velocity (vrad)**: Correct spectral distortions caused by stellar motion along the line of sight.
2. **Continuum Normalization**: Normalize spectra against a well-fitted continuum curve to isolate regions of interest.
3. **Calculating log R'HK and Comparing with Rotation**: Analyze the relationship between chromospheric activity and stellar rotation.
4. **Comparing Line Profiles of Chromospheric Lines**: Study variations in line profiles across a sample of young solar-type stars.
5. **Estimating v sin(i)**: Determine the projected rotational velocity of stars.

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

3. Download the spectral data:
   ```
   cd spectroscopy
   python download_fies_data.py
   ```

### Running the Pipeline

To process all the data and generate results:

```
cd spectroscopy
python pipeline.py
```

This will:
- Download FIES spectra (if not already present)
- Average multiple spectra for each target
- Calculate and correct for radial velocity
- Normalize the continuum
- Calculate chromospheric excess emission
- Estimate v sin(i) values
- Generate plots and save results

## Data

The project analyzes spectroscopic data from several young solar-type stars:
- V889Her
- V383Lac
- V453And
- BCet
- HNPeg
- EXCet
- V774Tau
- V834Tau

The raw data consists of FITS images that capture the spectra of light as recorded by the sensor, which are processed using the FIESTool pipeline.

## Analysis Components

The project includes several Jupyter notebooks that showcase each analysis step:

1. `1_spectrum_visualization.ipynb`: Visualization of the raw and processed spectra
2. `2_vrad_estimation.ipynb`: Estimation of radial velocity
3. `3_spectrum_normalization.ipynb`: Normalization of the continuum
4. `4_line_emission.ipynb`: Analysis of emission lines
5. `5_line_comparison.ipynb`: Comparison of chromospheric lines
6. `6_vsini_estimation.ipynb`: Estimation of v sin(i)

## Results

The analysis produces several outputs:
- Radial velocity corrections for each star
- Normalized spectra for analysis of specific spectral features
- Calculation of R'HK and Rossby number (Ro) for each star
- Estimation of v sin(i) values
- Comparative plots of chromospheric line features

All results are saved in the `data/out/` directory.
