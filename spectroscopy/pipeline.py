from pathlib import Path

import numpy as np
import pandas as pd

from download_fies_data import download_spectra
from spectroscopy.fies_utils import (
    estimate_v_rad,
    correct_spectrum_for_vrad,
    normalize_continuum_slice,
    slice_spectrum_data,
    estimate_vsini,
)
from spectroscopy.plotting_utils import save_rhk_ro_plot, save_chromospheric_lines_plot, save_vsini_comparison_plot
from spectroscopy.spectral_features import (
    get_spectral_features_vrad_correction,
    calculate_log_rhk,
    CONTINUUM_SLICE,
    calculate_log_ro,
    CaII_EMISSION_LINES,
)
from spectroscopy.spectrum_loader import SpectrumLoader

# Input data
_DATA_FOLDER = Path.cwd() / "data"
_FITS_DATA_FOLDER = _DATA_FOLDER / "fits"
_OBSERVATIONS_FILE = _DATA_FOLDER / "observations.csv"
_STAR_DATA = _DATA_FOLDER / "star_data.csv"

# Output data
_OUT_DATA_FOLDER = _DATA_FOLDER / "out"
_AVERAGED_FILE = _OUT_DATA_FOLDER / "1_averaged_spectra.npz"
_VRAD_CORRECTED_FILE = _OUT_DATA_FOLDER / "2_vrad_corrected_spectra.npz"
_CONTINUUM_NORMALIZED_FILE = _OUT_DATA_FOLDER / "3_continuum_corrected_spectra.npz"
_COMPARISON_PLOT = _OUT_DATA_FOLDER / "rhk_ro_comparison.png"
_VSINI_PLOT = _OUT_DATA_FOLDER / "vsini_comparison.png"
_DATAFRAME_FILE = _OUT_DATA_FOLDER / "results.csv"
_DATAFRAME_JSON_FILE = _OUT_DATA_FOLDER / "results.json"

# Fitting parameters
_fitting_default = {"degree": 4, "iterations": 8}
_fitting_settings = {"V834Tau": {"degree": 3, "iterations": 3}, "V383Lac": {"degree": 5, "iterations": 8}}


def average_multiple_spectra() -> dict[str, np.ndarray]:
    """
    For Each target we have multiple observations, with same exposure time.
    The final spectra is an average of all the observations.

    :return: Dictionary of {target: spectrum}, where the spectrum is a Nx2 numpy array.
    First column is wavelength in Å, second column is the intensity at that wavelength.
    """

    # These spectra were acquired on a different day (25/10/2023),
    # and the reduced spectra doesn't align with the one taken the next day (26/10/2023)
    excluded_spectra = ["FIGj250094", "FIGj250095"]

    loader = SpectrumLoader(fits_path=_FITS_DATA_FOLDER, observations_path=_OBSERVATIONS_FILE)

    spectra = {}
    for target in loader.target_list:
        w = []
        y = []

        for spectrum in loader.get_spectra(target):
            if spectrum.id not in excluded_spectra:
                w.append(spectrum.data[:, 0])
                y.append(spectrum.data[:, 1])

        w = np.array(w)

        # Assert that all spectra share the same wavelength values.
        assert np.all(w == w[0, :]), f"Spectra misalignment for target {target}"

        # Take the average of the spectra
        y = np.array(y).mean(axis=0)

        # Arrange the data in a two-column array
        spectra[target] = np.vstack([w[0, :], y]).T

    return spectra


def vrad_correction(spectra: dict[str, np.ndarray]) -> tuple[dict[str, float], dict[str, np.ndarray]]:
    """
    For each target, finds a suitable v_rad value to correct for doppler shifting due to the
    combined star's radial velocity and Earth's rotation.

    :return: A dictionary of corrected spectra, and the associated v_rad values
    """
    v_rads = {}
    corrections = {}
    spectral_features = get_spectral_features_vrad_correction(width=3, smoothing=None)

    for target, data in spectra.items():
        estimations = estimate_v_rad(data, spectral_features)

        # Use median of all the estimations for robustness
        v_rad = np.median(estimations)
        v_rads[target] = v_rad
        corrections[target] = correct_spectrum_for_vrad(data, v_rad)

    return v_rads, corrections


def continuum_normalization(spectra: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """
    Normalize the continuum on a small slice of the spectra
    """
    results = {}
    for target, data in spectra.items():
        normalized, _ = normalize_continuum_slice(
            data, **_fitting_settings.get(target, _fitting_default), w_min=CONTINUUM_SLICE[0], w_max=CONTINUUM_SLICE[1]
        )
        results[target] = normalized

    return results


def calculate_emission(
    normalized_spectra: dict[str, np.ndarray], star_data: pd.DataFrame
) -> tuple[dict[str, float], dict[str, float]]:
    """
    Calculates log R_{HK} and R_o for every target in the spectra.
    NOTE: The project asks to calculate the log of these quantities, but seems that R_hk is always negative, so the log
    can't be calculated. So I'll just leave it as is
    """
    rhk = {target: calculate_log_rhk(data, star_data.loc[target, "b_v"]) for target, data in normalized_spectra.items()}
    ro = {
        target: calculate_log_ro(star_data.loc[target, "prot"], star_data.loc[target, "b_v"])
        for target, data in normalized_spectra.items()
    }
    return rhk, ro


def calculate_vsini(corrected_spectra: dict[str, np.ndarray]) -> dict[str, float]:
    vsini = {}
    # Use a clean iron line as feature for fitting the breoadening profile curve
    feature = 6024.0568
    feature_line_width = 3
    spectrum_slice_width = 100
    slice_hw = spectrum_slice_width / 2
    feature_hw = feature_line_width / 2

    for target, data in corrected_spectra.items():
        normalized_slice, _ = normalize_continuum_slice(
            data=data,
            w_min=feature - slice_hw,
            w_max=feature + slice_hw,
            **_fitting_settings.get(target, _fitting_default),
        )

        normalized_feature = slice_spectrum_data(
            normalized_slice, w_min=feature - feature_hw, w_max=feature + feature_hw
        )

        vsini[target] = estimate_vsini(normalized_feature, feature)

    return vsini


def run():
    print("Running spectroscopy pipeline")

    # Download FIES spectra
    if not _FITS_DATA_FOLDER.exists():
        _FITS_DATA_FOLDER.mkdir(parents=True)
        download_spectra(_FITS_DATA_FOLDER)

    # Save average spectra in a compressed format
    averaged_spectra = average_multiple_spectra()

    # Calculate V_rad from the averaged spectrum of each target
    v_rads, corrected_spectra = vrad_correction(averaged_spectra)

    # Correct spectrum for continuum
    normalized_slices = continuum_normalization(corrected_spectra)

    # Calculate chromospheric excess emission in the CaII H&K lines
    star_data = pd.read_csv(_STAR_DATA, index_col="name")
    rhk, ro = calculate_emission(normalized_slices, star_data)

    # Calculate v sin i
    vsini_estimations = calculate_vsini(corrected_spectra)

    # Aggregate estimations:
    df = pd.DataFrame({"vrad": v_rads, "R_hk": rhk, "Ro": ro, "vsini": vsini_estimations})

    # Create output folder
    _OUT_DATA_FOLDER.mkdir(parents=True, exist_ok=True)

    # Persist results on disk
    print("Persisting results...")
    df.to_csv(_DATAFRAME_FILE, index_label="target")
    df.to_json(_DATAFRAME_JSON_FILE, indent=4, orient="index")

    # Save plots
    save_rhk_ro_plot(rhk, ro, _COMPARISON_PLOT)
    save_chromospheric_lines_plot(corrected_spectra, CaII_EMISSION_LINES, _OUT_DATA_FOLDER)
    save_vsini_comparison_plot(vsini_estimations, star_data, _VSINI_PLOT)

    # Used by jupyter notebooks for further analysis
    np.savez_compressed(_VRAD_CORRECTED_FILE, **corrected_spectra)
    np.savez_compressed(_CONTINUUM_NORMALIZED_FILE, **normalized_slices)

    print(f"Results saved to {_OUT_DATA_FOLDER}")


if __name__ == "__main__":
    run()
