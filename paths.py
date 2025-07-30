import os
from pathlib import Path

_CURRENT_DIR = Path(os.path.realpath(__file__)).parent
_DATA_FOLDER = _CURRENT_DIR / "data"

FITS_DATA_FOLDER = _DATA_FOLDER / "fits"
OUT_DATA_FOLDER = _DATA_FOLDER / "out"
OBSERVATIONS_FILE = _DATA_FOLDER / "observations.csv"
STAR_DATA = _DATA_FOLDER / "star_data.csv"

# Output files
AVERAGED_FILE = OUT_DATA_FOLDER / "1_averaged_spectra.npz"
VRAD_CORRECTED_FILE = OUT_DATA_FOLDER / "2_vrad_corrected_spectra.npz"
CONTINUUM_NORMALIZED_FILE = OUT_DATA_FOLDER / "3_continuum_corrected_spectra.npz"
COMPARISON_PLOT = OUT_DATA_FOLDER / "rhk_ro_comparison.png"
VSINI_PLOT = OUT_DATA_FOLDER / "vsini_comparison.png"
DATAFRAME_FILE = OUT_DATA_FOLDER / "results.csv"
DATAFRAME_JSON_FILE = OUT_DATA_FOLDER / "results.json"

FIES_REPO_URL = "https://kirnis.kapsi.fi/NOT2023/data/"
