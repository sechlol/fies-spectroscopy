from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits


@dataclass
class Spectrum:
    id: str
    target_name: str
    data: np.ndarray


class SpectrumLoader:
    def __init__(self, fits_path: Path, observations_path: Path):
        self._data_path = fits_path
        self._observations = self._load_observations(observations_path)
        self._loaded_spectra = {}

    @property
    def target_list(self) -> list[str]:
        return list(self._observations.keys())

    def get_available_spectra_count(self, target_name) -> int:
        return len(self._observations[target_name])

    def get_spectra(self, target_name: str) -> list[Spectrum]:
        if target_name not in self._loaded_spectra:
            spectra = [Spectrum(sid, target_name, self._read_1D_fits(sid)) for sid in self._observations[target_name]]

            # Ensures that all spectra have the same length
            min_size = np.min([len(s.data) for s in spectra])
            for s in spectra:
                s.data = s.data[:min_size]
            self._loaded_spectra[target_name] = spectra

        return self._loaded_spectra[target_name]

    @staticmethod
    def _load_observations(file_path: Path) -> dict[str, list[str]]:
        df = pd.read_csv(file_path)
        return df.groupby("target_name")["id"].apply(list).to_dict()

    def _read_1D_fits(self, spectrum_id: str) -> np.ndarray:
        path = self._data_path.joinpath(f"{spectrum_id}_step011_merge.fits")
        with fits.open(path) as hdu_list:
            data = hdu_list[0].data
            assert data.ndim == 1, f"Fits file has more than 1 dimension"

            header = hdu_list[0].header
            naxis = header["NAXIS1"]
            cdelt = header["CDELT1"]
            crval = header["CRVAL1"]
            # crpix = header['CRPIX1']
            # rawfile = header['FILENAME']
            # starname = header['OBJECT']

            # Create wl grid
            wl = np.arange(naxis)
            wl = wl * cdelt + crval

            # Correct negative outliers
            data[data < 0] = 0

            return np.array([wl, data]).T
