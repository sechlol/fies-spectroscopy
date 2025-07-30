from dataclasses import dataclass
from typing import Optional

import numpy as np

_SPECTRAL_LINES_VRAD = np.array([5633.9458, 6024.0568, 6411.6476, 6430.8446, 6439.0750, 6643.6304])
_SPECTRAL_LINES_NAMES_VRAD = ["Fe1", "Fe1", "Fe1", "Fe1", "Ca1", "Ni1"]

CaII_EMISSION_BANDS = {
    "V": (3891, 3911),
    "R": (3991, 4011),
    "H": (3968.47 - 1.09, 3968.47 + 1.09),
    "K": (3933.66 - 1.09, 3933.66 + 1.09),
}

CaII_EMISSION_LINES = {
    "H": 3968.47,
    "K": 3933.66,
    "Ha": 6562.8,
    "IRT_I": 8498,
    "IRT_R": 8542,
    "IRT_T": 8662,
}

CONTINUUM_SLICE = (3875, 4025)


@dataclass
class SpectralFeature:
    name: str
    w: float
    width: float
    smoothing: Optional[int] = None

    @property
    def half_width(self) -> float:
        return self.width / 2


def get_spectral_features_vrad_correction(width: Optional[float] = 3, smoothing: Optional[int] = None):
    return [
        SpectralFeature(name=n, w=f, width=df, smoothing=smoothing)
        for (n, f, df) in zip(
            _SPECTRAL_LINES_NAMES_VRAD, _SPECTRAL_LINES_VRAD, np.full_like(_SPECTRAL_LINES_VRAD, width)
        )
    ]


def calculate_log_rhk(data: np.ndarray, b_v: float) -> float:
    integrals = _get_integrals_of_hk_bands(data)
    s_index = 19.76 * (integrals["H"] + integrals["K"]) / (integrals["R"] + integrals["V"])
    cfc = 10 ** (0.25 * b_v**3 - 1.33 * b_v**2 + 0.43 * b_v + 0.24)
    r_phot = 10 ** (-4.898 + 1.918 * b_v**2 - 2.893 * b_v**3)
    r_hk_raw = 1.34 * (10**-4) * cfc * s_index
    r_hk_corrected = r_hk_raw - r_phot
    return np.log10(r_hk_corrected)


def calculate_log_ro(b_v: float, p_rot: float) -> float:
    x = 1 - b_v
    if x > 0:
        log_tc = 1.362 - 0.166 * x + 0.025 * x**2 - 5.323 * x**3
    else:
        log_tc = 1.362 - 0.14 * x
    return np.log10(p_rot * 10**log_tc)


def _get_triangular_function(n: int):
    a = 0
    b = n
    center = (a + b + 1) / 2
    x = np.arange(b)
    y = np.ones(b)

    slope_left = 1 / (center - a)
    slope_right = 1 / (center - b)

    y[x < center] = x[x < center] * slope_left
    y[x > center] = x[x > center] * slope_right + 2
    return y


def _get_integrals_of_hk_bands(data: np.ndarray):
    integrals = {}

    # For V and R bands, just calculate the integral of the (rectangular) curve
    for line in ["V", "R"]:
        x1, x2 = CaII_EMISSION_BANDS[line]
        i = (x1 < data[:, 0]) & (data[:, 0] < x2)
        integrals[line] = np.trapz(y=data[i, 1], x=data[i, 0])

    # H and K need to be multiplied by a triangular function
    for line in ["H", "K"]:
        x1, x2 = CaII_EMISSION_BANDS[line]
        i = (x1 <= data[:, 0]) & (data[:, 0] <= x2)
        x = data[i, 0]

        # For H and K, multiply the curve by a triangular function
        y_triangular_normalizer = _get_triangular_function(len(x))
        y = data[i, 1] * y_triangular_normalizer
        integrals[line] = np.trapz(y=y, x=x)

    return integrals
