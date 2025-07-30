import numpy as np

from typing import Optional
from scipy import optimize
from sklearn.preprocessing import minmax_scale

from spectroscopy.spectral_features import SpectralFeature

C = 299792.458  # km/s


def slice_spectrum_data(spectrum_data: np.ndarray, w_min, w_max) -> np.ndarray:
    return spectrum_data[(w_min < spectrum_data[:, 0]) & (spectrum_data[:, 0] < w_max)]


def estimate_v_rad(
    spectrum: np.ndarray, spectral_features: list[SpectralFeature], smoothing: Optional[int] = None
) -> np.ndarray:
    estimations = []
    for feature in spectral_features:
        smoothing = smoothing or feature.smoothing
        data = smooth_spectrum(spectrum, smoothing) if smoothing else spectrum
        result = optimize.minimize_scalar(
            _calculate_loss,
            bracket=(-feature.half_width, feature.half_width),
            bounds=(-feature.half_width, feature.half_width),
            args=(data, feature),
            method="bounded",
        )
        delta_f = result.x
        v_rad = delta_f * C / feature.w  # km/s
        estimations.append(v_rad)
    return np.array(estimations)


def correct_spectrum_for_vrad(spectrum: np.ndarray, v_rad: float) -> np.ndarray:
    w_corrected = spectrum[:, 0] + spectrum[:, 0] * v_rad / C
    return np.array([w_corrected, spectrum[:, 1]]).T


def normalize_continuum_slice(
    data: np.ndarray,
    degree: int,
    iterations: Optional[int] = 4,
    w_min: Optional[float] = 0,
    w_max: Optional[float] = np.inf,
) -> Tuple[np.ndarray, np.ndarray]:
    data_slice = slice_spectrum_data(data, w_min, w_max)
    return normalize_continuum(data_slice, degree=degree, iterations=iterations)


def normalize_continuum(data: np.ndarray, degree: int, iterations: Optional[int] = 4) -> Tuple[np.ndarray, np.ndarray]:
    x_full = data[:, 0]
    y_full = data[:, 1]
    x = x_full.copy()
    y = y_full.copy()
    fits = []

    # Continuum normalisation with an iterative linear fit
    # storing intermediate continuum fits in conts
    coefficients = np.zeros(degree)
    for _ in range(iterations):
        coefficients = np.polyfit(x, y, deg=degree)
        y_fit_full = np.polyval(coefficients, x_full)
        fits.append(y_fit_full)

        y_fit_subset = np.polyval(coefficients, x)
        y_normalized = y / y_fit_subset
        i_above_mean = y_normalized > np.mean(y_normalized)
        y = y[i_above_mean]
        x = x[i_above_mean]

    final_fit = np.polyval(coefficients, x_full)
    y_normalized = y_full / final_fit
    return np.vstack([x_full, y_normalized]).T, np.array(fits)


def _rotational_profile(x: np.ndarray, x0: float, vsini: float) -> np.ndarray:
    """
    Reference paper:
        - Accurate stellar rotational velocities using the Fourier transform of the cross correlation maximum
        - (Diaz et al)
        - https://www.aanda.org/articles/aa/pdf/2011/07/aa16386-10.pdf
    """
    epsilon = 0.6

    xl = np.log(x / x0) * C / vsini
    i = xl**2 < 1
    xi = xl[i]

    g = np.zeros_like(x)
    g1 = 2 * (1 - epsilon) * np.sqrt((1 - xi**2)) + np.pi * epsilon / 2 * (1 - xi**2)
    g2 = np.pi * (1 - epsilon / 3)
    g[i] = g1 / g2

    return 1 - g


def estimate_vsini(normalized_feature: np.ndarray, central_feature: float) -> float:
    y_minmax = minmax_scale(normalized_feature[:, 1])

    def loss_function(vsini: float):
        g = _rotational_profile(normalized_feature[:, 0], central_feature, vsini)
        mae = np.abs(minmax_scale(g) - y_minmax).sum()
        return mae

    fit_result = optimize.minimize(loss_function, x0=np.array([10]))
    return fit_result.x[0]


def _calculate_loss(delta_lambda: float, spectrum: np.ndarray, feature: SpectralFeature) -> float:
    interval = spectrum[np.abs(spectrum[:, 0] - feature.w + delta_lambda) <= feature.half_width][:, 1]
    residuals = np.abs(interval[::-1] - interval)
    # print(f"F = {feature.w}+{delta_lambda}, mean_residual: {residuals.mean()}")
    return residuals.mean()


def smooth_spectrum(spectrum: np.ndarray, periods: float):
    smoothed = spectrum.copy()
    smoothed[:, 1] = np.convolve(smoothed[:, 1], np.full(shape=periods, fill_value=1 / periods), mode="same")
    return smoothed
