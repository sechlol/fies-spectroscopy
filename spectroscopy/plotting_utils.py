from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from spectroscopy.fies_utils import smooth_spectrum, slice_spectrum_data
from spectroscopy.spectral_features import SpectralFeature
from spectroscopy.spectrum_loader import Spectrum


def save_rhk_ro_plot(rhk: dict[str, float], ro: dict[str, float], out_file: Path):
    plt.figure(figsize=(10, 6))

    for target in rhk.keys():
        plt.annotate(target, (ro[target], rhk[target]), textcoords="offset points", xytext=(0, 5), ha="center")

    plt.scatter(list(ro.values()), list(rhk.values()))
    plt.xlabel(r"$\log_{10}{R_o}$")
    plt.ylabel(r"$\log_{10}{R'_{HK}}$")
    plt.title(r"Scatter Plot of $R'_{HK}$ vs $R_o$")
    plt.tight_layout()
    plt.savefig(out_file)


def save_chromospheric_lines_plot(
    vrad_corrected_spectra: dict[str, np.ndarray],
    emission_lines: dict[str, tuple[str, float]],
    out_path: Path,
):
    for target, data in vrad_corrected_spectra.items():
        fig, axes = plt.subplots(2, 3, figsize=(10, 6))
        for ax, (line_name, w) in zip(axes.flatten(), emission_lines.items()):
            feature_slice = slice_spectrum_data(data, w_min=w - 1, w_max=w + 1)
            ax.plot(feature_slice[:, 0], feature_slice[:, 1])
            ax.axvline(w, color="red")
            ax.set(title=line_name, xlabel=r"$\lambda [Å]$")

        fig.suptitle(f"{target} Chromospheric line profiles")
        fig.tight_layout()
        plt.savefig(out_path / f"chromo_profiles_{target}.png")


def save_vsini_comparison_plot(vsini_estimates: dict[str, float], star_data: pd.DataFrame, out_file: Path):
    targets = list(vsini_estimates.keys())
    x_values = np.arange(len(targets))
    vsini_values = np.array(list(vsini_estimates.values()))

    plt.figure(figsize=(10, 6))
    plt.scatter(x_values, vsini_values, label="Estimate", color="red")

    vsini_min_values = star_data.loc[targets, "vsini_min"]
    vsini_max_values = star_data.loc[targets, "vsini_max"]
    bar_lengths = (vsini_max_values - vsini_min_values) / 2
    mid_points = (vsini_max_values + vsini_min_values) / 2
    plt.errorbar(
        x_values,
        mid_points,
        yerr=[bar_lengths, bar_lengths],
        label="Literature values",
        fmt="none",
        color="blue",
        capsize=4,
    )

    plt.title("Estimated vsini")
    plt.xlabel("Target")
    plt.ylabel("v sin(i) [km/s]")
    plt.xticks(x_values, targets, rotation=45, ha="right")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_file)


def plot_zoom(spectra: list[Spectrum], central_wavelengths: list[float], width: float, smoothing: Optional[int] = None):
    rows = len(spectra)
    cols = len(central_wavelengths)
    half_width = width / 2

    fig, axes = plt.subplots(rows, cols, figsize=(cols * 10, rows * 4))
    axes = axes.flatten()
    for i, spectrum in enumerate(spectra):
        for j, w in enumerate(central_wavelengths):
            data = smooth_spectrum(spectrum.data, smoothing) if smoothing else spectrum.data
            interval = data[np.abs(data[:, 0] - w) <= half_width]
            axes[i * cols + j].axvline(w, color="red")
            axes[i * cols + j].plot(interval[:, 0], interval[:, 1])
            axes[i * cols + j].set(
                title=rf"{spectrum.target_name} ({spectrum.id}): $\lambda = {w} \pm {half_width:.3f}$"
            )
    plt.tight_layout()
    plt.show()


def plot_spectral_features(
    spectrum: Spectrum, spectral_features: list[SpectralFeature], smoothing: Optional[int] = None
):
    rows = (len(spectral_features) + 1) // 2
    fig, axes = plt.subplots(rows, 2, figsize=(15, rows * 4))
    for i, (feature, ax) in enumerate(zip(spectral_features, axes.flatten())):
        smoothing = smoothing or feature.smoothing
        data = smooth_spectrum(spectrum.data, periods=smoothing) if smoothing else spectrum.data
        interval = data[np.abs(data[:, 0] - feature.w) < feature.half_width]
        x = interval[:, 0]
        mean_residual = np.abs(interval[::-1][:, 1] - interval[:, 1]).mean()
        ax.axvline(feature.w, color="red", label=f"w = {feature.w}")
        ax.plot(x, interval[:, 1], label="Original")
        ax.plot(x, interval[::-1][:, 1], label="Specular")
        ax.set(title=rf"{feature.name} - Residual $\mu = {mean_residual:.5f}$", xlabel=r"$\lambda [Å]$")
        ax.legend()

    plt.suptitle(f"{spectrum.target_name} - ({spectrum.id}) - spectral features")
    plt.tight_layout()
    plt.show()


def plot_raw_spectra(spectra: list[Spectrum]):
    n = len(spectra)
    fig, axes = plt.subplots(n, 1, figsize=(15, n * 3))
    axes = axes.flatten() if n > 1 else [axes]

    for spectrum, ax in zip(spectra, axes):
        ax.plot(spectrum.data[:, 0], spectrum.data[:, 1])
        ax.set(title=f"{spectrum.target_name} - {spectrum.id}", xlabel=r"$\lambda [Å]$")

    plt.tight_layout()
    plt.show()


def plot_residuals_from_mean(spectra: list[Spectrum]):
    epsilon = 1e-6
    all_spectra = np.vstack([s.data[:, 1] for s in spectra])
    residuals = np.abs(all_spectra - all_spectra.mean(axis=0))

    # epsilon needed to avoid division by 0
    rel_errors = residuals / (all_spectra + epsilon)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 8))
    for i, spectrum in enumerate(spectra):
        ax1.plot(spectrum.data[:, 0], residuals[i, :], label=spectrum.id, alpha=0.5)
        ax2.plot(spectrum.data[:, 0], rel_errors[i, :], label=spectrum.id, alpha=0.5)
        ax1.legend()
        ax2.legend()

    ax1.set(title=f"{spectra[0].target_name} Residuals from the mean spectrum", xlabel=r"$\lambda [Å]$")
    ax2.set(
        title=f"{spectra[0].target_name} Relative errors from the mean spectrum", xlabel=r"$\lambda [Å]$", ylim=(-2, 2)
    )
    plt.tight_layout()
    plt.show()
