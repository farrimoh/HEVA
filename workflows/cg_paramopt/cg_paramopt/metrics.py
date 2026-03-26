from __future__ import annotations

import logging
import time
from pathlib import Path

from .dependencies import kl_div, np, pd, plt


def load_aa_data(filepath: Path):
    return pd.read_csv(filepath)


def combine_params(df):
    phis = np.sort(np.concatenate([df["p01"], df["p20"], df["p33"], df["p12"]]))
    lengths = np.sort(np.concatenate([df["l0"], df["l1"]]))
    thetas = np.sort(np.concatenate([df["t0"], df["t1"]]))
    return {"l": lengths, "theta": thetas, "phi": phis}


def compute_histogram(data_dict, bins_list):
    result = []
    for i, (column, data) in enumerate(data_dict.items()):
        clean = data[~np.isnan(data)]
        hist, bins = np.histogram(clean, bins=bins_list[i])
        hist = hist / np.sum(hist)
        result.append({"param": column, "bins": bins, "hist": hist})
    return pd.DataFrame(result)


def kl_divergence(p, q):
    return np.sum(np.where(q != 0, kl_div(p, q), 0))


def compute_error(df_cg, df_aa, params, run_dir: Path):
    param_labels = [r"$err_l$", r"$err_\theta$", r"$err_\phi$"]
    plt.figure(figsize=[12, 4])
    errors = []

    for i, param in enumerate(df_cg["param"]):
        plt.subplot(1, 3, i + 1)
        hist_aa = df_aa[df_aa["param"] == param]["hist"].iloc[0]
        hist_cg = df_cg[df_cg["param"] == param]["hist"].iloc[0]
        bins = df_cg[df_cg["param"] == param]["bins"].iloc[0]
        bin_centers = [(bins[j] + bins[j + 1]) / 2 for j in range(len(bins) - 1)]

        err = 100 * kl_divergence(hist_aa, hist_cg)
        errors.append(err)

        plt.plot(bin_centers, hist_aa, label="AA")
        plt.plot(bin_centers, hist_cg, label="CG")
        plt.ylim([0, max(hist_aa) * 1.1])
        plt.yticks([])
        plt.title(f"{param_labels[i]}={err:.3f}")

    total_error = sum(errors)
    plt.suptitle(
        rf"$\kappa_l={params[0]:.0f}, \kappa_\theta={params[1]:.0f}, \kappa_\phi={params[2]:.0f}, error={total_error:.3f}$",
        y=0.98,
    )
    plt.tight_layout()
    timestamp = time.strftime("%H%M%S")
    filename = run_dir / f"{timestamp}_kappa_l_{params[0]:.0f}_theta_{params[1]:.0f}_phi_{params[2]:.0f}.png"
    plt.savefig(filename)
    logging.info("Saved plot: %s", filename)
    plt.close()
    return total_error

