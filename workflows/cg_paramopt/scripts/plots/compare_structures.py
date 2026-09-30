#!/usr/bin/env python3
"""Compare HEVA-style lengths, dihedrals and corner angles across AA patch trajectories."""
import csv, json, math, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import argparse
from pathlib import Path as _Path
import sys as _sys
_sys.path.insert(0, str(_Path(__file__).resolve().parent))
import _paths

_ap = argparse.ArgumentParser(description=__doc__)
_paths.add_arguments(_ap, "analysis")
_P = _paths.resolve(_ap.parse_args(), need=("patches", "capsid"))

ROOT = _P["observables"]
OUTDIR = _P["out"]
STRUCTS = ["Cp3", "Cp3_symm", "Cp5_flat", "Cp5_curved_as_capsid", "Cp5_too_curved", "Cp10", "Cp12"]
PLOT_STRUCTS = ["Cp3", "Cp5_curved_as_capsid", "Cp5_too_curved", "Cp10", "Cp12", "capsid"]
CAPSID_CSV = ROOT / "capsid" / "capsid_T4_stride100frames_observables.csv"
LABELS = {"capsid": "full T4 capsid", "Cp3": "Cp3 (AB,CD,BA)", "Cp3_symm": "Cp3 (DC,DC,DC)", "Cp5_flat": "Cp5 flat", "Cp5_curved_as_capsid": "Cp5 capsid-like",
          "Cp5_too_curved": "Cp5 too curved", "Cp10": "Cp10 fivefold", "Cp12": "Cp12 sixfold"}
COLORS = dict(zip(STRUCTS, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7"]))
COLORS["capsid"] = "#0b0b0b"
OBS = ["l0", "l1", "t0", "t1", "p33", "p12", "p01", "p20"]
TITLES = {"l0": "CD dimer length (Å)", "l1": "AB dimer length (Å)", "t0": "CD dihedral t0 (rad)", "t1": "AB dihedral t1 (rad)",
          "p33": "corner DC→DC p33 (rad)", "p12": "corner BA→AB p12 (rad)", "p01": "corner CD→BA p01 (rad)", "p20": "corner AB→CD p20 (rad)"}
LEGACY = {"l0": (1.0512, 0.0115), "l1": (0.9488, 0.0125), "t0": (0.2395, 0.0440), "t1": (0.4760, 0.0433),
          "p33": (1.0472, 0.0154), "p12": (1.1742, 0.0199), "p01": (0.9837, 0.0175), "p20": (0.9837, 0.0161)}
EST = ["patches", "calibrated", "centroid", "circumcenter"]
# figure estimator per structure: calibrated only where every incomplete vertex has 3 members (robust); centroid otherwise
FIG_EST = {s: "patches" for s in STRUCTS}

def load(est, s):
    p = CAPSID_CSV if s == "capsid" else ROOT / est / f"{s}_observables.csv"
    if not p.exists():
        return None
    data = {}
    frames = {}
    with p.open() as fh:
        for r in csv.DictReader(fh):
            data.setdefault(r["observable"], []).append(float(r["value"]))
            frames.setdefault(r["observable"], []).append(int(r["frame"]))
            if r["observable"].startswith("t"):
                data.setdefault(r["observable"] + "_class" + r["energy_class"], []).append(float(r["value"]))
    return {k: np.array(v) for k, v in data.items()}, {k: np.array(v) for k, v in frames.items()}

def normalized_lengths(d, fr):
    """Per-frame normalization by (mean l0 + mean l1)/2, as the importer does."""
    if "l0" not in d or "l1" not in d:
        return None
    nf = max(fr["l0"].max(), fr["l1"].max()) + 1
    m0 = np.bincount(fr["l0"], weights=d["l0"], minlength=nf) / np.bincount(fr["l0"], minlength=nf)
    m1 = np.bincount(fr["l1"], weights=d["l1"], minlength=nf) / np.bincount(fr["l1"], minlength=nf)
    scale = (m0 + m1) / 2
    return d["l0"] / scale[fr["l0"]], d["l1"] / scale[fr["l1"]], scale

def block_sem(vals, frames, nblocks=10):
    nf = frames.max() + 1
    edges = np.linspace(0, nf, nblocks + 1)
    means = [vals[(frames >= a) & (frames < b)].mean() for a, b in zip(edges[:-1], edges[1:])]
    return float(np.std(means, ddof=1) / math.sqrt(nblocks)), means

all_data = {est: {s: load(est, s) for s in STRUCTS} for est in EST}
for est in EST:
    all_data[est]["capsid"] = load(est, "capsid")
FIG_EST["capsid"] = "patches"
rows = []
table = {}
for s in STRUCTS:
    cal = all_data["patches"][s]
    if cal is None:
        continue
    d, fr = cal
    norm = normalized_lengths(d, fr)
    for o in OBS:
        if o not in d:
            continue
        sem, blocks = block_sem(d[o], fr[o])
        per_frame = len(d[o]) / (fr[o].max() + 1)
        alt = {e: (all_data[e][s][0][o].mean() if all_data[e][s] and o in all_data[e][s][0] else float("nan")) for e in EST[1:]}
        altsd_cal = None
        altsd = {e: (all_data[e][s][0][o].std() if all_data[e][s] and o in all_data[e][s][0] else float("nan")) for e in EST[1:]}
        extra = ""
        if o in ("l0", "l1") and norm is not None:
            nv = norm[0] if o == "l0" else norm[1]
            extra = f"{nv.mean():.3f} ± {nv.std():.3f}"
        classes = sorted(k for k in d if k.startswith(o + "_class"))
        cls = ",".join(f"{k.split('class')[1]}:{len(d[k]) // (fr[o].max() + 1)}" for k in classes)
        rows.append([s, o, f"{d[o].mean():.3f}", f"{d[o].std():.3f}", f"{sem:.3f}", f"{per_frame:.0f}", extra,
                     f"{alt['calibrated']:.3f}/{altsd['calibrated']:.3f}", f"{alt['centroid']:.3f}/{altsd['centroid']:.3f}", cls,
                     " ".join(f"{b:.3f}" for b in blocks[::3])])
        table[(s, o)] = (d[o].mean(), d[o].std(), sem)

with (OUTDIR / "comparison_table.csv").open("w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["structure", "observable", "mean_capsid_calibrated", "sd_capsid_calibrated", "sem_10blocks", "contacts_per_frame",
                "frame_normalized_length_mean_sd", "regular_polygon_mean/sd", "centroid_mean/sd", "dihedral_energy_class:count", "block_means_subset"])
    w.writerows(rows)

# Figure: one panel per observable, per-structure density curves, legacy capsid mean as dashed reference (angles only)
fig, axes = plt.subplots(2, 4, figsize=(16, 7.2), dpi=150)
plt.rcParams["font.size"] = 9
for ax, o in zip(axes.ravel(), OBS):
    handles = []
    for s in PLOT_STRUCTS:
        cal = all_data[FIG_EST[s]][s]
        if cal is None or o not in cal[0]:
            continue
        v = cal[0][o]
        lo, hi = np.percentile(v, [0.5, 99.5])
        bins = np.linspace(lo, hi, 80)
        h, e = np.histogram(v, bins=bins, density=True)
        c = 0.5 * (e[1:] + e[:-1])
        ax.plot(c, h, color=COLORS[s], lw=(2.2 if s == "capsid" else 1.6), label=LABELS[s])
    ax.set_title(TITLES[o], fontsize=10, loc="left")
    ax.set_yticks([])
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.grid(axis="x", color="#e5e4df", lw=0.6)
    ax.set_axisbelow(True)
h, l = axes[0, 0].get_legend_handles_labels()
seen = {}
for ax in axes.ravel():
    for hh, ll in zip(*ax.get_legend_handles_labels()):
        seen.setdefault(ll, hh)
fig.legend(seen.values(), seen.keys(), loc="lower center", ncol=6, frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.01))
fig.suptitle("AA intermediates vs full capsid: HEVA-style observables from residue-132 Cα landmarks (patch junctions corrected with capsid-measured offsets)",
             fontsize=11, x=0.01, ha="left")
fig.tight_layout(rect=(0, 0.05, 1, 0.96))
fig.savefig(OUTDIR / "observable_distributions.png", facecolor="#fcfcfb")

# second figure: dihedral energy-class split and lengths by structure (Å) with legacy ratio
fig2, ax2 = plt.subplots(1, 2, figsize=(12, 4.2), dpi=150)
names = [s for s in PLOT_STRUCTS if all_data[FIG_EST[s]][s]]
x = np.arange(len(names))
for k, (o, col) in enumerate([("l0", "#eb6834"), ("l1", "#2a78d6")]):
    means = [all_data[FIG_EST[s]][s][0][o].mean() if o in all_data[FIG_EST[s]][s][0] else np.nan for s in names]
    sds = [all_data[FIG_EST[s]][s][0][o].std() if o in all_data[FIG_EST[s]][s][0] else np.nan for s in names]
    ax2[0].errorbar(x + (k - 0.5) * 0.18, means, yerr=sds, fmt="o", color=col, ms=5, capsize=3, lw=1.2, label=f"{o} ({'CD' if o=='l0' else 'AB'} dimer)")
ax2[0].set_xticks(x); ax2[0].set_xticklabels([LABELS[s] for s in names], rotation=30, ha="right")
ax2[0].set_ylabel("vertex-to-vertex length (Å), mean ± SD"); ax2[0].legend(frameon=False)
ax2[0].set_title("Dimer lengths", loc="left", fontsize=10)
for o, col, mk in [("t0", "#eb6834", "s"), ("t1", "#2a78d6", "o")]:
    means = [all_data[FIG_EST[s]][s][0][o].mean() if o in all_data[FIG_EST[s]][s][0] else np.nan for s in names]
    sds = [all_data[FIG_EST[s]][s][0][o].std() if o in all_data[FIG_EST[s]][s][0] else np.nan for s in names]
    ax2[1].errorbar(x + (0.09 if o == "t1" else -0.09), means, yerr=sds, fmt=mk, color=col, ms=5, capsize=3, lw=1.2, label=f"{o} ({'CD' if o=='t0' else 'AB'} internal dimer)")
ax2[1].set_xticks(x); ax2[1].set_xticklabels([LABELS[s] for s in names], rotation=30, ha="right")
ax2[1].set_ylabel("dihedral between adjacent faces (rad), mean ± SD"); ax2[1].legend(frameon=False)
ax2[1].set_title("Dihedral angles", loc="left", fontsize=10)
for ax in ax2:
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.grid(axis="y", color="#e5e4df", lw=0.6); ax.set_axisbelow(True)
fig2.tight_layout()
fig2.savefig(OUTDIR / "lengths_and_dihedrals.png", facecolor="#fcfcfb")

print(f"{'structure':22s} {'obs':4s} {'mean':>7s} {'sd':>6s} {'sem':>6s} {'n/fr':>4s} {'norm len':>14s} {'reg-polygon':>14s} {'centroid':>14s} {'class':>8s}")
for r in rows:
    print(f"{r[0]:22s} {r[1]:4s} {r[2]:>7s} {r[3]:>6s} {r[4]:>6s} {r[5]:>4s} {r[6]:>14s} {r[7]:>14s} {r[8]:>14s} {r[9]:>8s}")
