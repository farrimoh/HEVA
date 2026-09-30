#!/usr/bin/env python3
"""One figure per structure: AA distributions of lengths (Å) and angles (rad), per contact and pooled."""
import csv, json, math
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import argparse
from pathlib import Path as _Path
import sys as _sys
_sys.path.insert(0, str(_Path(__file__).resolve().parent))
import _paths

_ap = argparse.ArgumentParser(description=__doc__)
_paths.add_arguments(_ap, "analysis/per_structure")
_P = _paths.resolve(_ap.parse_args(), need=("patches", "capsid"))

A = _P["observables"]
OUT = _P["out"]
STRUCTS = {"Cp3": ("patches/Cp3_observables.csv", "Cp3 triangle (AB,CD,BA)"),
           "Cp5_curved_as_capsid": ("patches/Cp5_curved_as_capsid_observables.csv", "Cp5 diamond, capsid-like"),
           "Cp5_too_curved": ("patches/Cp5_too_curved_observables.csv", "Cp5 diamond, too curved"),
           "Cp10": ("patches/Cp10_observables.csv", "Cp10 fivefold patch"),
           "Cp12": ("patches/Cp12_observables.csv", "Cp12 sixfold patch"),
           "capsid": ("capsid/capsid_T4_stride100frames_observables.csv", "Full T4 capsid (500 frames)")}
OBS = ["l0", "l1", "t0", "t1", "p33", "p12", "p01", "p20"]
TITLES = {"l0": "CD dimer length (Å)", "l1": "AB dimer length (Å)", "t0": "CD dihedral t0 (rad, signed)", "t1": "AB dihedral t1 (rad, signed)",
          "p33": "corner DC→DC p33 (rad)", "p12": "corner BA→AB p12 (rad)", "p01": "corner CD→BA p01 (rad)", "p20": "corner AB→CD p20 (rad)"}
CAPSID_REF = {"l0": 90.72, "l1": 81.89, "t0": 0.240, "t1": 0.476, "p33": 1.047, "p12": 1.174, "p01": 0.984, "p20": 0.984}
TAG_COLOR = {0: "#1c5cab", 1: "#5598e7", 2: "#9ec5f4", 3: "#cde2fb"}

def load(path):
    per = {}
    with open(path) as fh:
        r = csv.DictReader(fh)
        has_signed = "signed_value" in r.fieldnames; has_tag = "weak_junctions" in r.fieldnames
        for row in r:
            o = row["observable"]
            v = float(row["signed_value"]) if (has_signed and o.startswith("t")) else float(row["value"])
            tag = int(row["weak_junctions"]) if has_tag else 0
            per.setdefault(o, {}).setdefault(row["contact_label"], [tag, []])[1].append(v)
    return {o: {c: (t, np.array(v)) for c, (t, v) in d.items()} for o, d in per.items()}

summary = {}
for key, (rel, title) in STRUCTS.items():
    data = load(A / rel)
    fig, axes = plt.subplots(2, 4, figsize=(17, 7.4), dpi=140)
    summary[key] = {}
    for ax, o in zip(axes.ravel(), OBS):
        if o not in data:
            ax.axis("off"); ax.text(0.5, 0.5, f"{TITLES[o]}\nnot present in this structure", ha="center", va="center", fontsize=9, color="#52514e", transform=ax.transAxes); continue
        contacts = data[o]
        pooled = np.concatenate([v for _, v in contacts.values()])
        lo, hi = np.percentile(pooled, [0.2, 99.8]); pad = 0.08 * (hi - lo); bins = np.linspace(lo - pad, hi + pad, 70)
        many = len(contacts) > 12
        for c, (tag, v) in sorted(contacts.items()):
            h, e = np.histogram(v, bins=bins, density=True)
            ax.plot(0.5 * (e[1:] + e[:-1]), h, color=TAG_COLOR.get(tag, "#cde2fb"), lw=(0.6 if many else 1.1), alpha=(0.35 if many else 0.9))
        h, e = np.histogram(pooled, bins=bins, density=True)
        ax.plot(0.5 * (e[1:] + e[:-1]), h, color="#0b0b0b", lw=2.0, label="pooled")
        if key != "capsid":
            ax.axvline(CAPSID_REF[o], color="#52514e", ls="--", lw=1.0)
        m, s = pooled.mean(), pooled.std()
        cm = [v.mean() for _, v in contacts.values()]
        within = math.sqrt(np.mean([v.var() for _, v in contacts.values()]))
        between = float(np.std(cm))
        tags = sorted({t for t, _ in contacts.values()})
        ax.set_title(f"{TITLES[o]}\nmean {m:.3f}  sd {s:.3f}  |  {len(contacts)} contacts, tags {tags}\nwithin-contact sd {within:.3f}, between-contact sd {between:.3f}", fontsize=8, loc="left")
        ax.set_yticks([])
        for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
        ax.grid(axis="x", color="#e5e4df", lw=0.6); ax.set_axisbelow(True)
        summary[key][o] = {"mean": float(m), "sd": float(s), "within_sd": within, "between_sd": between,
                           "contacts": {c: {"tag": t, "mean": float(v.mean()), "sd": float(v.std())} for c, (t, v) in contacts.items()}}
    fig.suptitle(f"{title}: AA distributions per contact (blue shades by weak-junction tag: dark = 0, lighter = 1, 2) and pooled (black). Dashed grey: capsid mean.", fontsize=10, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT / f"{key}.png", facecolor="#fcfcfb"); plt.close(fig)
    print("wrote", key)
(OUT / "per_contact_summary.json").write_text(json.dumps(summary, indent=1))
