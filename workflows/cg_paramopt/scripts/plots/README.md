# AA patch extraction and distribution plots

Maps all-atom (AA) trajectories of HBV capsid intermediates onto HEVA's coarse-grained
observables, then plots their distributions against the full T4 capsid.

Inputs are AA trajectories with one landmark atom per protein: the residue-132 Cα of each
core-protein monomer. Segment IDs are the chain letter plus the asymmetric-unit index, so
dimers are A_i–B_i and C_i–D_i. Structures used so far: Cp3 triangles, the three Cp5
diamonds, the Cp10 fivefold and Cp12 sixfold patches, and the full T4 capsid.

## 1. Extract observables

`../extract_patch_observables.py` clusters landmarks into CG vertices, builds the patch
topology from segment IDs, orients faces to HEVA's face cycles, and writes every observable
with the same definitions as HEVA's `long_v1` sampler: dimer lengths (`l0` CD, `l1` AB),
dihedrals between adjacent faces (`t0`, `t1`), and in-face corner angles (`p33`, `p12`, `p01`,
`p20`). On the full capsid it reproduces the 2022 AA reference table.

A junction in an isolated patch usually has only two or three of its five or six proteins.
The `capsid` estimator places those junctions with offsets measured on the complete rings of
the full capsid, stored in `../../fixtures/aa_mapping/vertex_offset_calibration_T4.json` and
produced by `../calibrate_vertex_offsets.py`.

```bash
python ../extract_patch_observables.py --psf Cp10_protein_CA_1us.psf --dcd Cp10_protein_CA_1us.dcd \
  --name Cp10 --output-dir DATA/observables/patches --stride 5 --dt-ps 20 \
  --estimator capsid --calibration ../../fixtures/aa_mapping/vertex_offset_calibration_T4.json
# the complete capsid needs no correction
python ../extract_patch_observables.py --psf Ca132.psf --dcd prod.1us.stride100frames.Ca132.dcd \
  --name capsid_T4_stride100frames --output-dir DATA/observables/capsid --stride 1 --estimator centroid
```

Use `--relaxed-faces` for the CD-only flat diamond, which is not a T4 arrangement. Each run
writes `<name>_observables.csv` and `<name>_summary.json`. The CSV is long format, one row per
contact per frame:
`frame, source_frame, observable, contact_id, half_edge_type, observable_class, energy_class,
value, contact_label, signed_value, weak_junctions, end_tags`. Lengths are in Å and angles in
radians. `signed_value` carries the dihedral sign. `weak_junctions` counts the two-protein
junctions a contact touches; contacts with 0 are measured without any capsid-derived
correction. `end_tags` gives the occupancy of the junctions involved, for example `P3/5-H3/6`
for a dimer from a pentamer with 3 of 5 proteins present to a hexamer with 3 of 6.

`<name>_summary.json` records the mapping itself: which proteins form each CG junction, ring
size and occupancy, the oriented faces, the boundary half-edges, and the extraction settings.

## 2. Plot distributions

Both scripts read one data directory, via `--data-dir` or `HEVA_AA_DATA`:

```
DATA/observables/patches/   patch CSVs from step 1
DATA/observables/capsid/     capsid CSV from step 1
```

```bash
export HEVA_AA_DATA=/path/to/DATA
python compare_structures.py            # all structures overlaid per observable + comparison table
python per_structure_distributions.py   # one figure per structure: per-contact and pooled distributions
```

Figures go to `DATA/analysis/` unless `--out-dir` is given. If the optional folders
`observables/calibrated`, `centroid` or `circumcenter` exist, `compare_structures.py` also
fills the estimator-comparison columns of its table; otherwise those columns are `nan`.

Requirements: numpy, scipy, matplotlib, MDAnalysis (see `../../requirements.txt`).
