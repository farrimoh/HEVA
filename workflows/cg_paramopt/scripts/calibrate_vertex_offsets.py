#!/usr/bin/env python3
"""Measure partial-junction vertex offsets from a complete-capsid landmark trajectory.

For every complete ring (pentamer of A, hexamer of B/C/D) the true vertex is the
landmark centroid. For each contiguous 3-member subset and each 2-member subset
the script records where the true centre lies relative to the subset, in a
local frame built from the subset alone (plus the ring normal for pairs).
The output JSON is consumed by extract_patch_observables.py --calibration.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage


def kabsch(P, Q):
    """Rotation R, translation t with R @ p + t ~ q (least squares)."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    return R, qc - R @ pc


def pair_order(p, q, chain_p, chain_q, X, partner):
    """Canonical (p, q): lower chain letter first; for equal letters, by the cross-distance asymmetry to the dimer partners."""
    if chain_p != chain_q:
        return (p, q) if chain_p < chain_q else (q, p)
    # cross distances to the other landmark's dimer partner are strongly asymmetric around a ring
    # (about 24 A in the T4 capsid); order so that d(p, partner(q)) < d(q, partner(p))
    cross = np.linalg.norm(X[p] - X[partner[q]]) - np.linalg.norm(X[q] - X[partner[p]])
    return (p, q) if cross <= 0 else (q, p)


def ring_order(P):
    c = P.mean(0)
    n = np.linalg.svd(P - c)[2][2]
    e1 = P[0] - c
    e1 = e1 - n * np.dot(e1, n)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    ang = [math.atan2(np.dot(p - c, e2), np.dot(p - c, e1)) for p in P]
    return list(np.argsort(ang)), n


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--psf", required=True)
    ap.add_argument("--dcd", required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--stride", type=int, default=5)
    ap.add_argument("--cluster-cutoff", type=float, default=35.0)
    args = ap.parse_args()
    import MDAnalysis as mda
    u = mda.Universe(args.psf, args.dcd)
    chain = [a.segid[0] for a in u.atoms]
    segid = [a.segid for a in u.atoms]
    lookup = {sg: i for i, sg in enumerate(segid)}
    partner = {i: lookup[{"A": "B", "B": "A", "C": "D", "D": "C"}[sg[0]] + sg[1:]] for i, sg in enumerate(segid)}
    X0 = np.mean([u.atoms.positions.copy() for ts in u.trajectory[:20]], axis=0)
    lab = fcluster(linkage(X0, "single"), t=args.cluster_cutoff, criterion="distance")
    rings = {k: [i for i in range(len(chain)) if lab[i] == k] for k in set(lab)}
    capsid_centre = X0.mean(0)
    triples, pairs, radii, pairs_pf, templates = {}, {}, {}, {}, {}
    nframes = 0
    for ts in u.trajectory[::args.stride]:
        nframes += 1
        X = u.atoms.positions.astype(float)
        for mem in rings.values():
            P = X[mem]
            c = P.mean(0)
            order, n = ring_order(P)
            if np.dot(n, c - capsid_centre) < 0:
                n = -n  # outward
            m = len(mem)
            for i in mem:
                radii.setdefault(f"{m}:{chain[i]}", []).append(float(np.linalg.norm(X[i] - c)))
            for j in range(m):
                a, b, d = (mem[order[(j - 1) % m]], mem[order[j]], mem[order[(j + 1) % m]])
                Pt = X[[a, b, d]]
                c3 = Pt.mean(0)
                u1 = c3 - Pt[1]
                dist = np.linalg.norm(u1)
                u1 /= dist
                nn = np.cross(Pt[0] - Pt[1], Pt[2] - Pt[1])
                nn /= np.linalg.norm(nn)
                if np.dot(nn, n) < 0:
                    nn = -nn
                u2 = np.cross(nn, u1)
                v = c - c3
                triples.setdefault(f"{m}:{chain[b]}", []).append([np.dot(v, u1) / dist, np.dot(v, u2), np.dot(v, nn)])
                for step in (1, 2):
                    p, q = mem[order[j]], mem[order[(j + step) % m]]
                    mid = 0.5 * (X[p] + X[q])
                    chord = X[q] - X[p]
                    dvec = np.cross(n, chord)
                    dvec /= np.linalg.norm(dvec)
                    v = c - mid
                    if np.dot(v, dvec) < 0:
                        dvec = -dvec
                    key = f"{m}:{''.join(sorted(chain[p] + chain[q]))}:{step}"
                    pairs.setdefault(key, []).append([np.dot(v, dvec), np.dot(v, n), float(np.linalg.norm(chord))])
                    # partner frame: e1 along chord (lower chain letter -> higher), e2 toward the dimer partners, e3 = e1 x e2
                    pp, qq = pair_order(p, q, chain[p], chain[q], X, partner)
                    e1 = X[qq] - X[pp]
                    e1 /= np.linalg.norm(e1)
                    w = 0.5 * (X[partner[pp]] + X[partner[qq]]) - mid
                    e2 = w - np.dot(w, e1) * e1
                    e2 /= np.linalg.norm(e2)
                    e3 = np.cross(e1, e2)
                    pairs_pf.setdefault(key, []).append([np.dot(v, e1), np.dot(v, e2), np.dot(v, e3)])
                    # template: 4 present points (p, q, partner p, partner q) + centre, expressed in the first instance's frame
                    pts = np.array([X[pp], X[qq], X[partner[pp]], X[partner[qq]]])
                    tpl = templates.setdefault(key, {"ref": pts.copy(), "centres": [], "pts": []})
                    R, tt = kabsch(pts, tpl["ref"])
                    tpl["centres"].append(R @ c + tt)
                    tpl["pts"].append((R @ pts.T).T + tt)
    out = {"schema": "heva_vertex_offset_calibration_v1", "source_psf": args.psf, "source_dcd": args.dcd,
           "frames": nframes, "ring_count": {"5": sum(len(m) == 5 for m in rings.values()), "6": sum(len(m) == 6 for m in rings.values())},
           "radius": {k: {"mean": float(np.mean(v)), "sd": float(np.std(v))} for k, v in sorted(radii.items())},
           "triple": {k: {"f": float(np.mean(np.array(v)[:, 0])), "f_sd": float(np.std(np.array(v)[:, 0])),
                          "perp": float(np.mean(np.array(v)[:, 1])), "normal": float(np.mean(np.array(v)[:, 2])), "n": len(v)}
                      for k, v in sorted(triples.items())},
           "pair": {k: {"inplane": float(np.mean(np.array(v)[:, 0])), "inplane_sd": float(np.std(np.array(v)[:, 0])),
                        "normal": float(np.mean(np.array(v)[:, 1])), "normal_sd": float(np.std(np.array(v)[:, 1])),
                        "chord": float(np.mean(np.array(v)[:, 2])), "n": len(v)} for k, v in sorted(pairs.items())},
           "pair_partner_frame": {k: {"a": float(np.mean(np.array(v)[:, 0])), "b": float(np.mean(np.array(v)[:, 1])), "c": float(np.mean(np.array(v)[:, 2])),
                                      "sd": [float(x) for x in np.std(np.array(v), axis=0)], "n": len(v)} for k, v in sorted(pairs_pf.items())},
           "pair_template": {k: {"points": np.mean(np.array(v["pts"]), axis=0).tolist(), "centre": np.mean(np.array(v["centres"]), axis=0).tolist(),
                                 "centre_sd": float(np.mean(np.std(np.array(v["centres"]), axis=0))), "n": len(v["centres"]),
                                 "order": "p, q (chain(p) <= chain(q)), partner(p), partner(q)"} for k, v in sorted(templates.items())},
           "frames_convention": {"pair_template": "Kabsch-align template points onto the patch's 4 atoms, apply the same transform to template centre",
                                 "pair_partner_frame": "centre = mid + a*e1 + b*e2 + c*e3; e1 = unit(q - p) with chain(p) <= chain(q); e2 = unit component of (mean(partner(p), partner(q)) - mid) perpendicular to e1; e3 = e1 x e2","triple": "centre = c3 + f*(c3 - p_middle) + perp*u2 + normal*n_out; u2 = n_out x u1",
                                 "pair": "centre = midpoint + inplane*d + normal*n_out; d = unit(n_out x chord) oriented toward centre (away from patch neighbours)"}}
    args.output.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({k: out[k] for k in ("ring_count", "radius", "triple")}, indent=1))
    print(json.dumps({k: {"centre_sd": v["centre_sd"], "n": v["n"]} for k, v in out["pair_template"].items()}, indent=1))


if __name__ == "__main__":
    main()
