#!/usr/bin/env python3
"""Extract HEVA-style per-frame observables from AA patch trajectories.

Input: a PSF/DCD pair with one landmark atom per monomer, segid = chain letter
(A, B, C, D) followed by the asymmetric-unit index. Dimers are A_i-B_i and
C_i-D_i. Landmarks are clustered into CG vertices; half-edge XY runs from the
vertex containing X to the vertex containing Y with HEVA types CD=0, BA=1,
AB=2, DC=3. Faces are oriented to HEVA's two T=4 face cycles, (AB, CD, BA)
and (DC, DC, DC). Observables follow src/geometry.cpp's long_v1 sampler:
lengths per dimer, dihedrals per internal dimer (angle between the two face
normals), and corner angles per in-face contact with the same class rules.

Vertex positions from incomplete clusters are corrected geometrically
(estimator "calibrated"): a contiguous 3-of-k ring subset has its centroid
shifted along the middle-member axis by the regular-polygon offset; a 2-of-k
subset uses the chord midpoint shifted in the local tangent plane by
R*cos(pi/k) away from the patch, with R supplied per ring size.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage

TYPE_OF = {("A", "B"): 2, ("B", "A"): 1, ("C", "D"): 0, ("D", "C"): 3}
TYPE_NAME = {0: "CD", 1: "BA", 2: "AB", 3: "DC"}
FACE_CYCLES = {(2, 0, 1), (3, 3, 3)}
PHI_NAMES = ["p33", "p12", "p01", "p20"]


def get_dihedral_type(e, n, p, oe, on, op):
    cd = (0, 3)
    if e in cd and n == 1 and p == 2 and oe in cd and on == 1 and op == 2:
        return 0
    if e == 1 and n == 2 and p in cd and oe == 2 and on in cd and op == 1:
        return 1
    if e == 2 and n in cd and p == 1 and oe == 1 and on == 2 and op in cd:
        return 1
    if e == 0 and n == 1 and p == 2 and oe == 3 and on == 3 and op == 3:
        return 0
    if e == 3 and n == 3 and p == 3 and oe == 0 and on == 1 and op == 2:
        return 0
    if e == 3 and n == 1 and p == 2 and oe == 0 and on == 0 and op == 0:
        return 0
    if e == 0 and n == 0 and p == 0 and oe == 3 and on == 1 and op == 2:
        return 0
    return 3


def get_angle_type(e, n):
    if (e, n) in {(0, 0), (3, 3), (3, 0), (0, 3)}:
        return 0
    if (e, n) == (1, 2):
        return 1
    if (e, n) in {(0, 1), (3, 1)}:
        return 2
    if (e, n) in {(2, 0), (2, 3)}:
        return 3
    return 0


RELAXED = False


def phi_class(e, n):
    table = {(3, 3): 0, (1, 2): 1, (0, 1): 2, (2, 0): 3}
    if RELAXED and e in (0, 3) and n in (0, 3):
        return 0
    return table.get((e, n), -1)


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


def circumcenter(P):
    a, b, c = P
    ab, ac = b - a, c - a
    n = np.cross(ab, ac)
    d = 2.0 * np.dot(n, n)
    return a + (np.dot(ac, ac) * np.cross(n, ab) + np.dot(ab, ab) * np.cross(ac, n)) / d


def unit(v):
    return v / np.linalg.norm(v)


class PatchTopology:
    def __init__(self, segids, positions, cluster_cutoff, relaxed=False):
        self.relaxed = relaxed
        self.segids = list(segids)
        self.chain = [s[0] for s in self.segids]
        self.index = [int(s[1:]) for s in self.segids]
        for s in self.segids:
            if s[0] not in "ABCD" or not s[1:].isdigit():
                raise ValueError(f"Unexpected segid {s}")
        self.n = len(self.segids)
        self.dimers = self._pair_dimers()
        self.partner = {}
        for a, b in self.dimers:
            self.partner[a] = b
            self.partner[b] = a
        self.member_of = self._cluster(positions, cluster_cutoff)
        self.vertices = sorted(set(self.member_of))
        self.members = {v: [i for i in range(self.n) if self.member_of[i] == v] for v in self.vertices}
        self.ring_size = {}
        for v, mem in self.members.items():
            letters = {self.chain[i] for i in mem}
            if letters == {"A"}:
                self.ring_size[v] = 5
            elif letters <= {"B", "C", "D"}:
                self.ring_size[v] = 6
            else:
                raise ValueError(f"Vertex {v} mixes fivefold and sixfold chains: {[self.segids[i] for i in mem]}")
            if len(mem) > self.ring_size[v]:
                raise ValueError(f"Vertex {v} has too many members")
        for a, b in self.dimers:
            if self.member_of[a] == self.member_of[b]:
                raise ValueError(f"Dimer {self.segids[a]}-{self.segids[b]} landed in one vertex; lower the cutoff")
        self.he = self._half_edges()
        self.faces = self._faces()

    def _pair_dimers(self):
        lookup = {(c, i): k for k, (c, i) in enumerate(zip(self.chain, self.index))}
        dimers = []
        for k, (c, i) in enumerate(zip(self.chain, self.index)):
            if c in "AC":
                partner = lookup.get(("B" if c == "A" else "D", i))
                if partner is None:
                    raise ValueError(f"No partner for {self.segids[k]}")
                dimers.append((k, partner))
        if 2 * len(dimers) != self.n:
            raise ValueError("Unpaired monomers present")
        return dimers

    def _cluster(self, positions, cutoff):
        Z = linkage(positions, method="single")
        labels = fcluster(Z, t=cutoff, criterion="distance")
        return [int(x) for x in labels]

    def _half_edges(self):
        he = []
        for a, b in self.dimers:
            fwd = {"id": len(he), "vin": self.member_of[a], "vout": self.member_of[b],
                   "type": TYPE_OF[(self.chain[a], self.chain[b])], "op": len(he) + 1,
                   "next": -1, "prev": -1, "label": f"{self.segids[a]}{self.segids[b]}"}
            bwd = {"id": len(he) + 1, "vin": self.member_of[b], "vout": self.member_of[a],
                   "type": TYPE_OF[(self.chain[b], self.chain[a])], "op": len(he),
                   "next": -1, "prev": -1, "label": f"{self.segids[b]}{self.segids[a]}"}
            he += [fwd, bwd]
        return he

    def _faces(self):
        """Orient triangles to HEVA face cycles; resolve ambiguity by manifold consistency."""
        by_pair = {}
        for h in self.he:
            by_pair.setdefault((h["vin"], h["vout"]), []).append(h["id"])
        for pair, ids in by_pair.items():
            if len(ids) > 1:
                raise ValueError(f"Two dimers between vertices {pair}; topology not supported")
        allowed = set(FACE_CYCLES)
        if self.relaxed:
            # any CD/DC-only cycle (non-T4 sheets); orientation then follows manifold consistency
            allowed |= {(a, b, c) for a in (0, 3) for b in (0, 3) for c in (0, 3)}
        V = self.vertices
        candidates = []
        for i in range(len(V)):
            for j in range(i + 1, len(V)):
                for k in range(j + 1, len(V)):
                    tri = (V[i], V[j], V[k])
                    opts = []
                    for cyc in [(tri[0], tri[1], tri[2]), (tri[0], tri[2], tri[1])]:
                        ids = [by_pair.get((cyc[m], cyc[(m + 1) % 3]), [None])[0] for m in range(3)]
                        if None in ids:
                            continue
                        types = tuple(self.he[h]["type"] for h in ids)
                        rots = {types[r:] + types[:r] for r in range(3)}
                        if rots & allowed:
                            opts.append((0 if rots & FACE_CYCLES else 1, ids))
                    if opts:
                        candidates.append(sorted(opts))
        # Greedy assignment: unambiguous triangles first, then choose the orientation whose
        # half-edges are unused (adjacent faces traverse a shared dimer in opposite directions).
        used = set()
        faces = []
        pending = sorted(candidates, key=len)
        while pending:
            progress = False
            rest = []
            for opts in pending:
                viable = [ids for _, ids in opts if not (set(ids) & used)]
                if not viable:
                    raise ValueError("Half-edge used by two faces; orientation conflict")
                forced = len(opts) == 1 or len(viable) == 1
                touches = any(set(ids) & {self.he[h]["op"] for h in used} for ids in viable)
                if forced or touches:
                    ids = viable[0]
                    faces.append(ids)
                    used |= set(ids)
                    progress = True
                else:
                    rest.append(opts)
            pending = rest
            if pending and not progress:
                # disconnected or fully ambiguous: seed with the HEVA-preferred orientation
                opts = pending.pop(0)
                ids = opts[0][1]
                faces.append(ids)
                used |= set(ids)
        for f in faces:
            for m, h in enumerate(f):
                self.he[h]["next"] = f[(m + 1) % 3]
                self.he[h]["prev"] = f[(m - 1) % 3]
        return faces

    def summary(self):
        return {
            "n_monomers": self.n, "n_dimers": len(self.dimers), "n_vertices": len(self.vertices),
            "vertices": {str(v): {"members": [self.segids[i] for i in self.members[v]],
                                  "ring_size": self.ring_size[v], "present": len(self.members[v])}
                         for v in self.vertices},
            "faces": [[f"{TYPE_NAME[self.he[h]['type']]}:{self.he[h]['label']}" for h in f] for f in self.faces],
            "boundary_half_edges": [f"{TYPE_NAME[h['type']]}:{h['label']}" for h in self.he if h["next"] == -1],
        }


class VertexEstimator:
    def __init__(self, topo, mode, radius, calibration=None):
        self.t = topo
        self.mode = mode
        self.radius = radius
        self.cal = calibration
        if mode == "capsid" and not calibration:
            raise ValueError("--estimator capsid requires --calibration")
        self.chain = topo.chain
        self.triple_terms = "f"

    def _middle(self, P):
        d = np.linalg.norm(P[:, None] - P[None], axis=-1)
        return int(np.argmin(d.sum(1)))

    def _outward(self, n, point, X):
        return n if np.dot(n, point - X.mean(0)) >= 0 else -n

    def estimate(self, X):
        t = self.t
        out = {}
        pending2 = []
        for v in t.vertices:
            mem = t.members[v]
            P = X[mem]
            k, m = t.ring_size[v], len(P)
            if m == k or self.mode == "centroid":
                out[v] = P.mean(0)
            elif m == 3:
                if self.mode == "circumcenter":
                    out[v] = circumcenter(P)
                    continue
                mid = self._middle(P)
                c3 = P.mean(0)
                if self.mode == "capsid":
                    rec = self.cal["triple"][f"{k}:{self.chain[mem[mid]]}"]
                    u1 = c3 - P[mid]
                    dist = np.linalg.norm(u1)
                    u1 /= dist
                    others = [i for i in range(3) if i != mid]
                    nn = unit(np.cross(P[others[0]] - P[mid], P[others[1]] - P[mid]))
                    nn = self._outward(nn, c3, X)
                    u2 = np.cross(nn, u1)
                    if self.triple_terms == "full":
                        out[v] = c3 + rec["f"] * dist * u1 + rec["perp"] * u2 + rec["normal"] * nn
                    else:
                        out[v] = c3 + rec["f"] * dist * u1
                else:
                    c = math.cos(2 * math.pi / k)
                    f = (1 + 2 * c) / (2 - 2 * c)
                    out[v] = c3 + f * (c3 - P[mid])
            elif m == 2:
                out[v] = P.mean(0)
                pending2.append(v)
            else:
                raise ValueError(f"Vertex {v}: {m} of {k} members is not supported")
        if pending2 and self.mode != "centroid":
            normals = self._vertex_normals(out)
            for v in pending2:
                mem = t.members[v]
                P = X[mem]
                k = t.ring_size[v]
                chord = P[1] - P[0]
                n = self._outward(normals[v], P.mean(0), X)
                d = np.cross(n, chord)
                if np.linalg.norm(d) < 1e-9:
                    raise ValueError(f"Vertex {v}: degenerate tangent plane")
                d = unit(d)
                nbr = [out[h["vout"]] for h in t.he if h["vin"] == v]
                away = P.mean(0) - np.mean(nbr, axis=0)
                if np.dot(d, away) < 0:
                    d = -d
                chord_len = np.linalg.norm(chord)
                if self.mode == "capsid":
                    step = 1 if chord_len < 22.5 else 2
                    key = f"{k}:{''.join(sorted(self.chain[mem[0]] + self.chain[mem[1]]))}:{step}"
                    rec = self.cal["pair_template"][key]
                    pp, qq = pair_order(mem[0], mem[1], self.chain[mem[0]], self.chain[mem[1]], X, t.partner)
                    pts = np.array([X[pp], X[qq], X[t.partner[pp]], X[t.partner[qq]]])
                    R, tt = kabsch(np.array(rec["points"]), pts)
                    out[v] = R @ np.array(rec["centre"]) + tt
                else:
                    R = self.radius[k]
                    step = 1 if chord_len < 1.5 * R * math.sin(math.pi / k) * 1.4 else 2
                    out[v] = P.mean(0) + R * math.cos(step * math.pi / k) * d
        return out

    def _vertex_normals(self, V):
        t = self.t
        acc = {v: np.zeros(3) for v in t.vertices}
        for f in t.faces:
            a, b, c = (V[t.he[h]["vin"]] for h in f)
            n = np.cross(b - a, c - a)
            for h in f:
                acc[t.he[h]["vin"]] += n
        for v in acc:
            if np.linalg.norm(acc[v]) < 1e-9:
                # vertex outside any face: use global patch normal
                acc[v] = sum(np.cross(V[t.he[f[1]]["vin"]] - V[t.he[f[0]]["vin"]],
                                      V[t.he[f[2]]["vin"]] - V[t.he[f[0]]["vin"]]) for f in t.faces)
            acc[v] = unit(acc[v])
        return acc


def occupancy_labels(topo):
    """Vertex -> 'P<present>/5' for pentamer junctions, 'H<present>/6' for hexamer junctions."""
    return {v: f"{'P' if topo.ring_size[v] == 5 else 'H'}{len(topo.members[v])}/{topo.ring_size[v]}" for v in topo.vertices}


def frame_rows(topo, V):
    he = topo.he
    weak = {v for v in topo.vertices if len(topo.members[v]) == 2}
    occ = occupancy_labels(topo)
    vec = {h["id"]: V[h["vout"]] - V[h["vin"]] for h in he}
    length = {i: float(np.linalg.norm(v)) for i, v in vec.items()}
    normal = {}
    for h in he:
        if h["next"] != -1:
            normal[h["id"]] = unit(np.cross(vec[h["id"]], vec[h["next"]]))
    rows = []
    for h in he:
        op = he[h["op"]]
        lc = 0 if h["type"] in (0, 3) else 1
        if h["id"] < op["id"]:
            rows.append((f"l{lc}", h["id"], h["type"], lc, h["type"], length[h["id"]], h["label"], length[h["id"]],
                         len({h["vin"], h["vout"]} & weak), f"{occ[h['vin']]}-{occ[h['vout']]}"))
            if h["next"] != -1 and h["prev"] != -1 and op["next"] != -1 and op["prev"] != -1:
                theta = math.acos(max(-1.0, min(1.0, float(np.dot(normal[h["id"]], normal[op["id"]])))))
                ec = get_dihedral_type(h["type"], he[h["next"]]["type"], he[h["prev"]]["type"],
                                       op["type"], he[op["next"]]["type"], he[op["prev"]]["type"])
                # HEVA orientation test (bend_energy): the opposite face's far vertex must lie on the
                # negative side of this face's normal; otherwise the fold is on the "wrong" side.
                far_here = V[he[h["next"]]["vout"]]
                far_op = V[he[op["next"]]["vout"]]
                wrong_side = float(np.dot(far_op - far_here, normal[h["id"]])) > 0
                rows.append((f"t{lc}", h["id"], h["type"], lc, ec, theta, h["label"], -theta if wrong_side else theta,
                             len({h["vin"], h["vout"], he[h["next"]]["vout"], he[op["next"]]["vout"]} & weak),
                             f"{occ[h['vin']]}-{occ[h['vout']]}|{occ[he[h['next']]['vout']]},{occ[he[op['next']]['vout']]}"))
        if h["next"] == -1 or h["prev"] == -1:
            continue
        pc = phi_class(h["type"], he[h["next"]]["type"])
        if pc < 0:
            continue
        nx = he[h["next"]]
        cosine = float(np.dot(vec[op["id"]], vec[nx["id"]]) / (length[op["id"]] * length[nx["id"]]))
        phi = math.acos(max(-1.0, min(1.0, cosine)))
        rows.append((PHI_NAMES[pc], h["id"], h["type"], pc, get_angle_type(h["type"], nx["type"]), phi,
                     f"{h['label']}>{nx['label']}", phi, len({h["vout"], h["vin"], nx["vout"]} & weak),
                     f"{occ[h['vout']]}|{occ[h['vin']]},{occ[nx['vout']]}"))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--psf", type=Path, required=True)
    ap.add_argument("--dcd", type=Path, required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--start-frame", type=int, default=0)
    ap.add_argument("--dt-ps", type=float, default=None, help="Frame spacing in ps, recorded in the summary")
    ap.add_argument("--estimator", choices=["capsid", "calibrated", "centroid", "circumcenter"], default="capsid",
                    help="capsid: offsets measured on a complete capsid (needs --calibration); calibrated: regular-polygon geometry")
    ap.add_argument("--calibration", type=Path, default=None, help="JSON from calibrate_vertex_offsets.py")
    ap.add_argument("--triple-terms", choices=["f", "full"], default="f", help="capsid mode: use only the axial factor f, or also the perpendicular/normal offsets")
    ap.add_argument("--cluster-cutoff", type=float, default=35.0, help="Single-linkage cutoff in input units")
    ap.add_argument("--ring-radius-5", type=float, default=13.9)
    ap.add_argument("--ring-radius-6", type=float, default=16.4)
    ap.add_argument("--cluster-frames", type=int, default=50, help="Frames averaged before clustering")
    ap.add_argument("--relaxed-faces", action="store_true", help="Also accept CD,CD,CD faces (non-T4 sheets) and report their corners as p33")
    args = ap.parse_args()

    import MDAnalysis as mda
    u = mda.Universe(str(args.psf), str(args.dcd))
    segids = [a.segid for a in u.atoms]
    ref = np.zeros((len(segids), 3))
    nref = 0
    for ts in u.trajectory[args.start_frame:args.start_frame + args.cluster_frames]:
        ref += u.atoms.positions
        nref += 1
    ref /= nref
    global RELAXED
    RELAXED = args.relaxed_faces
    topo = PatchTopology(segids, ref, args.cluster_cutoff, relaxed=args.relaxed_faces)
    calibration = json.loads(args.calibration.read_text()) if args.calibration else None
    est = VertexEstimator(topo, args.estimator, {5: args.ring_radius_5, 6: args.ring_radius_6}, calibration)
    est.triple_terms = args.triple_terms

    args.output_dir.mkdir(parents=True, exist_ok=True)
    out_csv = args.output_dir / f"{args.name}_observables.csv"
    acc = {}
    nframes = 0
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["frame", "source_frame", "observable", "contact_id", "half_edge_type",
                    "observable_class", "energy_class", "value", "contact_label", "signed_value", "weak_junctions", "end_tags"])
        for fi, ts in enumerate(u.trajectory[args.start_frame::args.stride]):
            V = est.estimate(u.atoms.positions.astype(float))
            for name, cid, ht, oc, ec, val, label, sval, tag, ends in frame_rows(topo, V):
                w.writerow([fi, ts.frame, name, cid, ht, oc, ec, f"{val:.6f}", label, f"{sval:.6f}", tag, ends])
                acc.setdefault(name, []).append(val)
                if name.startswith("t"):
                    acc.setdefault(name + "_signed", []).append(sval)
            nframes += 1
    stats = {k: {"mean": float(np.mean(v)), "sd": float(np.std(v)), "n": len(v),
                 **({"fraction_negative": float(np.mean(np.array(v) < 0))} if k.endswith("_signed") else {})}
             for k, v in sorted(acc.items())}
    summary = {"name": args.name, "psf": str(args.psf), "dcd": str(args.dcd), "frames_used": nframes,
               "stride": args.stride, "start_frame": args.start_frame, "dt_ps": args.dt_ps,
               "estimator": args.estimator, "calibration": str(args.calibration) if args.calibration else None, "triple_terms": args.triple_terms,
               "cluster_cutoff": args.cluster_cutoff, "relaxed_faces": args.relaxed_faces,
               "ring_radius": {"5": args.ring_radius_5, "6": args.ring_radius_6},
               "topology": topo.summary(), "raw_stats": stats,
               "notes": ["Lengths are in input coordinate units (not frame-normalized).",
                         "Angles in radians, unsigned, HEVA long_v1 definitions."]}
    (args.output_dir / f"{args.name}_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"name": args.name, "frames": nframes, "topology": topo.summary(), "raw_stats": stats}, indent=1))


if __name__ == "__main__":
    sys.exit(main())
