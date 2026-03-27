#!/usr/bin/env python3
"""Optimize CG parameters against AA reference data using HEVA."""

from __future__ import annotations

import argparse
import logging
import sys
import time
from functools import partial
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
WORKFLOW_DIR = SCRIPT_DIR.parent
if str(WORKFLOW_DIR) not in sys.path:
    sys.path.insert(0, str(WORKFLOW_DIR))

from cg_paramopt.config import DEFAULT_PARAM_BOUNDS, build_config
from cg_paramopt.dependencies import Trials, fmin, hp, require_runtime_dependencies, tpe
from cg_paramopt.metrics import combine_params, compute_histogram, load_aa_data
from cg_paramopt.runner import initialize_run_folder, run_single_simulation
from cg_paramopt.verify import run_static_verification


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Optimize CG parameters against AA reference data using HEVA.")
    parser.add_argument("--skip-verify", action="store_true", help="Skip the static verification pass.")
    parser.add_argument("--verify-only", action="store_true", help="Run only the static verification pass.")
    parser.add_argument("--max-evals", type=int, default=10, help="Number of Hyperopt evaluations.")
    parser.add_argument("--frames", type=int, default=5000, help="Number of sampled frames per evaluation.")
    parser.add_argument("--relax-sweeps", type=int, default=50000, help="Relaxation sweeps before sampling begins.")
    parser.add_argument("--sample-every", type=int, default=1, help="Sweeps between sampled frames.")
    parser.add_argument("--runs-dir", default=None, help="Output directory for optimization runs.")
    parser.add_argument("--base-config", default=None, help="Override the workflow base config.")
    parser.add_argument("--assemble", default=None, help="Override the HEVA assemble binary path.")
    parser.add_argument("--prepare-initial-frame", default=None, help="Override the HEVA initial-frame preparation binary path.")
    parser.add_argument("--initial-structure", default=None, help="Path to an AA-translated Initial_frame.dat input or a HEVA restart file.")
    parser.add_argument("--aa-data", default=None, help="Path to AAdata.csv.")
    parser.add_argument("--input-format", choices=("initial_frame", "restart"), default="initial_frame", help="Interpret --initial-structure as an AA-translated initial frame or a HEVA restart file.")
    parser.add_argument("--profile", choices=("test", "run", "extended"), default="test", help="HEVA engine profile.")
    parser.add_argument("--index-capacity", type=int, default=None, help="Required only for profile=extended.")
    parser.add_argument("--dry-run", action="store_true", help="Generate trial configs without executing HEVA.")
    parser.add_argument(
        "--single-run",
        nargs=3,
        type=float,
        metavar=("LOG10_KAPPA_L", "LOG10_KAPPA_THETA", "LOG10_KAPPA_PHI"),
        help="Run one fixed parameter set instead of Hyperopt.",
    )
    return parser.parse_args()


def log_verification(config) -> None:
    logging.info("Running static verification for HEVA CG-ParamOpt")
    for item in run_static_verification(config):
        log = getattr(logging, item.level, logging.info)
        log("Verification %s: %s", item.label, item.detail)


def hyperopt_objective(space, df_hist_aa, run_dir, geometry_params, config):
    x = [space["kappa_l"], space["kappa_theta"], space["kappa_phi"]]
    result = run_single_simulation(x, df_hist_aa, run_dir, geometry_params, config, dry_run=False)
    logging.info("Trial completed with error: %.3f", result.error)
    return result.error


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s", handlers=[logging.StreamHandler()])

    config = build_config(
        WORKFLOW_DIR,
        runs_dir=args.runs_dir,
        base_config=args.base_config,
        assemble=args.assemble,
        prepare_initial_frame=args.prepare_initial_frame,
        initial_structure=args.initial_structure,
        aa_data=args.aa_data,
        max_evals=args.max_evals,
        frames=args.frames,
        relax_sweeps=args.relax_sweeps,
        sample_every=args.sample_every,
        input_format=args.input_format,
        profile=args.profile,
        index_capacity=args.index_capacity,
    )

    if not args.skip_verify:
        log_verification(config)
        if args.verify_only:
            return 0

    require_runtime_dependencies(require_search=not args.single_run and not args.dry_run)

    df_aa = load_aa_data(config.resolve_aa_data_path())
    geometry_params = list(df_aa.mean())
    df_aa_combined = combine_params(df_aa)
    df_hist_aa = compute_histogram(df_aa_combined, ["auto"] * len(df_aa_combined))

    timestamp = time.strftime("%Y%m%d-%H%M%S")
    run_dir = config.resolve_runs_dir() / f"{timestamp}_hyperopt"
    seed = int(timestamp.split("-")[1])
    initialize_run_folder(run_dir, seed)

    if args.single_run:
        result = run_single_simulation(list(args.single_run), df_hist_aa, run_dir, geometry_params, config, dry_run=args.dry_run)
        logging.info("Single run command: %s", result.command)
        logging.info("Single run trial dir: %s", result.trial_dir)
        logging.info("Single run error: %.3f", result.error)
        return 0

    if args.dry_run:
        midpoint = [(low + high) / 2 for low, high in DEFAULT_PARAM_BOUNDS.values()]
        result = run_single_simulation(midpoint, df_hist_aa, run_dir, geometry_params, config, dry_run=True)
        logging.info("Dry run command: %s", result.command)
        logging.info("Dry run trial dir: %s", result.trial_dir)
        return 0

    space = {
        "kappa_l": hp.uniform("kappa_l", *DEFAULT_PARAM_BOUNDS["kappa_l"]),
        "kappa_theta": hp.uniform("kappa_theta", *DEFAULT_PARAM_BOUNDS["kappa_theta"]),
        "kappa_phi": hp.uniform("kappa_phi", *DEFAULT_PARAM_BOUNDS["kappa_phi"]),
    }
    objective = partial(hyperopt_objective, df_hist_aa=df_hist_aa, run_dir=run_dir, geometry_params=geometry_params, config=config)

    trials = Trials()
    best = fmin(fn=objective, space=space, algo=tpe.suggest, max_evals=config.max_evals, trials=trials)
    logging.info("Optimization complete")
    logging.info("Best parameters (log10 scale):")
    logging.info("  kappa_l: %.3f -> %.0f", best["kappa_l"], 10 ** best["kappa_l"])
    logging.info("  kappa_theta: %.3f -> %.0f", best["kappa_theta"], 10 ** best["kappa_theta"])
    logging.info("  kappa_phi: %.3f -> %.0f", best["kappa_phi"], 10 ** best["kappa_phi"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
