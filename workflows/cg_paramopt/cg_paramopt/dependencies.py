from __future__ import annotations


try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_IMPORT_ERROR = None
except ModuleNotFoundError as exc:
    plt = None
    MATPLOTLIB_IMPORT_ERROR = exc

try:
    import numpy as np
    NUMPY_IMPORT_ERROR = None
except ModuleNotFoundError as exc:
    np = None
    NUMPY_IMPORT_ERROR = exc

try:
    import pandas as pd
    PANDAS_IMPORT_ERROR = None
except ModuleNotFoundError as exc:
    pd = None
    PANDAS_IMPORT_ERROR = exc

try:
    from scipy.special import kl_div
    SCIPY_IMPORT_ERROR = None
except ModuleNotFoundError as exc:
    kl_div = None
    SCIPY_IMPORT_ERROR = exc

try:
    from hyperopt import Trials, fmin, hp, tpe
    HYPEROPT_IMPORT_ERROR = None
except Exception as exc:
    Trials = None
    fmin = None
    hp = None
    tpe = None
    HYPEROPT_IMPORT_ERROR = exc


def _format_missing(name: str, exc: Exception | None) -> str:
    if exc is None:
        return name
    return f"{name} ({exc})"


def require_runtime_dependencies(*, require_search: bool) -> None:
    missing = []
    if plt is None:
        missing.append(_format_missing("matplotlib", MATPLOTLIB_IMPORT_ERROR))
    if np is None:
        missing.append(_format_missing("numpy", NUMPY_IMPORT_ERROR))
    if pd is None:
        missing.append(_format_missing("pandas", PANDAS_IMPORT_ERROR))
    if kl_div is None:
        missing.append(_format_missing("scipy", SCIPY_IMPORT_ERROR))
    if require_search and (Trials is None or fmin is None or hp is None or tpe is None):
        missing.append(_format_missing("hyperopt", HYPEROPT_IMPORT_ERROR))

    if missing:
        raise RuntimeError("Missing runtime dependencies for optimization: " + ", ".join(missing))
