from __future__ import annotations


def require_runtime_dependencies() -> None:
    try:
        import scipy  # noqa: F401
    except ImportError as exc:
        raise ImportError("scipy is required for the assembly_opt workflow") from exc

