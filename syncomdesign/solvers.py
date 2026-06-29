from __future__ import annotations


def configure_solver(model: object, solver_name: str | None = None, tolerance: float | None = None, threads: int | None = None) -> None:
    if not hasattr(model, "solver"):
        return
    if solver_name:
        try:
            model.solver = solver_name
        except Exception:
            pass
    if tolerance is not None:
        try:
            model.solver.configuration.tolerances.feasibility = tolerance
            model.solver.configuration.tolerances.optimality = tolerance
        except Exception:
            pass
    if threads is not None:
        try:
            model.solver.configuration.threads = threads
        except Exception:
            pass

