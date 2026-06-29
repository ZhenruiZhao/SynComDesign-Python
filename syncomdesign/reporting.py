from __future__ import annotations

import math
import time
from pathlib import Path
from typing import Mapping

from .io import write_tsv
from .objectives import total_biomass_flux


COMMUNITY_SUMMARY_FIELDS = [
    "combination_id",
    "community_size",
    "strain_names",
    "feasible",
    "solver_status",
    "objective_mode",
    "total_biomass",
    "active_strains",
    "nitrate_uptake",
    "nitrite_uptake",
    "nitrite_secretion",
    "no_uptake",
    "no_secretion",
    "n2o_uptake",
    "n2o_secretion",
    "n2o_net_flux",
    "n2_secretion",
    "nitrate_uptake_per_biomass",
    "n2o_uptake_per_biomass",
    "n2_production_per_biomass",
    "minimum_growth_satisfied",
    "runtime_seconds",
    "warning_message",
]


def result_row(
    combination_id: str,
    strain_names: list[str],
    model: object,
    solution: object | Mapping[str, float],
    objective_mode: str,
    flux_values: Mapping[str, Mapping[str, float]],
    started_at: float,
) -> dict[str, object]:
    biomass_by_strain = _strain_biomass(model, solution)
    total_biomass = total_biomass_flux(model, solution)
    if math.isnan(total_biomass):
        total_biomass = _objective_value(solution)
    active = [strain for strain in strain_names if biomass_by_strain.get(strain, 0.0) > 1e-9]
    row = {
        "combination_id": combination_id,
        "community_size": len(strain_names),
        "strain_names": ";".join(strain_names),
        "feasible": _status(solution) in {"optimal", "1"},
        "solver_status": _status(solution),
        "objective_mode": objective_mode,
        "total_biomass": total_biomass,
        "active_strains": ";".join(active),
        "nitrate_uptake": flux_values["no3"]["uptake"],
        "nitrite_uptake": flux_values["no2"]["uptake"],
        "nitrite_secretion": flux_values["no2"]["secretion"],
        "no_uptake": flux_values["no"]["uptake"],
        "no_secretion": flux_values["no"]["secretion"],
        "n2o_uptake": flux_values["n2o"]["uptake"],
        "n2o_secretion": flux_values["n2o"]["secretion"],
        "n2o_net_flux": flux_values["n2o"]["net_flux"],
        "n2_secretion": flux_values["n2"]["secretion"],
        "nitrate_uptake_per_biomass": _safe_divide(flux_values["no3"]["uptake"], total_biomass),
        "n2o_uptake_per_biomass": _safe_divide(flux_values["n2o"]["uptake"], total_biomass),
        "n2_production_per_biomass": _safe_divide(flux_values["n2"]["secretion"], total_biomass),
        "minimum_growth_satisfied": all(biomass_by_strain.get(strain, 0.0) > 1e-9 for strain in strain_names),
        "runtime_seconds": time.perf_counter() - started_at,
        "warning_message": "",
    }
    for strain in strain_names:
        row[f"strain_biomass_{_matlab_make_valid_name(strain)}"] = biomass_by_strain.get(strain, 0.0)
    return row


def write_run_outputs(outdir: str | Path, tables: Mapping[str, list[dict[str, object]]]) -> None:
    outdir = Path(outdir)
    write_tsv(outdir / "community_summary.tsv", tables.get("community_summary", []))
    write_tsv(outdir / "model_validation.tsv", tables.get("model_validation", []), ["strain", "model_path", "valid", "biomass_rxn", "warning_message"])
    write_tsv(outdir / "all_combinations.tsv", tables.get("all_combinations", []))
    write_tsv(outdir / "objective_trace.tsv", tables.get("objective_trace", []))
    write_tsv(outdir / "reaction_classification.tsv", tables.get("reaction_classification", []))
    write_tsv(outdir / "community_build_trace.tsv", tables.get("community_build_trace", []))


def _strain_biomass(model: object, solution: object | Mapping[str, float]) -> dict[str, float]:
    out = {}
    meta = getattr(model, "syncomdesign", {})
    biomass_map = meta.get("biomassMap", []) if isinstance(meta, dict) else []
    for row in biomass_map:
        out[str(row["strain"])] = _flux(solution, str(row["biomass_rxn"]))
    return out


def _flux(solution: object | Mapping[str, float], reaction_id: str) -> float:
    if isinstance(solution, Mapping):
        return float(solution.get(reaction_id, 0.0))
    if hasattr(solution, "fluxes"):
        try:
            return float(solution.fluxes[reaction_id])
        except Exception:
            return 0.0
    return 0.0


def _status(solution: object | Mapping[str, float]) -> str:
    if isinstance(solution, Mapping):
        return str(solution.get("status", "unknown"))
    return str(getattr(solution, "status", getattr(solution, "stat", "unknown")))


def _objective_value(solution: object | Mapping[str, float]) -> float:
    if isinstance(solution, Mapping):
        return float(solution.get("objective_value", math.nan))
    return float(getattr(solution, "objective_value", getattr(solution, "f", math.nan)))


def _safe_divide(numerator: float, denominator: float) -> float:
    if denominator in {0, None} or math.isnan(float(denominator)):
        return math.nan
    return numerator / denominator


def _matlab_make_valid_name(value: str) -> str:
    import re

    value = re.sub(r"[^A-Za-z0-9_]", "_", str(value))
    if not value or not re.match(r"[A-Za-z]", value[0]):
        value = "x" + value
    return value

