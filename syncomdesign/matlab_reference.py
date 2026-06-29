from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable

from .io import read_tsv, write_tsv

TABLES_TO_COMPARE = [
    ("models_detected", "reference_models_detected.tsv", "models_detected.tsv"),
    ("all_combinations", "reference_all_combinations.tsv", "all_combinations.tsv"),
    ("medium_mapping", "reference_medium_to_shared_exchange_mapping.tsv", "medium_to_shared_exchange_mapping.tsv"),
    ("reaction_classification", "reference_reaction_classification.tsv", "reaction_classification.tsv"),
    ("external_medium_bounds", "reference_external_medium_bounds.tsv", "external_medium_bounds.tsv"),
    ("interface_bounds", "reference_interface_bounds.tsv", "interface_bounds.tsv"),
    ("internal_transport_bounds", "reference_internal_transport_bounds.tsv", "internal_transport_bounds.tsv"),
    ("community_summary", "reference_community_summary.tsv", "community_summary.tsv"),
    ("flux_mapping", "reference_flux_mapping.tsv", "flux_mapping.tsv"),
    ("objective_trace", "reference_objective_trace.tsv", "objective_trace.tsv"),
]

NUMERIC_COLUMNS = {
    "lower_bound",
    "upper_bound",
    "before_lower_bound",
    "before_upper_bound",
    "lb",
    "ub",
    "primary_value",
    "secondary_value",
    "growth_fraction",
    "total_biomass",
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
    "runtime_seconds",
    "net_flux",
    "uptake",
    "secretion",
}

PRIORITY_TABLES = {"external_medium_bounds", "interface_bounds", "internal_transport_bounds", "medium_mapping"}


def compare_matlab_reference(
    python_outdir: str | Path,
    reference_dir: str | Path,
    outdir: str | Path,
    tolerance: float = 1e-6,
) -> dict[str, object]:
    python_outdir = Path(python_outdir)
    reference_dir = Path(reference_dir)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    diff_rows = []
    for table_name, ref_name, py_name in TABLES_TO_COMPARE:
        ref_path = reference_dir / ref_name
        py_path = python_outdir / py_name
        if not ref_path.exists() or not py_path.exists():
            status = "MISSING"
            summary_rows.append(
                {
                    "table": table_name,
                    "status": status,
                    "reference_rows": _row_count(ref_path),
                    "python_rows": _row_count(py_path),
                    "differences": 1,
                }
            )
            diff_rows.append({"table": table_name, "key": "", "column": "", "reference": str(ref_path.exists()), "python": str(py_path.exists()), "difference": "missing_file"})
            continue
        ref_rows = read_tsv(ref_path)
        py_rows = read_tsv(py_path)
        table_diffs = _compare_rows(table_name, ref_rows, py_rows, tolerance)
        summary_rows.append(
            {
                "table": table_name,
                "status": "PASS" if not table_diffs else "FAIL",
                "reference_rows": len(ref_rows),
                "python_rows": len(py_rows),
                "differences": len(table_diffs),
            }
        )
        diff_rows.extend(table_diffs)

    write_tsv(outdir / "matlab_alignment_summary.tsv", summary_rows, ["table", "status", "reference_rows", "python_rows", "differences"])
    write_tsv(outdir / "matlab_alignment_differences.tsv", diff_rows, ["table", "key", "column", "reference", "python", "difference"])
    status = _overall_status(summary_rows)
    _write_report(outdir / "matlab_alignment_report.md", status, summary_rows, diff_rows)
    return {"status": status, "summary": summary_rows, "differences": diff_rows}


def _compare_rows(table_name: str, ref_rows: list[dict[str, str]], py_rows: list[dict[str, str]], tolerance: float) -> list[dict[str, str]]:
    ref_sorted = sorted(ref_rows, key=_row_key)
    py_sorted = sorted(py_rows, key=_row_key)
    diffs = []
    if len(ref_sorted) != len(py_sorted):
        diffs.append(
            {
                "table": table_name,
                "key": "row_count",
                "column": "",
                "reference": str(len(ref_sorted)),
                "python": str(len(py_sorted)),
                "difference": "row_count_mismatch",
            }
        )
    for idx, (ref_row, py_row) in enumerate(zip(ref_sorted, py_sorted), start=1):
        columns = sorted(set(ref_row) | set(py_row))
        key = _display_key(ref_row, idx)
        for column in columns:
            ref_value = ref_row.get(column, "")
            py_value = py_row.get(column, "")
            if column in {"runtime_seconds"}:
                continue
            if column in NUMERIC_COLUMNS:
                if not _numeric_close(ref_value, py_value, tolerance):
                    diffs.append(_diff(table_name, key, column, ref_value, py_value, "numeric_mismatch"))
            elif _normalize_scalar(ref_value) != _normalize_scalar(py_value):
                diffs.append(_diff(table_name, key, column, ref_value, py_value, "value_mismatch"))
    return diffs


def _overall_status(summary_rows: Iterable[dict[str, object]]) -> str:
    rows = list(summary_rows)
    if all(row["status"] == "PASS" for row in rows):
        return "PASS"
    if any(row["table"] in PRIORITY_TABLES and row["status"] != "PASS" for row in rows):
        return "FAIL"
    return "PARTIAL"


def _write_report(path: Path, status: str, summary_rows: list[dict[str, object]], diff_rows: list[dict[str, str]]) -> None:
    lines = [
        "# MATLAB Alignment Report",
        "",
        f"Overall status: {status}",
        "",
        "## Summary",
        "",
        "| Table | Status | Reference rows | Python rows | Differences |",
        "| --- | --- | ---: | ---: | ---: |",
    ]
    for row in summary_rows:
        lines.append(f"| {row['table']} | {row['status']} | {row['reference_rows']} | {row['python_rows']} | {row['differences']} |")
    lines.extend(["", "## First Differences", ""])
    if not diff_rows:
        lines.append("No differences found within tolerance.")
    else:
        for row in diff_rows[:100]:
            lines.append(
                f"- {row['table']} `{row['key']}` column `{row['column']}`: reference `{row['reference']}` vs python `{row['python']}` ({row['difference']})"
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _row_count(path: Path) -> int:
    return len(read_tsv(path)) if path.exists() else 0


def _row_key(row: dict[str, str]) -> tuple[str, ...]:
    preferred = [
        "combination_id",
        "reaction_id",
        "medium_exchange_rxn",
        "medium_metabolite",
        "canonical_id",
        "strain",
        "model_path",
    ]
    values = [str(row.get(key, "")) for key in preferred if key in row]
    if values:
        return tuple(values)
    return tuple(f"{key}={row[key]}" for key in sorted(row))


def _display_key(row: dict[str, str], idx: int) -> str:
    keys = [key for key in ["combination_id", "reaction_id", "medium_exchange_rxn", "canonical_id", "strain"] if key in row]
    if not keys:
        return str(idx)
    return "|".join(str(row.get(key, "")) for key in keys)


def _numeric_close(left: str, right: str, tolerance: float) -> bool:
    try:
        a = float(left)
        b = float(right)
    except ValueError:
        return str(left) == str(right)
    if math.isnan(a) and math.isnan(b):
        return True
    return abs(a - b) <= tolerance


def _normalize_scalar(value: str) -> str:
    normalized = str(value).strip().lower()
    if normalized in {"true", "1"}:
        return "true"
    if normalized in {"false", "0"}:
        return "false"
    return str(value)


def _diff(table: str, key: str, column: str, ref: str, py: str, kind: str) -> dict[str, str]:
    return {"table": table, "key": key, "column": column, "reference": str(ref), "python": str(py), "difference": kind}
