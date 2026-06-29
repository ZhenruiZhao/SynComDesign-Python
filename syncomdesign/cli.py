from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from .combinations import enumerate_all, write_all_combinations
from .community import build_community_model, build_community_trace, classify_community_reactions
from .config import apply_overrides, load_config
from .diagnostics import Diagnostics
from .fluxes import extract_fluxes, read_metabolite_aliases, write_flux_outputs
from .io import read_tsv, resolve_path, write_tsv
from .matlab_reference import compare_matlab_reference
from .medium import apply_community_medium, read_medium_file, write_medium_outputs
from .objectives import add_all_species_active_constraint, solve_objective
from .reporting import result_row, write_run_outputs
from .solvers import configure_solver
from .validation import detect_model_files, load_models, validate_project_inputs
from .zero_diagnostics import write_zero_biomass_diagnostics


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "run":
            return run_command(args)
        if args.command == "validate":
            return validate_command(args)
        if args.command == "compare-matlab":
            return compare_matlab_command(args)
        if args.command == "diagnose-zero":
            return diagnose_zero_command(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    parser.print_help()
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="syncomdesign", description="MATLAB-aligned SynComDesign Python CLI")
    sub = parser.add_subparsers(dest="command")

    run = sub.add_parser("run", help="Run SynComDesign")
    add_common_config_args(run)
    run.add_argument("--models", help="Override models directory")
    run.add_argument("--medium", help="Override medium file")
    run.add_argument("--outdir", help="Override output directory")
    run.add_argument("--solver", help="Override solver name")
    run.add_argument("--threads", type=int, default=None)
    run.add_argument("--tmpdir", default=None)

    validate = sub.add_parser("validate", help="Validate config and input paths")
    add_common_config_args(validate)

    compare = sub.add_parser("compare-matlab", help="Compare Python outputs to MATLAB reference exports")
    add_common_config_args(compare)
    compare.add_argument("--reference", required=True, help="MATLAB python_reference_exports directory")
    compare.add_argument("--outdir", default="results_compare_matlab", help="Alignment report output directory")
    compare.add_argument("--python-outdir", default=None, help="Python run output directory; defaults to config project.output_dir")
    compare.add_argument("--tolerance", type=float, default=1e-6)

    diagnose = sub.add_parser("diagnose-zero", help="Write zero-biomass diagnostic tables")
    add_common_config_args(diagnose)
    diagnose.add_argument("--models", help="Override models directory")
    diagnose.add_argument("--medium", help="Override medium file")
    diagnose.add_argument("--outdir", default="results", help="Diagnostic output directory")
    diagnose.add_argument("--solver", help="Override solver name")
    return parser


def add_common_config_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--config", required=True, help="Path to syncomdesign_config.yml")


def run_command(args: argparse.Namespace) -> int:
    cfg = apply_overrides(
        load_config(args.config),
        models=args.models,
        medium=args.medium,
        outdir=args.outdir,
        solver=args.solver,
        threads=args.threads,
        tmpdir=args.tmpdir,
    )
    base_dir = Path(cfg["_base_dir"])
    outdir = _resolve_output_dir(base_dir, cfg["project"]["output_dir"])
    diag = Diagnostics(outdir)
    diag.log("Starting SynComDesign Python run.")

    model_dir = Path(resolve_path(base_dir, cfg["models"]["directory"]))
    detected = detect_model_files(model_dir, str(cfg["models"].get("file_pattern", "*.xml")))
    write_tsv(outdir / "models_detected.tsv", detected, ["strain", "model_path"])
    model_infos, validation_rows = load_models(cfg, base_dir)
    strain_names = [str(info["name"]) for info in model_infos]
    combos = enumerate_all(
        strain_names,
        min_size=int(cfg["combinations"].get("min_size") or 1),
        max_size=cfg["combinations"].get("max_size"),
        target_strain=cfg["objective"].get("target_strain"),
        objective_mode=cfg["objective"].get("scenario_id"),
    )
    combo_rows = write_all_combinations(outdir / "all_combinations.tsv", combos)

    medium = read_medium_file(Path(resolve_path(base_dir, cfg["medium"]["file"])))
    aliases = read_metabolite_aliases(Path(resolve_path(base_dir, cfg["models"]["metabolite_aliases_file"])))
    tables: dict[str, list[dict[str, object]]] = {
        "model_validation": validation_rows,
        "all_combinations": combo_rows,
        "community_summary": [],
        "objective_trace": [],
        "reaction_classification": [],
        "community_build_trace": [],
    }
    medium_tables = {
        "medium_to_shared_exchange_mapping": [],
        "medium_mapping_warnings": [],
        "external_medium_bounds": [],
        "interface_bounds": [],
        "internal_transport_bounds": [],
    }
    flux_mapping_rows: list[dict[str, object]] = []
    flux_value_rows: list[dict[str, object]] = []

    by_name = {str(info["name"]): info for info in model_infos}
    for combo in combos:
        combo_id = "+".join(combo)
        started = time.perf_counter()
        try:
            diag.log(f"Starting combination {combo_id}.")
            combo_infos = [by_name[strain] for strain in combo]
            community = build_community_model(combo_infos, cfg["community"])
            diag.log(f"Built community {combo_id}.")
            configure_solver(community, cfg["solver"].get("name"), cfg["solver"].get("tolerance"), cfg["solver"].get("threads"))
            classes = classify_community_reactions(community)
            for row in classes:
                row["combination_id"] = combo_id
            tables["reaction_classification"].extend(classes)
            tables["community_build_trace"].extend(build_community_trace(combo_id, community))

            medium_options = dict(cfg["medium"])
            medium_options["shared_environment_compartment"] = cfg["community"].get("shared_environment_compartment", "u")
            medium_out = apply_community_medium(community, medium, classes, medium_options, combination_id=combo_id)
            for key in medium_tables:
                medium_tables[key].extend(medium_out[key])
            if cfg["community"].get("require_all_species_active", False):
                add_all_species_active_constraint(community, float(cfg["community"].get("minimum_biomass_flux") or 1e-6))
            solution, objective_trace = solve_objective(community, cfg["objective"], aliases, combo_id)
            diag.log(f"Solved combination {combo_id}: {objective_trace.status}.")
            tables["objective_trace"].append(objective_trace.as_row())
            flux_values, mapping_rows, value_rows = extract_fluxes(community, solution, aliases, combo_id)
            flux_mapping_rows.extend(mapping_rows)
            flux_value_rows.extend(value_rows)
            tables["community_summary"].append(result_row(combo_id, list(combo), community, solution, cfg["objective"]["type"], flux_values, started))
        except Exception as exc:
            diag.fail_combination(combo_id, str(exc), time.perf_counter() - started)
            continue

    write_run_outputs(outdir, tables)
    write_medium_outputs(outdir, medium_tables)
    write_flux_outputs(outdir, flux_mapping_rows, flux_value_rows)
    diag.write()
    diag.log("Finished SynComDesign Python run.")
    print(f"SynComDesign Python run completed: {outdir}")
    print(f"Combinations: {len(combos)}")
    print(f"Failed combinations: {len(diag.failed_combinations)}")
    return 0 if not diag.failed_combinations else 1


def validate_command(args: argparse.Namespace) -> int:
    cfg = load_config(args.config)
    base_dir = Path(cfg["_base_dir"])
    outdir = _resolve_output_dir(base_dir, cfg["project"]["output_dir"])
    checks = validate_project_inputs(cfg, base_dir)
    write_tsv(outdir / "validation_checks.tsv", checks, ["check", "path", "exists"])
    ok = all(bool(row["exists"]) for row in checks)
    print(f"Validation {'PASS' if ok else 'FAIL'}: {outdir / 'validation_checks.tsv'}")
    return 0 if ok else 1


def compare_matlab_command(args: argparse.Namespace) -> int:
    cfg = load_config(args.config)
    base_dir = Path(cfg["_base_dir"])
    python_outdir = Path(args.python_outdir) if args.python_outdir else _infer_python_outdir(base_dir, cfg, Path(args.outdir))
    report = compare_matlab_reference(python_outdir, args.reference, args.outdir, args.tolerance)
    print(f"MATLAB alignment: {report['status']}")
    print(f"Report: {Path(args.outdir) / 'matlab_alignment_report.md'}")
    return 0 if report["status"] == "PASS" else 1


def diagnose_zero_command(args: argparse.Namespace) -> int:
    cfg = apply_overrides(
        load_config(args.config),
        models=args.models,
        medium=args.medium,
        outdir=args.outdir,
        solver=args.solver,
    )
    base_dir = Path(cfg["_base_dir"])
    outdir = _resolve_output_dir(base_dir, cfg["project"]["output_dir"])
    debug_dir = write_zero_biomass_diagnostics(cfg, outdir)
    print(f"Zero-biomass diagnostics written: {debug_dir}")
    return 0


def _resolve_output_dir(base_dir: Path, value: object) -> Path:
    path = Path(str(value))
    return path if path.is_absolute() else base_dir / path


def _infer_python_outdir(base_dir: Path, cfg: dict[str, object], report_outdir: Path) -> Path:
    if report_outdir.name.endswith("_compare"):
        candidate = report_outdir.with_name(report_outdir.name[: -len("_compare")])
        if candidate.exists():
            return candidate
    return _resolve_output_dir(base_dir, cfg["project"]["output_dir"])


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
