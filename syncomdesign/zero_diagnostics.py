from __future__ import annotations

from pathlib import Path
from typing import Mapping

from .community import build_community_model, classify_community_reactions
from .combinations import enumerate_all
from .fluxes import read_metabolite_aliases
from .io import resolve_path, write_tsv
from .medium import apply_community_medium, read_medium_file
from .objectives import add_all_species_active_constraint, solve_objective, total_biomass_flux
from .solvers import configure_solver
from .validation import load_models


def write_zero_biomass_diagnostics(cfg: Mapping[str, object], outdir: str | Path) -> Path:
    base_dir = Path(str(cfg["_base_dir"]))
    debug_dir = Path(outdir) / "debug_zero_fix"
    debug_dir.mkdir(parents=True, exist_ok=True)

    model_infos, validation_rows = load_models(cfg, base_dir)
    medium = read_medium_file(Path(resolve_path(base_dir, cfg["medium"]["file"])))
    aliases = read_metabolite_aliases(Path(resolve_path(base_dir, cfg["models"]["metabolite_aliases_file"])))
    names = [str(info["name"]) for info in model_infos]
    combos = enumerate_all(names)
    by_name = {str(info["name"]): info for info in model_infos}

    raw_rows = []
    singleton_rows = []
    shared_rows = []
    objective_rows = []
    biomass_rows = []
    constraint_rows = []
    external_stoich_rows = []
    interface_stoich_rows = []
    mass_balance_rows = []
    nutrient_rows = []

    for info in model_infos:
        raw_model = info["model"].copy()
        configure_solver(raw_model, cfg["solver"].get("name"), cfg["solver"].get("tolerance"), cfg["solver"].get("threads"))
        raw_model.objective = raw_model.reactions.get_by_id(str(info["biomass_rxn"]))
        raw_sol = raw_model.optimize()
        raw_rows.append({"strain_id": info["name"], "status": raw_sol.status, "biomass": raw_sol.objective_value})

    selected = combos[: min(5, len(combos))]
    for combo in selected:
        combo_id = "+".join(combo)
        community = build_community_model([by_name[strain] for strain in combo], cfg["community"])
        configure_solver(community, cfg["solver"].get("name"), cfg["solver"].get("tolerance"), cfg["solver"].get("threads"))
        pre_solution, _pre_trace = solve_objective(community, cfg["objective"], aliases, combo_id)
        singleton_rows.append({"combination_id": combo_id, "stage": "before_medium", "status": pre_solution.status, "total_biomass": total_biomass_flux(community, pre_solution)})

        classes = classify_community_reactions(community)
        medium_options = dict(cfg["medium"])
        medium_options["shared_environment_compartment"] = cfg["community"].get("shared_environment_compartment", "u")
        medium_out = apply_community_medium(community, medium, classes, medium_options, combination_id=combo_id)
        constraints = []
        if cfg["community"].get("require_all_species_active", False):
            constraints = add_all_species_active_constraint(community, float(cfg["community"].get("minimum_biomass_flux") or 1e-6))
        solution, objective_trace = solve_objective(community, cfg["objective"], aliases, combo_id)
        shared_rows.append({"combination_id": combo_id, "stage": "after_medium", "status": solution.status, "total_biomass": total_biomass_flux(community, solution)})
        objective_rows.append(objective_trace.as_row())
        constraint_rows.append({"combination_id": combo_id, "expected_constraints": len(combo), "actual_constraints": len(constraints), "constraint_names": ";".join(constraints)})

        for row in community.syncomdesign.get("biomassMap", []):
            rxn = community.reactions.get_by_id(row["biomass_rxn"])
            biomass_rows.append({"combination_id": combo_id, "strain_id": row["strain"], "biomass_rxn": rxn.id, "lower_bound": rxn.lower_bound, "upper_bound": rxn.upper_bound})
        for row in medium_out["external_medium_bounds"]:
            rxn = community.reactions.get_by_id(row["reaction_id"])
            external_stoich_rows.append(_stoich_row(combo_id, rxn, row["reaction_class"]))
        class_map = {row["reaction_id"]: row["classification"] for row in classes}
        for rxn_id, reaction_class in class_map.items():
            if reaction_class == "strain_shared_interface":
                interface_stoich_rows.append(_stoich_row(combo_id, community.reactions.get_by_id(rxn_id), reaction_class))

        for met in community.metabolites:
            if str(met.id).endswith("[u]"):
                linked = [rxn.id for rxn in met.reactions]
                mass_balance_rows.append({"combination_id": combo_id, "shared_metabolite": met.id, "linked_reactions": ";".join(sorted(linked)), "linked_count": len(linked)})
        for target in ["no3", "no2", "no", "n2o", "n2", "glu__L", "h2o", "h"]:
            nutrient_rows.append(_nutrient_path_row(combo_id, community, target))

    write_tsv(debug_dir / "01_raw_single_model_growth.tsv", raw_rows)
    write_tsv(debug_dir / "02_prefixed_singleton_growth.tsv", singleton_rows)
    write_tsv(debug_dir / "03_shared_singleton_growth.tsv", shared_rows)
    write_tsv(debug_dir / "04_multistrain_growth.tsv", shared_rows)
    write_tsv(debug_dir / "05_objective_integrity.tsv", objective_rows)
    write_tsv(debug_dir / "06_biomass_reaction_bounds.tsv", biomass_rows)
    write_tsv(debug_dir / "07_growth_constraints.tsv", constraint_rows)
    write_tsv(debug_dir / "08_external_exchange_stoichiometry.tsv", external_stoich_rows)
    write_tsv(debug_dir / "09_interface_stoichiometry.tsv", interface_stoich_rows)
    write_tsv(debug_dir / "10_shared_metabolite_mass_balance.tsv", mass_balance_rows)
    write_tsv(debug_dir / "11_nutrient_path_trace.tsv", nutrient_rows)
    write_tsv(debug_dir / "12_matlab_reference_diff.tsv", validation_rows)
    return debug_dir


def _stoich_row(combination_id: str, reaction: object, reaction_class: str) -> dict[str, object]:
    return {
        "combination_id": combination_id,
        "reaction_id": reaction.id,
        "reaction_class": reaction_class,
        "metabolites": ";".join(f"{met.id}:{coef:g}" for met, coef in reaction.metabolites.items()),
        "stoich_count": len(reaction.metabolites),
        "lower_bound": reaction.lower_bound,
        "upper_bound": reaction.upper_bound,
    }


def _nutrient_path_row(combination_id: str, model: object, target: str) -> dict[str, object]:
    external = f"R_EX_{target}_u"
    interfaces = [rxn.id for rxn in model.reactions if rxn.id.endswith(f"__R_EX_{target}_e")]
    shared = f"{target}[u]"
    return {
        "combination_id": combination_id,
        "target": target,
        "shared_metabolite": shared,
        "external_exchange_exists": external in {rxn.id for rxn in model.reactions},
        "interface_reactions": ";".join(sorted(interfaces)),
    }
