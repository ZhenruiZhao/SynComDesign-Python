from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping

from .fluxes import aliases_for_target, find_metabolite_exchanges


@dataclass
class ObjectiveTrace:
    combination_id: str
    scenario_id: int
    primary_objective: str
    primary_value: float
    secondary_objective: str = ""
    secondary_value: float = math.nan
    growth_fraction: float = math.nan
    target_strain: str = ""
    status: str = "unknown"

    def as_row(self) -> dict[str, object]:
        return self.__dict__.copy()


def scenario_type(scenario_id: int | str | None, explicit: str | None = None) -> str:
    if explicit:
        return explicit
    return {
        1: "total_biomass",
        2: "target_strain_biomass",
        3: "equal_composition",
        4: "fixed_composition",
        5: "growth_then_n2o_consumption",
    }[int(scenario_id or 1)]


def configure_objective(model: object, objective_config: Mapping[str, object]) -> str:
    mode = scenario_type(objective_config.get("scenario_id"), str(objective_config.get("type", "")) or None)
    biomass_rows = _biomass_map(model)
    if mode in {"total_biomass", "growth_then_n2o_consumption", "weighted_biomass_function", "fixed_composition", "equal_composition"}:
        _set_linear_objective(model, [row["biomass_rxn"] for row in biomass_rows])
    elif mode == "target_strain_biomass":
        target = str(objective_config.get("target_strain") or "")
        rows = [row for row in biomass_rows if str(row["strain"]).lower() == target.lower()]
        if not rows:
            raise ValueError(f'target_strain "{target}" is not present in this community')
        _set_linear_objective(model, [row["biomass_rxn"] for row in rows])
    else:
        raise ValueError(f"unknown objective mode: {mode}")
    return mode


def add_fixed_composition_constraint(model: object, ratios: list[float] | None = None) -> None:
    biomass_rxns = [row["biomass_rxn"] for row in _biomass_map(model) if row.get("biomass_rxn")]
    n = len(biomass_rxns)
    if n <= 1:
        return
    if ratios is None or len(ratios) == 0:
        ratios = [1.0] * n
    if len(ratios) != n:
        raise ValueError(f"expected {n} composition ratios, got {len(ratios)}")
    if any(value <= 0 for value in ratios):
        raise ValueError("composition ratios must be positive")
    total = sum(ratios)
    normalized = [value / total for value in ratios]
    try:
        for idx in range(1, n):
            first = _get_reaction(model, biomass_rxns[0])
            current = _get_reaction(model, biomass_rxns[idx])
            constraint = model.problem.Constraint(normalized[0] * current.flux_expression - normalized[idx] * first.flux_expression, lb=0, ub=0)
            model.add_cons_vars([constraint])
        model.syncomdesign["compositionRatios"] = normalized
    except AttributeError:
        meta = getattr(model, "syncomdesign", {})
        if isinstance(meta, dict):
            meta["compositionRatios"] = normalized


def add_all_species_active_constraint(model: object, minimum_biomass_flux: float) -> list[str]:
    constraint_names: list[str] = []
    for row in _biomass_map(model):
        reaction = _get_reaction(model, str(row.get("biomass_rxn", "")))
        if reaction is None:
            continue
        try:
            name = f"syncomdesign_min_biomass_{_safe_constraint_name(reaction.id)}"
            if hasattr(model.solver, "constraints") and name in model.solver.constraints:
                model.remove_cons_vars([model.solver.constraints[name]])
            constraint = model.problem.Constraint(reaction.flux_expression, lb=float(minimum_biomass_flux), name=name)
            model.add_cons_vars([constraint])
            constraint_names.append(name)
        except AttributeError:
            continue
    _update_solver(model)
    return constraint_names


def solve_objective(
    model: object,
    objective_config: Mapping[str, object],
    alias_table: Mapping[str, list[Mapping[str, str]]] | None = None,
    combination_id: str = "",
):
    mode = configure_objective(model, objective_config)
    scenario_id = int(objective_config.get("scenario_id") or 1)
    target_strain = str(objective_config.get("target_strain") or "")
    if mode in {"equal_composition", "fixed_composition"}:
        ratios = objective_config.get("composition_ratio")
        add_fixed_composition_constraint(model, list(ratios) if ratios else None)
    if mode == "growth_then_n2o_consumption":
        return _solve_growth_then_n2o(model, objective_config, alias_table, combination_id, scenario_id, target_strain)
    solution = _optimize(model)
    primary_value = _solution_objective_value(solution)
    status = _solution_status(solution)
    return solution, ObjectiveTrace(
        combination_id=combination_id,
        scenario_id=scenario_id,
        primary_objective=mode,
        primary_value=primary_value,
        target_strain=target_strain,
        status=status,
    )


def total_biomass_flux(model: object, solution: object | Mapping[str, float]) -> float:
    total = 0.0
    found = False
    for row in _biomass_map(model):
        reaction_id = str(row.get("biomass_rxn", ""))
        if not reaction_id:
            continue
        value = _solution_flux(model, solution, reaction_id)
        if not math.isnan(value):
            total += value
            found = True
    return total if found else math.nan


def _solve_growth_then_n2o(model, objective_config, alias_table, combination_id, scenario_id, target_strain):
    configure_objective(model, {"type": "total_biomass"})
    growth_solution = _optimize(model)
    max_growth = total_biomass_flux(model, growth_solution)
    if math.isnan(max_growth):
        max_growth = _solution_objective_value(growth_solution)
    growth_fraction = float(objective_config.get("growth_fraction") or 0.9)
    aliases = aliases_for_target(alias_table, "n2o", "nitrous_oxide")
    n2o_rxns = find_metabolite_exchanges(model, aliases)
    if not n2o_rxns:
        return growth_solution, ObjectiveTrace(
            combination_id, scenario_id, "total_biomass", max_growth, "", math.nan, growth_fraction, target_strain, _solution_status(growth_solution)
        )
    try:
        biomass_expr = sum(_get_reaction(model, row["biomass_rxn"]).flux_expression for row in _biomass_map(model))
        constraint = model.problem.Constraint(biomass_expr, lb=growth_fraction * max_growth)
        model.add_cons_vars([constraint])
        _set_linear_objective(model, n2o_rxns, coefficient=-1.0)
        functional_solution = _optimize(model)
        secondary = sum(max(0.0, -_solution_flux(model, functional_solution, rxn)) for rxn in n2o_rxns)
        return functional_solution, ObjectiveTrace(
            combination_id,
            scenario_id,
            "total_biomass",
            max_growth,
            ";".join(n2o_rxns),
            secondary,
            growth_fraction,
            target_strain,
            _solution_status(functional_solution),
        )
    except AttributeError:
        return growth_solution, ObjectiveTrace(
            combination_id,
            scenario_id,
            "total_biomass",
            max_growth,
            ";".join(n2o_rxns),
            math.nan,
            growth_fraction,
            target_strain,
            _solution_status(growth_solution),
        )


def _set_linear_objective(model: object, reaction_ids: list[str], coefficient: float = 1.0) -> None:
    reactions = [_get_reaction(model, reaction_id) for reaction_id in reaction_ids]
    reactions = [reaction for reaction in reactions if reaction is not None]
    if hasattr(model, "objective"):
        try:
            expression = sum(float(coefficient) * reaction.flux_expression for reaction in reactions)
            model.objective = expression
            model.objective_direction = "max"
            _update_solver(model)
            return
        except Exception:
            pass
    setattr(model, "objective_reactions", [(reaction.id, coefficient) for reaction in reactions])


def _optimize(model: object):
    if hasattr(model, "optimize"):
        return model.optimize()
    return {}


def _update_solver(model: object) -> None:
    try:
        model.solver.update()
    except Exception:
        pass


def _safe_constraint_name(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch == "_" else "_" for ch in str(value))


def _biomass_map(model: object) -> list[dict[str, str]]:
    meta = getattr(model, "syncomdesign", {})
    if isinstance(meta, dict):
        return list(meta.get("biomassMap", []))
    return []


def _get_reaction(model: object, reaction_id: str):
    reactions = getattr(model, "reactions", [])
    if hasattr(reactions, "get_by_id"):
        try:
            return reactions.get_by_id(reaction_id)
        except KeyError:
            return None
    for reaction in reactions:
        if getattr(reaction, "id", None) == reaction_id:
            return reaction
    return None


def _solution_flux(model: object, solution: object | Mapping[str, float], reaction_id: str) -> float:
    if isinstance(solution, Mapping):
        return float(solution.get(reaction_id, math.nan))
    if hasattr(solution, "fluxes"):
        try:
            return float(solution.fluxes[reaction_id])
        except Exception:
            return math.nan
    return math.nan


def _solution_status(solution: object) -> str:
    if isinstance(solution, Mapping):
        return str(solution.get("status", "unknown"))
    return str(getattr(solution, "status", getattr(solution, "stat", "unknown")))


def _solution_objective_value(solution: object) -> float:
    if isinstance(solution, Mapping):
        return float(solution.get("objective_value", math.nan))
    return float(getattr(solution, "objective_value", getattr(solution, "f", math.nan)))
