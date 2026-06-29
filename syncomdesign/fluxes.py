from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Iterable, Mapping

from .io import read_tsv, write_tsv

TARGETS = [
    ("no3", "nitrate"),
    ("no2", "nitrite"),
    ("no", "nitric_oxide"),
    ("n2o", "nitrous_oxide"),
    ("n2", "dinitrogen"),
]


def read_metabolite_aliases(path: str | Path) -> dict[str, list[dict[str, str]]]:
    aliases: dict[str, list[dict[str, str]]] = {}
    for row in read_tsv(path):
        canonical = str(row.get("canonical_id", ""))
        aliases.setdefault(canonical, []).append(
            {
                "alias": str(row.get("alias", "")),
                "category": str(row.get("category", "")),
            }
        )
    return aliases


def aliases_for_target(alias_table: Mapping[str, list[Mapping[str, str]]] | None, canonical: str, category: str) -> list[str]:
    aliases = [canonical, f"{canonical}_e", f"EX_{canonical}_e", category]
    for row in (alias_table or {}).get(canonical, []):
        alias = str(row.get("alias", ""))
        row_category = str(row.get("category", ""))
        if alias and alias not in aliases:
            aliases.append(alias)
        if row_category == category and row_category not in aliases:
            aliases.append(row_category)
    return aliases


def find_metabolite_exchanges(model: object, aliases: Iterable[str]) -> list[str]:
    normalized_aliases = {_normalize_token(alias) for alias in aliases if str(alias) != ""}
    hits = []
    for reaction in _detect_exchange_reactions(model):
        rxn_tokens = {
            _normalize_token(reaction["reaction_id"]),
            _normalize_token(_strip_sbml_prefix(reaction["reaction_id"])),
            _normalize_token(reaction["metabolite_id"]),
            _normalize_token(_strip_sbml_prefix(reaction["metabolite_id"])),
            _normalize_token(_strip_compartment(reaction["metabolite_id"])),
            _normalize_token(_strip_sbml_prefix(_strip_compartment(reaction["metabolite_id"]))),
        }
        if normalized_aliases.intersection(rxn_tokens):
            hits.append(reaction["reaction_id"])
    return hits


def extract_fluxes(
    model: object,
    solution: object | Mapping[str, float],
    alias_table: Mapping[str, list[Mapping[str, str]]] | None = None,
    combination_id: str = "",
) -> tuple[dict[str, dict[str, float]], list[dict[str, object]], list[dict[str, object]]]:
    values: dict[str, dict[str, float]] = {}
    mapping_rows: list[dict[str, object]] = []
    value_rows: list[dict[str, object]] = []
    for canonical, category in TARGETS:
        aliases = aliases_for_target(alias_table, canonical, category)
        reaction_ids = find_metabolite_exchanges(model, aliases)
        flux_values = [_solution_flux(model, solution, reaction_id) for reaction_id in reaction_ids]
        net_flux = sum(value for value in flux_values if not math.isnan(value)) if flux_values else math.nan
        uptake = max(0.0, -net_flux) if not math.isnan(net_flux) else math.nan
        secretion = max(0.0, net_flux) if not math.isnan(net_flux) else math.nan
        values[canonical] = {"net_flux": net_flux, "uptake": uptake, "secretion": secretion}
        mapping_rows.append(
            {
                "combination_id": combination_id,
                "canonical_id": canonical,
                "category": category,
                "aliases": ";".join(aliases),
                "reaction_id": ";".join(reaction_ids),
            }
        )
        value_rows.append(
            {
                "combination_id": combination_id,
                "canonical_id": canonical,
                "category": category,
                "reaction_id": ";".join(reaction_ids),
                "net_flux": net_flux,
                "uptake": uptake,
                "secretion": secretion,
            }
        )
    return values, mapping_rows, value_rows


def write_flux_outputs(outdir: str | Path, mapping_rows: list[dict[str, object]], value_rows: list[dict[str, object]]) -> None:
    write_tsv(outdir / Path("flux_mapping.tsv"), mapping_rows, ["combination_id", "canonical_id", "category", "aliases", "reaction_id"])
    write_tsv(
        outdir / Path("flux_values.tsv"),
        value_rows,
        ["combination_id", "canonical_id", "category", "reaction_id", "net_flux", "uptake", "secretion"],
    )


def _detect_exchange_reactions(model: object) -> list[dict[str, str]]:
    out = []
    for reaction in getattr(model, "reactions", []):
        rxn_id = reaction.id
        metabolites = getattr(reaction, "metabolites", {})
        stoich_count = sum(1 for coeff in metabolites.values() if abs(float(coeff)) > 0)
        met_id = _first_metabolite_id(reaction)
        single_external = stoich_count == 1 and re.search(r"(\[[eu]\]$|_[eu]$)", met_id or "") is not None
        if rxn_id.startswith("EX_") or rxn_id.startswith("R_EX_") or single_external:
            out.append({"reaction_id": rxn_id, "metabolite_id": met_id})
    return out


def _solution_flux(model: object, solution: object | Mapping[str, float], reaction_id: str) -> float:
    if isinstance(solution, Mapping):
        return float(solution.get(reaction_id, math.nan))
    if hasattr(solution, "fluxes"):
        try:
            return float(solution.fluxes[reaction_id])
        except Exception:
            return math.nan
    if hasattr(solution, "x"):
        reactions = list(getattr(model, "reactions", []))
        for idx, reaction in enumerate(reactions):
            if reaction.id == reaction_id:
                return float(solution.x[idx])
    return math.nan


def _first_metabolite_id(reaction: object) -> str:
    for metabolite in getattr(reaction, "metabolites", {}):
        return getattr(metabolite, "id", str(metabolite))
    return ""


def _strip_sbml_prefix(value: str) -> str:
    return re.sub(r"^[RM]_", "", str(value))


def _strip_compartment(value: str) -> str:
    out = re.sub(r"\[[^\]]+\]$", "", str(value))
    return re.sub(r"_[A-Za-z][A-Za-z0-9]*$", "", out)


def _normalize_token(value: str) -> str:
    return str(value).lower()

