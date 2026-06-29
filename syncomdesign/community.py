from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Mapping

from .io import write_tsv


@dataclass
class SimpleMetabolite:
    id: str


@dataclass
class SimpleReaction:
    id: str
    lower_bound: float = -1000.0
    upper_bound: float = 1000.0
    metabolites: dict[str, float] = field(default_factory=dict)


@dataclass
class SimpleModel:
    reactions: list[SimpleReaction] = field(default_factory=list)
    metabolites: list[SimpleMetabolite] = field(default_factory=list)
    syncomdesign: dict[str, object] = field(default_factory=dict)


def matlab_make_valid_name(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_]", "_", str(value))
    if value == "":
        return "x"
    if not re.match(r"[A-Za-z]", value[0]):
        value = "x" + value
    return value


def build_community_model(model_infos: Iterable[Mapping[str, object]], config: Mapping[str, object] | None = None):
    config = dict(config or {})
    shared_compartment = str(config.get("shared_environment_compartment", "u"))
    infos = list(model_infos)
    try:
        import cobra
    except Exception as exc:  # pragma: no cover - depends on optional COBRApy install
        raise RuntimeError("COBRApy is required to build community models from SBML") from exc

    community = cobra.Model("SynComDesign_community")
    reaction_map: list[dict[str, object]] = []
    metabolite_map: list[dict[str, object]] = []
    biomass_map: list[dict[str, object]] = []
    shared_metabolites: dict[str, object] = {}
    community_reactions = []

    for info in infos:
        strain_id = str(info["name"])
        safe = matlab_make_valid_name(strain_id)
        source_model = info["model"]
        biomass_rxn = str(info.get("biomass_rxn") or info.get("biomassRxn") or "")
        exchange = boundary_exchange_reactions(source_model)
        exchange_rxns = {rxn_id for rxn_id, _met_id in exchange}
        exchange_mets = {met_id for _rxn_id, met_id in exchange}

        met_lookup = {}
        for met in source_model.metabolites:
            new_met = cobra.Metabolite(f"{safe}__{met.id}", name=getattr(met, "name", ""), compartment=getattr(met, "compartment", ""))
            met_lookup[met.id] = new_met
            role = "strain_exchange_metabolite" if met.id in exchange_mets else "strain_internal"
            metabolite_map.append({"community_met": new_met.id, "strain": strain_id, "source_met": met.id, "role": role})

        for rxn in source_model.reactions:
            new_rxn = cobra.Reaction(f"{safe}__{rxn.id}")
            new_rxn.lower_bound = float(rxn.lower_bound)
            new_rxn.upper_bound = float(rxn.upper_bound)
            new_rxn.name = getattr(rxn, "name", "")
            metabolites = {met_lookup[met.id]: coeff for met, coeff in rxn.metabolites.items()}
            if rxn.id in exchange_rxns:
                source_met = dict(exchange)[rxn.id]
                source_coeff = float(rxn.metabolites.get(source_model.metabolites.get_by_id(source_met), 0.0))
                shared_coeff = 1.0 if source_coeff == 0 else -source_coeff
                shared_id = f"{canonical_met_id(source_met)}[{shared_compartment}]"
                shared = shared_metabolites.get(shared_id)
                if shared is None:
                    shared = cobra.Metabolite(shared_id, compartment=shared_compartment)
                    shared_metabolites[shared_id] = shared
                metabolites[shared] = metabolites.get(shared, 0.0) + shared_coeff
            new_rxn.add_metabolites(metabolites)
            community_reactions.append(new_rxn)
            role = "strain_shared_interface" if rxn.id in exchange_rxns else "strain_internal"
            reaction_map.append({"community_rxn": new_rxn.id, "strain": strain_id, "source_rxn": rxn.id, "role": role})

        if biomass_rxn:
            biomass_map.append({"strain": strain_id, "biomass_rxn": f"{safe}__{biomass_rxn}"})

    community.add_reactions(community_reactions)
    existing_rxn_ids = {rxn.id for rxn in community.reactions}
    external_map: dict[str, str] = {}
    external_reactions = []
    for shared_id in shared_metabolites.keys():
        shared = shared_metabolites[shared_id]
        rxn_id = shared_exchange_rxn_id(shared_id, shared_compartment)
        if rxn_id not in existing_rxn_ids:
            rxn = cobra.Reaction(rxn_id)
            rxn.lower_bound = -1000.0
            rxn.upper_bound = 1000.0
            rxn.add_metabolites({shared: -1.0})
            external_reactions.append(rxn)
            existing_rxn_ids.add(rxn_id)
        external_map[shared_id] = rxn_id
        reaction_map.append({"community_rxn": rxn_id, "strain": "external", "source_rxn": rxn_id, "role": "external_medium_exchange"})
    if external_reactions:
        community.add_reactions(external_reactions)

    community.syncomdesign = {
        "strainNames": [str(info["name"]) for info in infos],
        "reactionMap": reaction_map,
        "metaboliteMap": metabolite_map,
        "biomassMap": biomass_map,
        "externalExchangeMap": external_map,
        "externalSharedExchangeRxns": list(external_map.values()),
        "sharedCompartment": shared_compartment,
    }
    return community


def classify_community_reactions(model: object) -> list[dict[str, object]]:
    role_map = _role_map(model)
    rows = []
    for reaction in _iter_reactions(model):
        rxn_id = reaction.id
        role = role_map.get(rxn_id, "")
        if role == "external_medium_exchange" or _is_shared_external_exchange(reaction):
            cls = "external_medium_exchange"
            medium_applies = True
        elif _is_transport_reaction(rxn_id):
            cls = "internal_transport"
            medium_applies = False
        elif role == "strain_shared_interface" or "__R_EX_" in rxn_id or "__EX_" in rxn_id:
            cls = "strain_shared_interface"
            medium_applies = False
        elif role in {"strain_internal", "exchange"} or _stoich_count(reaction) > 1:
            cls = "metabolic_reaction"
            medium_applies = False
        else:
            cls = "unknown"
            medium_applies = False
        rows.append({"reaction_id": rxn_id, "classification": cls, "medium_applies": medium_applies})
    return rows


def build_community_trace(combination_id: str, model: object) -> list[dict[str, object]]:
    class_map = {row["reaction_id"]: row["classification"] for row in classify_community_reactions(model)}
    meta = _meta(model)
    reaction_map = meta.get("reactionMap", [])
    external_map = meta.get("externalExchangeMap", {})
    reverse_external = {rxn: met for met, rxn in external_map.items()}
    rows = []
    for row in reaction_map:
        rxn_id = str(row["community_rxn"])
        reaction = _get_reaction(model, rxn_id)
        rows.append(
            {
                "combination_id": combination_id,
                "strain_id": row.get("strain", ""),
                "original_reaction": row.get("source_rxn", ""),
                "new_reaction": rxn_id,
                "reaction_class": class_map.get(rxn_id, "unknown"),
                "original_metabolite": "",
                "shared_metabolite": reverse_external.get(rxn_id, ""),
                "lb": getattr(reaction, "lower_bound", ""),
                "ub": getattr(reaction, "upper_bound", ""),
                "note": _trace_note(str(row.get("role", ""))),
            }
        )
    return rows


def write_community_debug_tables(outdir: str | Path, combination_id: str, model: object) -> None:
    outdir = Path(outdir)
    write_tsv(outdir / "reaction_classification.tsv", classify_community_reactions(model))
    write_tsv(outdir / "community_build_trace.tsv", build_community_trace(combination_id, model))


def boundary_exchange_reactions(model: object) -> list[tuple[str, str]]:
    out = []
    for reaction in _iter_reactions(model):
        rxn_id = reaction.id
        if not (rxn_id.startswith("EX_") or rxn_id.startswith("R_EX_")):
            continue
        if _is_transport_reaction(rxn_id):
            continue
        met_id = _first_metabolite_id(reaction)
        if met_id:
            out.append((rxn_id, met_id))
        else:
            out.append((rxn_id, exchange_met_from_rxn_id(rxn_id)))
    return out


def canonical_met_id(met_id: str) -> str:
    canonical = re.sub(r"\[[^\]]+\]$", "", str(met_id))
    canonical = re.sub(r"_[A-Za-z][A-Za-z0-9]*$", "", canonical)
    canonical = re.sub(r"^M_", "", canonical)
    return canonical


def exchange_met_from_rxn_id(rxn_id: str) -> str:
    met = re.sub(r"^R_EX_", "", str(rxn_id))
    met = re.sub(r"^EX_", "", met)
    met = re.sub(r"_(e|u)$", "_e", met)
    return met


def shared_exchange_rxn_id(shared_met: str, shared_compartment: str = "u") -> str:
    base = re.sub(rf"\[{re.escape(shared_compartment)}\]$", "", str(shared_met))
    base = re.sub(r"[^A-Za-z0-9_]", "_", base)
    return f"R_EX_{base}_{shared_compartment}"


def _trace_note(role: str) -> str:
    if role == "external_medium_exchange":
        return "external shared exchange for medium"
    if role == "strain_shared_interface":
        return "strain exchange connected to shared pool"
    return ""


def _role_map(model: object) -> dict[str, str]:
    return {str(row["community_rxn"]): str(row["role"]) for row in _meta(model).get("reactionMap", [])}


def _meta(model: object) -> dict[str, object]:
    meta = getattr(model, "syncomdesign", {})
    return meta if isinstance(meta, dict) else {}


def _is_shared_external_exchange(reaction: object) -> bool:
    return str(reaction.id).startswith("R_EX_") and _stoich_count(reaction) == 1


def _is_transport_reaction(rxn_id: str) -> bool:
    return re.search(r"(tex|tpp|t2pp|t3pp|tipp|tppi)$", str(rxn_id)) is not None


def _stoich_count(reaction: object) -> int:
    metabolites = getattr(reaction, "metabolites", {})
    return sum(1 for coeff in metabolites.values() if abs(float(coeff)) > 0)


def _first_metabolite_id(reaction: object) -> str:
    metabolites = getattr(reaction, "metabolites", {})
    for met in metabolites:
        return getattr(met, "id", str(met))
    return ""


def _iter_reactions(model: object):
    return list(getattr(model, "reactions", []))


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
