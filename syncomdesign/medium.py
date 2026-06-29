from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, MutableMapping

from .io import as_float, read_tsv, write_tsv


@dataclass(frozen=True)
class MediumEntry:
    metabolite: str
    exchange_rxn: str
    lower_bound: float
    upper_bound: float
    description: str = ""


class MediumTable(list[MediumEntry]):
    pass


def read_medium_file(path: str | Path) -> MediumTable:
    rows = read_tsv(path)
    out = MediumTable()
    for row in rows:
        out.append(
            MediumEntry(
                metabolite=str(row.get("metabolite", "")),
                exchange_rxn=str(row.get("exchange_rxn", "")),
                lower_bound=as_float(row.get("lower_bound"), 0.0),
                upper_bound=as_float(row.get("upper_bound"), 1000.0),
                description=str(row.get("description", "")),
            )
        )
    return out


def complete_external_medium(
    community_model: object,
    medium_table: Iterable[MediumEntry],
    external_exchange_map: Mapping[str, str] | Iterable[Mapping[str, str]],
) -> list[dict[str, object]]:
    external = _normalise_external_map(external_exchange_map)
    listed = {
        map_medium_to_shared_exchange(entry, community_model, aliases={}).get("shared_external_exchange")
        for entry in medium_table
    }
    rows = []
    for shared_metabolite, reaction_id in external.items():
        rows.append(
            {
                "shared_metabolite": shared_metabolite,
                "shared_external_exchange": reaction_id,
                "listed_in_medium": reaction_id in listed,
            }
        )
    return rows


def apply_community_medium(
    community_model: object,
    medium_table: Iterable[MediumEntry],
    reaction_classes: Mapping[str, str] | Iterable[Mapping[str, object]],
    options: Mapping[str, object] | None = None,
    outdir: str | Path | None = None,
    combination_id: str | None = None,
) -> dict[str, list[dict[str, object]]]:
    options = dict(options or {})
    condition = str(options.get("condition", "anaerobic")).lower()
    close_unlisted = bool(options.get("close_unlisted_external_medium_uptakes", True))
    shared_compartment = str(options.get("shared_environment_compartment", "u"))
    aliases = options.get("aliases", {}) or {}
    class_map = _normalise_class_map(reaction_classes)

    before = _bounds_by_reaction(community_model)
    external_rxns = [rid for rid, cls in class_map.items() if cls == "external_medium_exchange"]
    if close_unlisted:
        for reaction_id in external_rxns:
            reaction = _get_reaction(community_model, reaction_id)
            if reaction is not None and _get_lb(reaction) < 0:
                _set_lb(reaction, 0.0)

    mapping_rows: list[dict[str, object]] = []
    warning_rows: list[dict[str, object]] = []
    applied_rxns: set[str] = set()
    for entry in medium_table:
        mapped = map_medium_to_shared_exchange(
            entry,
            community_model,
            aliases=aliases,
            shared_compartment=shared_compartment,
            external_rxns=external_rxns,
        )
        lb = entry.lower_bound
        ub = entry.upper_bound
        if _is_oxygen_medium(entry.exchange_rxn, entry.metabolite):
            if condition == "anaerobic":
                lb = 0.0
            elif condition == "microaerobic":
                lb = max(lb, -1.0)
        mapped.update({"lower_bound": lb, "upper_bound": ub, "combination_id": combination_id or ""})
        if mapped["found"]:
            reaction_id = str(mapped["shared_external_exchange"])
            reaction = _get_reaction(community_model, reaction_id)
            if reaction is not None:
                _set_lb(reaction, lb)
                _set_ub(reaction, ub)
                applied_rxns.add(reaction_id)
        else:
            warning = {
                "combination_id": combination_id or "",
                "medium_exchange_rxn": entry.exchange_rxn,
                "medium_metabolite": entry.metabolite,
                "expected_shared_metabolite": mapped["shared_metabolite"],
                "warning": mapped["warning"],
            }
            warning_rows.append(warning)
        mapping_rows.append(mapped)

    external_bounds = _bounds_for_class(community_model, class_map, "external_medium_exchange", before, combination_id)
    interface_bounds = _bounds_for_class(community_model, class_map, "strain_shared_interface", before, combination_id)
    internal_bounds = _bounds_for_class(community_model, class_map, "internal_transport", before, combination_id)

    changed_forbidden = [
        row
        for row in interface_bounds + internal_bounds
        if row["lower_bound"] != row["before_lower_bound"] or row["upper_bound"] != row["before_upper_bound"]
    ]
    if changed_forbidden:
        raise AssertionError("community medium modified strain-interface or internal transport reactions")

    outputs = {
        "medium_to_shared_exchange_mapping": mapping_rows,
        "medium_mapping_warnings": warning_rows,
        "external_medium_bounds": external_bounds,
        "interface_bounds": interface_bounds,
        "internal_transport_bounds": internal_bounds,
    }
    if outdir is not None:
        write_medium_outputs(outdir, outputs)
    return outputs


def map_medium_to_shared_exchange(
    medium_entry: MediumEntry | Mapping[str, object],
    community_model: object,
    aliases: Mapping[str, str] | None = None,
    shared_compartment: str = "u",
    external_rxns: Iterable[str] | None = None,
) -> dict[str, object]:
    entry = _coerce_entry(medium_entry)
    base = _canonical_medium_base(entry, aliases or {})
    shared_metabolite = f"{base}[{shared_compartment}]"
    external_map = _external_map_from_model(community_model)
    reaction_id = external_map.get(shared_metabolite, "")
    if external_rxns is not None and reaction_id not in set(external_rxns):
        reaction_id = ""
    found = reaction_id != ""
    warning = "" if found else (
        f"No external shared exchange found for {entry.exchange_rxn}/{entry.metabolite}; "
        f"expected shared metabolite {shared_metabolite}."
    )
    return {
        "medium_exchange_rxn": entry.exchange_rxn,
        "medium_metabolite": entry.metabolite,
        "shared_metabolite": shared_metabolite,
        "shared_external_exchange": reaction_id,
        "found": found,
        "warning": warning,
    }


def write_medium_outputs(outdir: str | Path, outputs: Mapping[str, list[dict[str, object]]]) -> None:
    outdir = Path(outdir)
    write_tsv(
        outdir / "medium_to_shared_exchange_mapping.tsv",
        outputs.get("medium_to_shared_exchange_mapping", []),
        [
            "combination_id",
            "medium_exchange_rxn",
            "medium_metabolite",
            "shared_metabolite",
            "shared_external_exchange",
            "found",
            "lower_bound",
            "upper_bound",
            "warning",
        ],
    )
    write_tsv(
        outdir / "medium_mapping_warnings.tsv",
        outputs.get("medium_mapping_warnings", []),
        ["combination_id", "medium_exchange_rxn", "medium_metabolite", "expected_shared_metabolite", "warning"],
    )
    for key, filename in [
        ("external_medium_bounds", "external_medium_bounds.tsv"),
        ("interface_bounds", "interface_bounds.tsv"),
        ("internal_transport_bounds", "internal_transport_bounds.tsv"),
    ]:
        write_tsv(
            outdir / filename,
            outputs.get(key, []),
            ["combination_id", "reaction_id", "lower_bound", "upper_bound", "reaction_class"],
        )


def _coerce_entry(entry: MediumEntry | Mapping[str, object]) -> MediumEntry:
    if isinstance(entry, MediumEntry):
        return entry
    return MediumEntry(
        metabolite=str(entry.get("metabolite", "")),
        exchange_rxn=str(entry.get("exchange_rxn", "")),
        lower_bound=as_float(entry.get("lower_bound"), 0.0),
        upper_bound=as_float(entry.get("upper_bound"), 1000.0),
        description=str(entry.get("description", "")),
    )


def _canonical_medium_base(entry: MediumEntry, aliases: Mapping[str, str]) -> str:
    base = entry.metabolite
    if base == "" or base.lower() == "missing":
        base = entry.exchange_rxn
        if base.startswith("R_EX_"):
            base = base[5:]
        elif base.startswith("EX_"):
            base = base[3:]
    if base.startswith("M_"):
        base = base[2:]
    for suffix in (f"_{comp}" for comp in ["e", "u"]):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
    if "[" in base and base.endswith("]"):
        base = base.rsplit("[", 1)[0]
    return str(aliases.get(base, base))


def _is_oxygen_medium(exchange_rxn: str, metabolite: str) -> bool:
    rxn = exchange_rxn[2:] if exchange_rxn.startswith("R_") else exchange_rxn
    met = metabolite
    for suffix in ["_e", "_u"]:
        if met.endswith(suffix):
            met = met[: -len(suffix)]
    return rxn.lower() in {"ex_o2_e", "ex_o2_u", "ex_o2", "o2"} or met.lower() in {"ex_o2_e", "ex_o2_u", "ex_o2", "o2"}


def _normalise_external_map(external_exchange_map: Mapping[str, str] | Iterable[Mapping[str, str]]) -> dict[str, str]:
    if isinstance(external_exchange_map, Mapping):
        return {str(k): str(v) for k, v in external_exchange_map.items()}
    out = {}
    for row in external_exchange_map:
        out[str(row["shared_metabolite"])] = str(row["external_exchange_rxn"])
    return out


def _external_map_from_model(model: object) -> dict[str, str]:
    meta = getattr(model, "syncomdesign", None)
    if isinstance(meta, MutableMapping):
        raw = meta.get("externalExchangeMap", {})
    else:
        raw = getattr(meta, "externalExchangeMap", {}) if meta is not None else {}
    return _normalise_external_map(raw)


def _normalise_class_map(reaction_classes: Mapping[str, str] | Iterable[Mapping[str, object]]) -> dict[str, str]:
    if isinstance(reaction_classes, Mapping):
        return {str(k): str(v) for k, v in reaction_classes.items()}
    out = {}
    for row in reaction_classes:
        out[str(row["reaction_id"])] = str(row["classification"])
    return out


def _bounds_by_reaction(model: object) -> dict[str, tuple[float, float]]:
    return {reaction.id: (_get_lb(reaction), _get_ub(reaction)) for reaction in _iter_reactions(model)}


def _bounds_for_class(
    model: object,
    class_map: Mapping[str, str],
    wanted: str,
    before: Mapping[str, tuple[float, float]],
    combination_id: str | None,
) -> list[dict[str, object]]:
    rows = []
    for reaction_id, reaction_class in class_map.items():
        if reaction_class != wanted:
            continue
        reaction = _get_reaction(model, reaction_id)
        if reaction is None:
            continue
        before_lb, before_ub = before.get(reaction_id, (_get_lb(reaction), _get_ub(reaction)))
        rows.append(
            {
                "combination_id": combination_id or "",
                "reaction_id": reaction_id,
                "reaction_class": reaction_class,
                "before_lower_bound": before_lb,
                "before_upper_bound": before_ub,
                "lower_bound": _get_lb(reaction),
                "upper_bound": _get_ub(reaction),
            }
        )
    return rows


def _iter_reactions(model: object):
    reactions = getattr(model, "reactions", [])
    return list(reactions)


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


def _get_lb(reaction: object) -> float:
    return float(getattr(reaction, "lower_bound"))


def _get_ub(reaction: object) -> float:
    return float(getattr(reaction, "upper_bound"))


def _set_lb(reaction: object, value: float) -> None:
    setattr(reaction, "lower_bound", float(value))


def _set_ub(reaction: object, value: float) -> None:
    setattr(reaction, "upper_bound", float(value))
