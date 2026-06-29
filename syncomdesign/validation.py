from __future__ import annotations

from pathlib import Path
from typing import Mapping
import xml.etree.ElementTree as ET
import logging

from .io import read_tsv


def load_biomass_reactions(path: str | Path) -> dict[str, str]:
    rows = read_tsv(path)
    out = {}
    for row in rows:
        strain = str(row.get("strain", ""))
        biomass = str(row.get("biomass_rxn", ""))
        if strain:
            out[strain] = biomass
    return out


def detect_model_files(directory: str | Path, pattern: str = "*.xml") -> list[dict[str, str]]:
    directory = Path(directory)
    return [{"strain": path.stem, "model_path": str(path)} for path in sorted(directory.glob(pattern))]


def load_models(config: Mapping[str, object], base_dir: str | Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    try:
        from cobra.io import read_sbml_model
    except Exception as exc:  # pragma: no cover - depends on optional COBRApy install
        raise RuntimeError("COBRApy is required to load SBML models") from exc
    logging.getLogger("cobra").setLevel(logging.ERROR)

    base_dir = Path(base_dir)
    models_cfg = config["models"]
    model_dir = _resolve(base_dir, models_cfg["directory"])
    biomass_file = _resolve(base_dir, models_cfg["biomass_reactions_file"])
    biomass_map = load_biomass_reactions(biomass_file)
    model_infos = []
    validation_rows = []
    for row in detect_model_files(model_dir, str(models_cfg.get("file_pattern", "*.xml"))):
        strain = str(row["strain"])
        path = Path(row["model_path"])
        try:
            original_reactions = _sbml_reaction_ids(path)
            boundary_species = _sbml_boundary_species_ids(path)
            model = read_sbml_model(str(path))
            _remove_boundary_species_stoichiometry(model, boundary_species)
            _remove_non_original_reactions(model, original_reactions)
            _restore_original_reaction_ids(model, original_reactions)
            biomass_rxn = biomass_map.get(strain) or detect_biomass_reaction(model)
            model_infos.append({"name": strain, "path": str(path), "model": model, "biomass_rxn": biomass_rxn})
            validation_rows.append({"strain": strain, "model_path": str(path), "valid": True, "biomass_rxn": biomass_rxn, "warning_message": ""})
        except Exception as exc:
            validation_rows.append({"strain": strain, "model_path": str(path), "valid": False, "biomass_rxn": "", "warning_message": str(exc)})
    return model_infos, validation_rows


def detect_biomass_reaction(model: object) -> str:
    candidates = []
    for reaction in getattr(model, "reactions", []):
        rid = reaction.id.lower()
        if rid == "growth" or "biomass" in rid or "growth" in rid:
            candidates.append(reaction.id)
    if not candidates:
        raise ValueError("no biomass reaction found; provide config/biomass_reactions.tsv")
    return candidates[0]


def validate_project_inputs(config: Mapping[str, object], base_dir: str | Path) -> list[dict[str, object]]:
    base_dir = Path(base_dir)
    checks = []
    models_cfg = config["models"]
    medium_cfg = config["medium"]
    for label, path in [
        ("models.directory", _resolve(base_dir, models_cfg["directory"])),
        ("models.biomass_reactions_file", _resolve(base_dir, models_cfg["biomass_reactions_file"])),
        ("models.metabolite_aliases_file", _resolve(base_dir, models_cfg["metabolite_aliases_file"])),
        ("medium.file", _resolve(base_dir, medium_cfg["file"])),
    ]:
        checks.append({"check": label, "path": str(path), "exists": path.exists()})
    return checks


def _resolve(base_dir: Path, value: object) -> Path:
    path = Path(str(value))
    return path if path.is_absolute() else base_dir / path


def _sbml_reaction_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    for _event, elem in ET.iterparse(path, events=("start",)):
        if elem.tag.endswith("reaction"):
            reaction_id = elem.attrib.get("id")
            if reaction_id:
                ids.add(reaction_id)
        elem.clear()
    return ids


def _sbml_boundary_species_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    for _event, elem in ET.iterparse(path, events=("start",)):
        if elem.tag.endswith("species") and elem.attrib.get("boundaryCondition", "").lower() == "true":
            species_id = elem.attrib.get("id")
            if species_id:
                ids.update(_species_id_variants(species_id))
        elem.clear()
    return ids


def _species_id_variants(species_id: str) -> set[str]:
    variants = {species_id}
    if species_id.startswith("M_"):
        variants.add(species_id[2:])
    else:
        variants.add(f"M_{species_id}")
    return variants


def _remove_boundary_species_stoichiometry(model: object, boundary_species_ids: set[str]) -> None:
    """Match MATLAB readCbModel behavior for SBML boundaryCondition species.

    COBRApy materializes boundary species inside exchange and transport
    reactions. The MATLAB reference readCbModel drops those species from the
    stoichiometric matrix, so the Python model must remove the same rows before
    community construction.
    """
    if not boundary_species_ids:
        return
    changed = False
    for reaction in list(getattr(model, "reactions", [])):
        removals = {}
        for metabolite, coefficient in list(getattr(reaction, "metabolites", {}).items()):
            met_id = getattr(metabolite, "id", str(metabolite))
            if met_id in boundary_species_ids:
                removals[metabolite] = -coefficient
        if removals:
            reaction.add_metabolites(removals)
            changed = True
    if changed and hasattr(model, "repair"):
        model.repair()


def _remove_non_original_reactions(model: object, original_ids: set[str]) -> None:
    if not original_ids:
        return
    to_remove = []
    for reaction in list(getattr(model, "reactions", [])):
        rid = reaction.id
        variants = {rid}
        if not rid.startswith("R_"):
            variants.add(f"R_{rid}")
        if rid.startswith("R_"):
            variants.add(rid[2:])
        if variants.isdisjoint(original_ids):
            to_remove.append(reaction)
    if to_remove:
        model.remove_reactions(to_remove, remove_orphans=False)


def _restore_original_reaction_ids(model: object, original_ids: set[str]) -> None:
    changed = False
    for reaction in list(getattr(model, "reactions", [])):
        rid = reaction.id
        if rid in original_ids:
            continue
        prefixed = f"R_{rid}"
        if prefixed in original_ids:
            reaction.id = prefixed
            changed = True
    if changed and hasattr(model, "repair"):
        model.repair()
