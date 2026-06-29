from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any

try:
    import yaml
except Exception:  # pragma: no cover - exercised when PyYAML is unavailable
    yaml = None


DEFAULT_CONFIG: dict[str, Any] = {
    "project": {"name": "SynComDesign", "output_dir": "results", "tmp_dir": None},
    "models": {
        "directory": "models",
        "file_pattern": "*.xml",
        "biomass_reactions_file": "config/biomass_reactions.tsv",
        "metabolite_aliases_file": "config/metabolite_aliases.tsv",
        "preserve_strain_id_as_string": True,
    },
    "medium": {
        "file": "media/medium.tsv",
        "condition": "anaerobic",
        "community_medium_mode": "external_shared_only",
        "close_unlisted_external_medium_uptakes": True,
        "allow_cross_feeding": True,
        "close_strain_interface_uptakes": False,
        "close_internal_transport": False,
        "legacy_all_exchange": False,
    },
    "combinations": {"mode": "all", "min_size": 1, "max_size": None},
    "objective": {"scenario_id": 1, "growth_fraction": 0.9, "target_strain": None, "biomass_weights": "equal"},
    "community": {"require_all_species_active": True, "minimum_biomass_flux": 1e-6, "shared_environment_compartment": "u"},
    "solver": {"name": "glpk", "tolerance": 1e-9},
    "output": {"write_debug_tables": True, "write_failed_combinations": True, "write_community_models": False},
}

SCENARIO_TYPES = {
    1: "total_biomass",
    2: "target_strain_biomass",
    3: "equal_composition",
    4: "fixed_composition",
    5: "growth_then_n2o_consumption",
}


def load_config(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    if yaml is None:
        data = _parse_simple_yaml(path.read_text(encoding="utf-8-sig"))
    else:
        data = yaml.safe_load(path.read_text(encoding="utf-8-sig")) or {}
    cfg = deepcopy(DEFAULT_CONFIG)
    _deep_update(cfg, data)
    scenario = cfg.get("objective", {}).get("scenario_id")
    if scenario is not None:
        cfg["objective"]["type"] = SCENARIO_TYPES[int(scenario)]
    else:
        cfg["objective"].setdefault("type", "total_biomass")
    cfg["_config_path"] = str(path)
    cfg["_base_dir"] = str(path.parent.parent if path.parent.name == "config" else path.parent)
    return cfg


def apply_overrides(cfg: dict[str, Any], **overrides: Any) -> dict[str, Any]:
    cfg = deepcopy(cfg)
    if overrides.get("models"):
        cfg["models"]["directory"] = overrides["models"]
    if overrides.get("medium"):
        cfg["medium"]["file"] = overrides["medium"]
    if overrides.get("outdir"):
        cfg["project"]["output_dir"] = overrides["outdir"]
    if overrides.get("solver"):
        cfg["solver"]["name"] = overrides["solver"]
    if overrides.get("tmpdir"):
        cfg["project"]["tmp_dir"] = overrides["tmpdir"]
    if overrides.get("threads") is not None:
        cfg["solver"]["threads"] = overrides["threads"]
    return cfg


def _deep_update(target: dict[str, Any], source: dict[str, Any]) -> None:
    for key, value in source.items():
        if isinstance(value, dict) and isinstance(target.get(key), dict):
            _deep_update(target[key], value)
        else:
            target[key] = value


def _parse_simple_yaml(text: str) -> dict[str, Any]:
    data: dict[str, Any] = {}
    current: dict[str, Any] | None = None
    for raw_line in text.splitlines():
        line = raw_line.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        if not line.startswith(" ") and line.endswith(":"):
            section = line[:-1].strip()
            current = data.setdefault(section, {})
            continue
        if ":" not in line:
            continue
        key, raw = line.strip().split(":", 1)
        value = _parse_scalar(raw.strip())
        if current is None:
            data[key] = value
        else:
            current[key] = value
    return data


def _parse_scalar(raw: str) -> Any:
    if raw == "" or raw.lower() == "null":
        return None
    if raw.lower() == "true":
        return True
    if raw.lower() == "false":
        return False
    if raw.startswith("[") and raw.endswith("]"):
        body = raw[1:-1].strip()
        if not body:
            return []
        return [_parse_scalar(item.strip()) for item in body.split(",")]
    try:
        if any(ch in raw for ch in [".", "e", "E"]):
            return float(raw)
        return int(raw)
    except ValueError:
        return raw.strip("'\"")

