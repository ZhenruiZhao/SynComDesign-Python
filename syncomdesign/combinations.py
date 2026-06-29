from __future__ import annotations

from itertools import combinations
from pathlib import Path
from typing import Iterable

from .io import write_tsv


def enumerate_all(
    strains: Iterable[str],
    min_size: int = 1,
    max_size: int | None = None,
    target_strain: str | None = None,
    objective_mode: str | int | None = None,
) -> list[tuple[str, ...]]:
    strain_list = [str(strain) for strain in strains]
    if len(set(strain_list)) != len(strain_list):
        raise ValueError("strain IDs must be unique")
    if max_size is None:
        max_size = len(strain_list)
    min_size = max(1, int(min_size))
    max_size = min(int(max_size), len(strain_list))
    require_target = _is_id2(objective_mode) and target_strain not in {None, ""}
    out: list[tuple[str, ...]] = []
    for size in range(min_size, max_size + 1):
        for combo in combinations(strain_list, size):
            if require_target and str(target_strain) not in combo:
                continue
            out.append(combo)
    return out


def write_all_combinations(path: str | Path, combos: Iterable[tuple[str, ...]]) -> list[dict[str, object]]:
    rows = []
    for index, combo in enumerate(combos, start=1):
        rows.append(
            {
                "combination_index": index,
                "combination_id": "+".join(combo),
                "community_size": len(combo),
                "strain_names": ";".join(combo),
            }
        )
    write_tsv(path, rows, ["combination_index", "combination_id", "community_size", "strain_names"])
    return rows


def _is_id2(mode: str | int | None) -> bool:
    if mode is None:
        return False
    if isinstance(mode, int):
        return mode == 2
    normalized = str(mode).lower()
    return normalized in {"2", "id2", "target_strain_biomass"}

