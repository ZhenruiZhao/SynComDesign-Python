from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable, Mapping, Sequence


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    path = Path(path)
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        dialect = csv.excel_tab if path.suffix.lower() in {".tsv", ".txt"} else csv.excel
        return [dict(row) for row in csv.DictReader(handle, dialect=dialect)]


def write_tsv(path: str | Path, rows: Iterable[Mapping[str, object]], fieldnames: Sequence[str] | None = None) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    materialized = [dict(row) for row in rows]
    if fieldnames is None:
        fieldnames = []
        for row in materialized:
            for key in row:
                if key not in fieldnames:
                    fieldnames.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), dialect=csv.excel_tab, extrasaction="ignore")
        writer.writeheader()
        for row in materialized:
            writer.writerow({key: _format_value(row.get(key, "")) for key in fieldnames})


def _format_value(value: object) -> object:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    return value


def as_float(value: object, default: float = 0.0) -> float:
    if value is None or value == "":
        return default
    return float(value)


def resolve_path(base_dir: str | Path, path: str | Path | None) -> Path | None:
    if path is None or str(path) == "":
        return None
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return Path(base_dir) / candidate

