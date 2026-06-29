from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

from .io import write_tsv


@dataclass
class Diagnostics:
    outdir: Path
    warnings: list[dict[str, object]] = field(default_factory=list)
    failed_combinations: list[dict[str, object]] = field(default_factory=list)

    def log(self, message: str) -> None:
        self.outdir.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with (self.outdir / "run.log").open("a", encoding="utf-8") as handle:
            handle.write(f"[{stamp}] {message}\n")

    def warn(self, table: str, message: str, **extra: object) -> None:
        row = {"table": table, "warning": message}
        row.update(extra)
        self.warnings.append(row)
        self.log(f"WARNING {table}: {message}")

    def fail_combination(self, combination_id: str, error_message: str, runtime_seconds: float = 0.0) -> None:
        self.failed_combinations.append(
            {
                "combination_id": combination_id,
                "error_message": error_message,
                "runtime_seconds": runtime_seconds,
            }
        )
        self.log(f"Combination failed: {combination_id} :: {error_message}")

    def write(self) -> None:
        write_tsv(self.outdir / "warnings.tsv", self.warnings, ["table", "warning", "combination_id", "reaction_id"])
        write_tsv(
            self.outdir / "failed_combinations.tsv",
            self.failed_combinations,
            ["combination_id", "error_message", "runtime_seconds"],
        )

