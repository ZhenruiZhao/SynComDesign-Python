from pathlib import Path
import os

import pytest

from syncomdesign.io import read_tsv
from syncomdesign.config import apply_overrides, load_config
from syncomdesign.validation import load_models


REFERENCE_ENV = os.environ.get("SYNCOMDESIGN_MATLAB_REFERENCE")
MODEL_DIR_ENV = os.environ.get("SYNCOMDESIGN_MODEL_DIR")
REFERENCE = Path(REFERENCE_ENV) if REFERENCE_ENV else None
MODEL_DIR = Path(MODEL_DIR_ENV) if MODEL_DIR_ENV else None


pytestmark = pytest.mark.skipif(
    REFERENCE is None or MODEL_DIR is None or not REFERENCE.exists() or not MODEL_DIR.exists(),
    reason="set SYNCOMDESIGN_MATLAB_REFERENCE and SYNCOMDESIGN_MODEL_DIR to run MATLAB alignment tests",
)


def test_matlab_reference_bounds_match():
    rows = read_tsv(REFERENCE / "reference_external_medium_bounds.tsv")
    assert rows
    assert all(row["reaction_class"] == "external_medium_exchange" for row in rows[:100])
    assert any(row["reaction_id"].startswith("R_EX_no3_") for row in rows)


def test_matlab_reference_summary_match():
    rows = read_tsv(REFERENCE / "reference_community_summary.tsv")
    assert len(rows) == 31
    first = rows[0]
    assert first["combination_id"] == "005"
    assert "total_biomass" in first


def test_python_loader_removes_sbml_boundary_species_like_matlab():
    project = Path(__file__).resolve().parents[1]
    cfg = apply_overrides(
        load_config(project / "config" / "syncomdesign_config.yml"),
        models=str(MODEL_DIR),
    )
    infos, _rows = load_models(cfg, project)
    model_005 = next(info["model"] for info in infos if info["name"] == "005")

    assert len(model_005.reactions.get_by_id("R_EX_no3_e").metabolites) == 0
