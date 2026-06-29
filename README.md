# SynComDesign Python

MATLAB-aligned Python CLI for SynComDesign community metabolic design.

This repository is the minimal GitHub version. It contains source code, configuration templates, tests, scripts, and documentation. Large user model files and run outputs are intentionally excluded.

## What Is Aligned

The Python implementation follows the fixed MATLAB reference logic:

- Medium only applies to `external_medium_exchange`.
- Unlisted external shared uptake is closed.
- Strain-interface reactions are not closed by medium.
- Internal transport reactions are not closed by medium.
- Cross-feeding is allowed through shared-pool mass balance.
- All non-empty strain combinations are enumerated by default.
- ID2 target-strain mode is the only mode that filters by target strain.
- NO, NO2, NO3, N2O, and N2 fluxes use explicit aliases, not substring matching.
- SBML `boundaryCondition=true` species are handled to match MATLAB `readCbModel` behavior.

## Install

```bash
conda env create -f environment.yml
conda activate syncomdesign
pip install -e .
```

Check:

```bash
syncomdesign --help
pytest -q
```

## Inputs

Prepare a user project with:

```text
project/
  config/
    syncomdesign_config.yml
    biomass_reactions.tsv
    metabolite_aliases.tsv
  media/
    medium.tsv
  models/
    *.xml
```

Real model files are not committed to this GitHub repository. Put them in `models/` or pass them with `--models`.

## Run

```bash
syncomdesign run \
  --config config/syncomdesign_config.yml \
  --models models \
  --medium media/medium.tsv \
  --outdir results_id1 \
  --solver glpk \
  --threads 1
```

## Objective Modes (ID1-ID5)

Choose the design objective in `config/syncomdesign_config.yml`:

```yaml
objective:
  scenario_id: 1
  growth_fraction: 0.9
  target_strain: null
  biomass_weights: equal
```

Set `objective.scenario_id` to one of:

| ID | Mode | Use case | Extra settings |
| --- | --- | --- | --- |
| ID1 | Total community biomass | Find combinations with maximum total community growth. | None |
| ID2 | Target-strain biomass | Find combinations that maximize one selected strain. Only combinations containing the target strain are evaluated. | Set `target_strain`, for example `"005"` |
| ID3 | Equal community composition | Force strain biomass fluxes to equal proportions while optimizing growth. | None |
| ID4 | Fixed community composition | Force strain biomass fluxes to configured proportions. If no ratio is provided, equal proportions are used to match MATLAB behavior. | Optional `composition_ratio` |
| ID5 | Growth first, N2O consumption second | First maximize growth, then keep growth above `growth_fraction * max_growth` and maximize N2O uptake. | Set `growth_fraction` |

Examples:

ID1:

```yaml
objective:
  scenario_id: 1
  growth_fraction: 0.9
  target_strain: null
  biomass_weights: equal
```

ID2 for strain `005`:

```yaml
objective:
  scenario_id: 2
  target_strain: "005"
  growth_fraction: 0.9
  biomass_weights: equal
```

ID5:

```yaml
objective:
  scenario_id: 5
  growth_fraction: 0.9
  target_strain: null
  biomass_weights: equal
```

A practical server pattern is to keep one config file per objective:

```bash
cp config/syncomdesign_config.yml config/id1.yml
cp config/syncomdesign_config.yml config/id2_005.yml
cp config/syncomdesign_config.yml config/id5.yml
```

Then edit `scenario_id`, `target_strain`, and `growth_fraction`, and write results to separate output directories:

```bash
syncomdesign run --config config/id2_005.yml --outdir results_id2_005
```

Validate inputs:

```bash
syncomdesign validate --config config/syncomdesign_config.yml
```

Diagnose unexpected zero biomass:

```bash
syncomdesign diagnose-zero \
  --config config/syncomdesign_config.yml \
  --models models \
  --medium media/medium.tsv \
  --outdir results_id1
```

## Outputs

Main output tables include:

- `community_summary.tsv`
- `objective_trace.tsv`
- `all_combinations.tsv`
- `failed_combinations.tsv`
- `medium_to_shared_exchange_mapping.tsv`
- `medium_mapping_warnings.tsv`
- `external_medium_bounds.tsv`
- `interface_bounds.tsv`
- `internal_transport_bounds.tsv`
- `reaction_classification.tsv`
- `community_build_trace.tsv`
- `flux_mapping.tsv`
- `flux_values.tsv`

## MATLAB Alignment

If MATLAB reference exports are available:

```bash
syncomdesign compare-matlab \
  --config config/syncomdesign_config.yml \
  --python-outdir results_id1 \
  --reference /path/to/python_reference_exports \
  --outdir results_compare_matlab
```

Inspect:

```text
results_compare_matlab/matlab_alignment_report.md
results_compare_matlab/matlab_alignment_summary.tsv
results_compare_matlab/matlab_alignment_differences.tsv
```

## Server Use

See [docs/SERVER_INSTALLATION.md](docs/SERVER_INSTALLATION.md) for a detailed multi-user Linux/PBS installation guide.

## Repository Hygiene

Do not commit:

- real SBML/COBRA model files,
- `results*/`,
- `python_reference_exports/`,
- solver logs,
- Python caches.

These are already excluded in `.gitignore`.
