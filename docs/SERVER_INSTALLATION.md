# Server Installation Guide

This guide installs SynComDesign once on a Linux server so multiple users can run their own configs, model folders, medium files, and output directories.

## 1. Prepare Conda

Install Miniconda or use an existing cluster Conda installation.

```bash
source /share/software/miniconda3/etc/profile.d/conda.sh
```

Create the environment:

```bash
cd /share/software/SynComDesign-Python
conda env create -f environment.yml
conda activate syncomdesign
pip install -e .
```

If the environment already exists:

```bash
conda activate syncomdesign
pip install -e .
```

Check the command:

```bash
syncomdesign --help
syncomdesign validate --config config/syncomdesign_config.yml
```

## 2. Recommended Shared Layout

Install the software once:

```text
/share/software/SynComDesign-Python/
  syncomdesign/
  config/
  media/
  docs/
  scripts/
  pyproject.toml
  environment.yml
```

Each user keeps their own project:

```text
/share/home/$USER/syncomdesign_project/
  config/
    syncomdesign_config.yml
    biomass_reactions.tsv
    metabolite_aliases.tsv
  media/
    medium.tsv
  models/
    005.xml
    016.xml
    ...
  results/
```

Do not write user results into the installed package directory.

## 3. User Project Setup

```bash
mkdir -p /share/home/$USER/syncomdesign_project/{config,media,models,results}
cp /share/software/SynComDesign-Python/config/*.tsv /share/home/$USER/syncomdesign_project/config/
cp /share/software/SynComDesign-Python/config/syncomdesign_config.yml /share/home/$USER/syncomdesign_project/config/
cp /share/software/SynComDesign-Python/media/medium.tsv /share/home/$USER/syncomdesign_project/media/
```

Copy SBML models into:

```text
/share/home/$USER/syncomdesign_project/models/
```

Edit:

```text
/share/home/$USER/syncomdesign_project/config/biomass_reactions.tsv
/share/home/$USER/syncomdesign_project/media/medium.tsv
/share/home/$USER/syncomdesign_project/config/syncomdesign_config.yml
```

The config paths are resolved relative to the config file location unless command-line overrides are used.

## 4. Choosing Objective Modes (ID1-ID5)

Edit the `objective` section in:

```text
/share/home/$USER/syncomdesign_project/config/syncomdesign_config.yml
```

Base template:

```yaml
objective:
  scenario_id: 1
  growth_fraction: 0.9
  target_strain: null
  biomass_weights: equal
```

Available modes:

| ID | Meaning | Server use |
| --- | --- | --- |
| ID1 | Maximize total community biomass. | Default screening mode. |
| ID2 | Maximize target-strain biomass. Only combinations containing `target_strain` are evaluated. | Use when one strain must be present or promoted. |
| ID3 | Equal community composition. | Use when all selected strains should grow at equal biomass proportions. |
| ID4 | Fixed community composition. | Use when a known composition ratio should be enforced. |
| ID5 | Growth-first, N2O-consumption-second. | Use when growth must be retained while prioritizing N2O uptake. |

Example ID1 config:

```yaml
objective:
  scenario_id: 1
  growth_fraction: 0.9
  target_strain: null
  biomass_weights: equal
```

Example ID2 config for strain `005`:

```yaml
objective:
  scenario_id: 2
  target_strain: "005"
  growth_fraction: 0.9
  biomass_weights: equal
```

Example ID5 config:

```yaml
objective:
  scenario_id: 5
  growth_fraction: 0.9
  target_strain: null
  biomass_weights: equal
```

Recommended project layout for multiple scenarios:

```bash
cp $PROJECT/config/syncomdesign_config.yml $PROJECT/config/id1.yml
cp $PROJECT/config/syncomdesign_config.yml $PROJECT/config/id2_005.yml
cp $PROJECT/config/syncomdesign_config.yml $PROJECT/config/id5.yml
```

Run each scenario into a separate output directory:

```bash
syncomdesign run \
  --config $PROJECT/config/id2_005.yml \
  --models $PROJECT/models \
  --medium $PROJECT/media/medium.tsv \
  --outdir $PROJECT/results_id2_005 \
  --solver glpk \
  --threads 1
```

## 5. Run Locally on a Login Node

For a small smoke test:

```bash
conda activate syncomdesign
syncomdesign validate \
  --config /share/home/$USER/syncomdesign_project/config/syncomdesign_config.yml
```

Run ID1:

```bash
syncomdesign run \
  --config /share/home/$USER/syncomdesign_project/config/syncomdesign_config.yml \
  --models /share/home/$USER/syncomdesign_project/models \
  --medium /share/home/$USER/syncomdesign_project/media/medium.tsv \
  --outdir /share/home/$USER/syncomdesign_project/results_id1 \
  --solver glpk \
  --threads 1
```

## 6. PBS Job Script

Save as `run_syncomdesign.pbs`:

```bash
#!/bin/bash
#PBS -N syncomdesign
#PBS -l nodes=node1:ppn=20
#PBS -j oe

source /share/software/miniconda3/etc/profile.d/conda.sh
conda activate syncomdesign

PROJECT=/share/home/$USER/syncomdesign_project

syncomdesign run \
  --config $PROJECT/config/syncomdesign_config.yml \
  --models $PROJECT/models \
  --medium $PROJECT/media/medium.tsv \
  --outdir $PROJECT/results_id1 \
  --solver glpk \
  --threads 1
```

Submit:

```bash
qsub run_syncomdesign.pbs
```

Monitor:

```bash
qstat -u $USER
tail -f syncomdesign.o*
```

## 7. MATLAB Alignment Check

If MATLAB reference exports are available:

```bash
syncomdesign compare-matlab \
  --config /share/home/$USER/syncomdesign_project/config/syncomdesign_config.yml \
  --python-outdir /share/home/$USER/syncomdesign_project/results_id1 \
  --reference /share/home/$USER/python_reference_exports \
  --outdir /share/home/$USER/syncomdesign_project/results_compare_matlab
```

Check:

```text
results_compare_matlab/matlab_alignment_report.md
results_compare_matlab/matlab_alignment_summary.tsv
results_compare_matlab/matlab_alignment_differences.tsv
```

Priority checks:

1. `medium_mapping`
2. `external_medium_bounds`
3. `interface_bounds`
4. `internal_transport_bounds`
5. `community_summary`

If bounds differ, fix mapping/bounds first. Do not explain bounds mismatches as solver differences.

## 8. Zero-Biomass Diagnostics

If all biomass values are unexpectedly zero:

```bash
syncomdesign diagnose-zero \
  --config /share/home/$USER/syncomdesign_project/config/syncomdesign_config.yml \
  --models /share/home/$USER/syncomdesign_project/models \
  --medium /share/home/$USER/syncomdesign_project/media/medium.tsv \
  --outdir /share/home/$USER/syncomdesign_project/results_id1
```

Diagnostics are written to:

```text
results_id1/debug_zero_fix/
```

## 9. Important Medium Rules

- Medium only changes `external_medium_exchange`.
- Unlisted external shared uptake is closed.
- Strain-interface reactions are not closed by medium.
- Internal transport reactions are not closed by medium.
- Cross-feeding is allowed through shared-pool mass balance.
- COBRApy `model.medium` is not used for the community medium.

## 10. Common Problems

No models detected:

```bash
syncomdesign validate --config /path/to/config/syncomdesign_config.yml
```

Check `models.directory`, `models.file_pattern`, and `--models`.

All combinations failed:

Check `failed_combinations.tsv` and `run.log` in the output directory.

All biomass values are zero:

Run `syncomdesign diagnose-zero` and inspect `debug_zero_fix/`.

MATLAB alignment is PARTIAL:

First check whether medium mapping and bounds pass. Fluxes can differ under degenerate optima even when objective and bounds match.
