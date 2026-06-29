#!/usr/bin/env bash
set -euo pipefail

syncomdesign run \
  --config config/syncomdesign_config.yml \
  --models models \
  --medium media/medium.tsv \
  --outdir results \
  --solver glpk \
  --threads 1
