# Models Directory

Put user SBML/COBRA model files here, or pass another directory with:

```bash
syncomdesign run --config config/syncomdesign_config.yml --models /path/to/models
```

Large real model files are intentionally not committed to GitHub.

Supported file pattern is configured in `config/syncomdesign_config.yml`:

```yaml
models:
  directory: models
  file_pattern: "*.xml"
```
