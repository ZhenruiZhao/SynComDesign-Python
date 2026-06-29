# GitHub Upload Checklist

Run these commands from the minimal repository directory:

```bash
cd /path/to/SynComDesign-Python-github-minimal
pytest -q
python -m compileall syncomdesign
```

Remove generated caches before committing:

```bash
find . -type d -name "__pycache__" -prune -exec rm -rf {} +
rm -rf .pytest_cache
```

Initialize and push:

```bash
git init
git add .
git status
git commit -m "Initial MATLAB-aligned SynComDesign Python CLI"
git branch -M main
git remote add origin https://github.com/YOUR_ORG/YOUR_REPO.git
git push -u origin main
```

Before pushing, confirm that no large model or result files are staged:

```bash
git status --short
git ls-files | grep -E '(^results|python_reference_exports|\\.xml$|\\.mat$|\\.sbml$)' || true
```

Real model files should stay outside GitHub and be supplied by each user through `--models`.
