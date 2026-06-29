# SynComDesign Python 中文说明

这是严格对齐 SynComDesign MATLAB 最终修复版的 Python CLI 最小 GitHub 版本。

本仓库只保留源码、配置模板、测试、脚本和文档。真实模型文件和运行结果不上传 GitHub。

## 对齐原则

Python 版本遵循 MATLAB reference 的 community medium 逻辑：

- medium 只作用于 `external_medium_exchange`；
- medium 未列出的 external shared uptake 会被关闭；
- strain-interface 不被 medium 关闭；
- internal transport 不被 medium 关闭；
- cross-feeding 由 shared pool 质量守恒允许；
- 默认枚举所有非空菌株组合；
- 只有 ID2 target-strain mode 会按 target strain 过滤组合；
- NO、NO2、NO3、N2O、N2 使用 alias 表严格映射，不用模糊 substring；
- SBML 中 `boundaryCondition=true` 的外界 species 按 MATLAB `readCbModel` 行为处理。

## 安装

```bash
conda env create -f environment.yml
conda activate syncomdesign
pip install -e .
```

检查：

```bash
syncomdesign --help
pytest -q
```

## 输入目录

建议每个用户准备自己的项目目录：

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

真实模型文件请放入 `models/`，或运行时用 `--models` 指定。

## 运行

```bash
syncomdesign run \
  --config config/syncomdesign_config.yml \
  --models models \
  --medium media/medium.tsv \
  --outdir results_id1 \
  --solver glpk \
  --threads 1
```

检查输入：

```bash
syncomdesign validate --config config/syncomdesign_config.yml
```

如果出现异常的全部 biomass 为 0：

```bash
syncomdesign diagnose-zero \
  --config config/syncomdesign_config.yml \
  --models models \
  --medium media/medium.tsv \
  --outdir results_id1
```

诊断表会写入：

```text
results_id1/debug_zero_fix/
```

## MATLAB 对齐

如果已有 MATLAB reference exports：

```bash
syncomdesign compare-matlab \
  --config config/syncomdesign_config.yml \
  --python-outdir results_id1 \
  --reference /path/to/python_reference_exports \
  --outdir results_compare_matlab
```

重点查看：

```text
results_compare_matlab/matlab_alignment_report.md
results_compare_matlab/matlab_alignment_summary.tsv
results_compare_matlab/matlab_alignment_differences.tsv
```

如果 medium bounds 或 interface bounds 不一致，优先修复 mapping/bounds，不要直接归因于 solver 差异。

## 服务器安装

详细教程见：

```text
docs/SERVER_INSTALLATION.md
```

包括多人共享安装、用户项目目录、PBS 脚本、MATLAB alignment 和 zero-biomass 诊断。
