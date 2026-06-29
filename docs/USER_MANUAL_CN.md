# SynComDesign Python 用户手册

## 目标

SynComDesign Python 是基于 SynComDesign MATLAB 最终修复版实现的命令行工具，适合在无图形界面的服务器上运行。

Python 版以 MATLAB reference export 为金标准，不使用 COBRApy `model.medium` 直接替代 community medium 逻辑。

## 输入准备

每个用户准备自己的项目目录：

```text
project/
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

菌株 ID 按文件名字符串处理，例如 `005.xml` 的菌株 ID 是 `005`，不会变成数字 `5`。

## 安装

```bash
conda env create -f environment.yml
conda activate syncomdesign
pip install -e .
```

## 检查输入

```bash
syncomdesign validate --config config/syncomdesign_config.yml
```

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

服务器多人使用时，每个用户应使用自己的 config、models、medium 和 outdir，不要把结果写入软件安装目录。

## Objective Modes

- ID1：最大化 community total biomass。
- ID2：最大化 target strain biomass，只评估包含 target strain 的组合。
- ID3：equal community composition。
- ID4：fixed community composition；未提供比例时按 MATLAB 行为使用等比例。
- ID5：先最大化 growth，再在 growth fraction 约束下最大化 N2O uptake。

## Medium 与 Cross-Feeding

- medium 只限制 external shared exchange。
- medium 给出的物质才允许从外界摄取。
- 未列入 medium 的 external shared uptake 设为 0。
- strain-interface 不被 medium 关闭。
- internal transport 不被 medium 关闭。
- 菌株分泌到 shared pool 的物质可以被其他菌株摄取。
- cross-feeding 由 shared pool 质量守恒保证。

## 结果文件

主要结果：

```text
community_summary.tsv
objective_trace.tsv
flux_values.tsv
failed_combinations.tsv
```

调试和对齐文件：

```text
all_combinations.tsv
community_build_trace.tsv
reaction_classification.tsv
medium_to_shared_exchange_mapping.tsv
medium_mapping_warnings.tsv
external_medium_bounds.tsv
interface_bounds.tsv
internal_transport_bounds.tsv
flux_mapping.tsv
```

如果 `failed_combinations.tsv` 有内容，先查看具体组合错误；程序会尽量继续运行后续组合。


## Zero Biomass 诊断

如果全部组合 biomass 异常为 0：

```bash
syncomdesign diagnose-zero \
  --config config/syncomdesign_config.yml \
  --models models \
  --medium media/medium.tsv \
  --outdir results_id1
```

诊断表位于：

```text
results_id1/debug_zero_fix/
```
