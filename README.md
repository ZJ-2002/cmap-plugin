# cmap

LINCS/CMap 连通性评分（主计划 §18.2）。查询签名对 GCT 参考谱矩阵的
全量评分排名，喂下游逆转候选筛选。

## 节点

| kind | 输入 | 输出 |
|---|---|---|
| `cmap_connectivity` | 查询签名 TSV（gene,value）+ GCT | `cmap_scores.tsv`（signature_id/score/n_common/a_up/a_down）、log |

参数：`metric`（spearman|pearson|wcs，默认 spearman）、`gene_col`/
`value_col`、`up_n`/`down_n`（默认 150，仅 wcs 用）、`chunk`（默认 2000）。

## 布局

```text
cmap/
├── manifest.toml
├── scripts/cmap_connectivity.py
├── Dockerfile
└── README.md
```

镜像：`localhost/autonomics/biotools-py@sha256:21e14c1582a9d4c261c0a23a11864293313ee2b49b9b7d85bcabc89210b48c7e`。

## 口径（§18.2 硬约束）

- **完整排名输出，不做条件筛选**：细胞类型/剂量/时间过滤在下游注册
  表完成，本节点不挑"最强逆转"（挑峰值=多重比较里偷看答案）。
- 三种度量并列提供：spearman/pearson 全基因相关；wcs = `a_down − a_up`
  （up/down 集按查询**带符号值**取最正/最负前 up_n/down_n，v1 按 |value|
  排是错误口径已修；查询 up 基因在逆转参考谱中位次靠后 → 正分），
  [0,1] 位次缩放；up/down 集重叠时 WARN。
- 参考矩阵按魔数自动识别：GCT 文本（#1.2 头）/ GCTx（/0/DATA/0）/ 各自
  的 .gz 形态；gctx 形状与元数据不符 → 退出 2。
- 共同基因 <10 → 退出 2（拒绝在零交集上算相关；LINCS gctx 行是 Entrez
  ID，先用 gene_info 映射查询符号，geo_suppl 可取）。
- 按签名【列】分块（完整基因集 × 列块，v2 修正：v1 按基因行分块会把
  同一签名拆成多个子集相关的行），不整块稠密化。

## 解释边界（log 明示）

表达逆转≠药效；非肌肉细胞系为跨模型外推，需在解读层降级（§18.2）。

## 冒烟记录（2026-10-04，biotools-py 容器内实跑）

50 基因 × 7 谱合成矩阵（sig0=查询完全同向、sig1=完全逆向），四格全验：

- gct/spearman：sig0=+1.0 居首、sig1=−1.0 垫底；gctx/spearman 同
  （`source=gctx`，HDF5 /0/DATA/0 路径）。
- wcs（up_n=down_n=10，带符号取集）：sig1=+0.816 居首（a_up=0.09，
  查询上调基因在逆向谱中位次垫底）、sig0=−0.816 垫底；gct/gctx 一致。
- v1 缺陷复现记录：按 |value| 排且 up_n+down_n ≥ 共同基因数时两集
  全重叠、全部得 0 分（退化），此为 v2 重写动机。
