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
  （查询 up 基因在逆转参考谱中位次靠后 → 正分），[0,1] 位次缩放。
- GCT 必须 `#1.2` 头；共同基因 <10 → 退出 2（拒绝在零交集上算相关）。
- 分块读 GCT（chunk 行×全列），不整块稠密化。

## 解释边界（log 明示）

表达逆转≠药效；非肌肉细胞系为跨模型外推，需在解读层降级（§18.2）。

## 冒烟记录

合成同向谱 spearman=+1.0、反向谱 −1.0；wcs 反向 +0.394 / 同向 −0.394 /
随机 ≈0；随机谱 n_common 一致。
