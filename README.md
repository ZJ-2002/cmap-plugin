# cmap

LINCS/CMap 连通性评分（主计划 §18.2）。查询签名对 GCT/GCTx 参考谱矩阵的
全量评分排名，喂下游逆转候选筛选。

## 节点

| kind | 输入 | 输出 |
|---|---|---|
| `cmap_connectivity` | 查询签名 TSV（gene,value）+ GCT/GCTx | `cmap_scores.tsv`（signature_id/metric/score/direction/n_common/a_up/a_down）、log |

参数：`metric`（spearman|pearson|wcs，默认 spearman）、`gene_col`/
`value_col`、`up_n`/`down_n`（默认 150，仅 wcs 用）、`chunk`（默认 2000，
按签名列分块）。

## 布局

```text
cmap/
├── manifest.toml
├── scripts/cmap_connectivity.py
├── tests/test_cmap_connectivity.py  # 离线：官方结构 GCTx fixture + 已知答案
├── Dockerfile
└── README.md
```

镜像：`localhost/autonomics/biotools-py@sha256:21e14c1582a9d4c261c0a23a11864293313ee2b49b9b7d85bcabc89210b48c7e`。

## 口径（§18.2 硬约束）

- **完整排名输出，不做条件筛选**：细胞类型/剂量/时间过滤在下游注册
  表完成，本节点不挑"最强逆转"（挑峰值=多重比较里偷看答案）。
- 三种度量并列提供：spearman/pearson 全基因相关（正=同向、负=逆转）；
  `wcs` = **探索性位次逆转分**（exploratory rank-reversal score）=
  `a_down − a_up`（up/down 集按查询**带符号值**取最正/最负前 up_n/down_n，
  v1 按 |value| 排是错误口径已修；查询 up 基因在逆转参考谱中位次靠后 →
  正分），[0,1] 位次缩放；up/down 集重叠时 WARN。**wcs 是本节点自定义
  的 up/down 平均相对位次之差，不是官方 CMap WTCS/NCS/tau**——metric
  标识仅为向后兼容保留，引用时不得称官方连通性分。
- **GCTx 官方布局（v3，审计 F14）**：矩阵在 `/0/DATA/0/matrix`（v2 直接
  读 `/0/DATA/0`，那是 group——官方结构上 `.shape` 直接 AttributeError，
  v2 冒烟之所以通过是因为 fixture 照 bug 造）；官方存储方向
  (n_signature, n_gene)——COL（签名）轴 0、ROW（基因）轴 1，读入转置成
  gene×signature（cmapPy `parse_gctx` 同语义）；行/列 meta 取
  `/0/META/ROW/id`、`/0/META/COL/id`；形状与元数据不符 → 退出 2。
- **方向语义（v3 修正）**：结果表 `direction` 列逐行标注——相关类
  正=`same_direction`、负=`reversed`、零=`flat`；wcs 正=`reversed`、
  负=`concordant`。排序随 metric 逆转方向：spearman/pearson **升序**
  （最负=最强逆转在前）、wcs **降序**（最正=最强逆转在前）。v2 对全部
  metric 统一降序是错的（降序≠统一逆转——相关类的最强逆转在最负端），
  已废止。
- 参考矩阵按魔数自动识别：GCT 文本（#1.2 头）/ GCTx / 各自的 .gz 形态；
  gctx 形状与元数据不符 → 退出 2。
- 共同基因 <10 → 退出 2（拒绝在零交集上算相关；LINCS gctx 行是 Entrez
  ID，先用 gene_info 映射查询符号，geo_suppl 可取）。
- 按签名【列】分块（完整基因集 × 列块，v2 修正：v1 按基因行分块会把
  同一签名拆成多个子集相关的行），不整块稠密化。

## 解释边界（log 明示）

表达逆转≠药效；非肌肉细胞系为跨模型外推，需在解读层降级（§18.2）；
wcs 为探索性分数，不等价官方 WTCS/NCS/tau。

## 冒烟记录（2026-10-04，biotools-py 容器内实跑）

50 基因 × 7 谱合成矩阵（sig0=查询完全同向、sig1=完全逆向），四格全验：

- gct/spearman：sig0=+1.0 居首、sig1=−1.0 垫底；gctx/spearman 同
  （`source=gctx`）。
- wcs（up_n=down_n=10，带符号取集）：sig1=+0.816 居首（a_up=0.09，
  查询上调基因在逆向谱中位次垫底）、sig0=−0.816 垫底；gct/gctx 一致。
- v1 缺陷复现记录：按 |value| 排且 up_n+down_n ≥ 共同基因数时两集
  全重叠、全部得 0 分（退化），此为 v2 重写动机。
- **v3 勘误**：上条 gctx 冒烟的 fixture 把数据直接写在 `/0/DATA/0`
  dataset 上（照 v2 的 bug 造的假官方结构），未按官方 group 布局；
  真官方结构 v2 必 AttributeError（见 tests/ 复现用例）。数值结论
  （gct 侧与 wcs 口径）不受影响，gctx 侧以 v3 官方结构 fixture 重验。

## 测试（2026-10-05，离线 venv 实跑，tests/test_cmap_connectivity.py）

- 官方结构 GCTx fixture（/0/META/ROW|COL/id + /0/DATA/0/matrix，23 基因
  × 7 谱**非方阵**，存储方向 (7,23)）：v2 路径在其上 `.shape`
  AttributeError 复现锁定；v3 读出转置后与同数据 GCT 文本逐分一致。
- 已知答案：sig0=查询同向量 → ρ=+1、direction=same_direction；
  sig1=精确负向量 → ρ=−1、direction=reversed；与 scipy
  spearmanr/pearsonr 全谱对拍。
- 方向/排序：spearman/pearson 升序（最强逆转 sig1 居首）、wcs 降序；
  v2 的统一降序在相关类上会把同向最强排最前（测试锁定不再发生）。
- wcs 手算金标准：a_up/a_down 位次缩放均值与 score 逐位对拍；
  log 含 "NOT the official CMap WTCS/NCS/tau" 诚实标注。
- 契约锁定：输出列 == manifest "输出：" 行声明（tomllib）；gz 两形态
  （gct.gz / gctx.gz）嗅探与分块（chunk=3）不分块一致性。
