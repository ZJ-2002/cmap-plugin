# 中文注释：LINCS/CMap 连通性评分。
# 输入构念：查询签名（gene,value）+ GCT 参考谱矩阵。规则来源：方案 §18.2
# ——方向相容比较；输出完整排名，条件筛选（细胞类型/剂量/时间）在下游
# 注册表完成，不在此挑最强逆转。统计口径：spearman/pearson 全基因相关；
# wcs 按查询 |value| 排名取 up_n/down_n，参考谱位次缩放 [0,1]，
# score = a_down - a_up（高 = 逆转）。失败处理：GCT 头/共同基因不足 → 退出 2。
# 解释边界：表达逆转≠药效；非肌肉细胞系为跨模型外推，需降级（§18.2）。
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata

METRIC = os.environ.get("CMAP_METRIC", "spearman").lower()
GENE_COL = os.environ.get("CMAP_GENE_COL", "gene")
VALUE_COL = os.environ.get("CMAP_VALUE_COL", "value")
UP_N = int(os.environ.get("CMAP_UP_N", "150"))
DOWN_N = int(os.environ.get("CMAP_DOWN_N", "150"))
CHUNK = int(os.environ.get("CMAP_CHUNK", "2000"))
QUERY = os.environ["AUTONOMICS_INPUT0"]
GCT = os.environ["AUTONOMICS_INPUT1"]
OUT_TSV = os.environ["AUTONOMICS_OUTPUT0"]
OUT_LOG = os.environ["AUTONOMICS_OUTPUT1"]


def fail(message):
    print(message, file=sys.stderr)
    sys.exit(2)


if METRIC not in ("spearman", "pearson", "wcs"):
    fail(f"metric must be spearman|pearson|wcs, got {METRIC}")

query = pd.read_csv(QUERY, sep="\t")
for column in (GENE_COL, VALUE_COL):
    if column not in query.columns:
        fail(f"query TSV missing column: {column}")
query = query[[GENE_COL, VALUE_COL]].dropna()
query = query.groupby(GENE_COL, as_index=False)[VALUE_COL].mean()
query_values = query[VALUE_COL].to_numpy(dtype=float)

with open(GCT, "r", encoding="utf-8") as handle:
    first = handle.readline()
    if not first.startswith("#1.2"):
        fail("GCT does not start with #1.2 header line")
    dims = handle.readline().split()
    n_rows, n_cols = int(dims[0]), int(dims[1])

reader = pd.read_csv(GCT, sep="\t", skiprows=2, index_col=0, chunksize=CHUNK)
results = []
n_total = 0
query_genes = query[GENE_COL].astype(str)

if METRIC == "spearman":
    query_rank = rankdata(query_values)
    query_centered = query_rank - query_rank.mean()
    query_ss = float(np.sum(query_centered ** 2))
elif METRIC == "pearson":
    query_centered = query_values - query_values.mean()
    query_ss = float(np.sum(query_centered ** 2))

for chunk in reader:
    chunk.index = chunk.index.astype(str)
    if "Description" in chunk.columns:
        chunk = chunk.drop(columns=["Description"])
    common = pd.Index(query_genes).intersection(chunk.index)
    if len(common) < 10:
        if n_total == 0:
            fail(f"fewer than 10 common genes between query and GCT ({len(common)})")
    q = query.set_index(GENE_COL).loc[common, VALUE_COL].to_numpy(dtype=float)
    block = chunk.loc[common].to_numpy(dtype=float)
    if METRIC == "spearman":
        qr = rankdata(q)
        qc = qr - qr.mean()
        ranks = np.apply_along_axis(rankdata, 0, block).astype(np.float64)
        ranks = ranks - ranks.mean(axis=0, keepdims=True)
        denom = np.sqrt(float((qc ** 2).sum()) * (ranks ** 2).sum(axis=0))
        scores = (qc[:, None] * ranks).sum(axis=0) / np.where(denom == 0, np.nan, denom)
        a_up = np.full(block.shape[1], np.nan)
        a_down = np.full(block.shape[1], np.nan)
    elif METRIC == "pearson":
        qc = q - q.mean()
        centered = block - block.mean(axis=0, keepdims=True)
        denom = np.sqrt((qc ** 2).sum() * (centered ** 2).sum(axis=0))
        scores = (qc[:, None] * centered).sum(axis=0) / np.where(denom == 0, np.nan, denom)
        a_up = np.full(block.shape[1], np.nan)
        a_down = np.full(block.shape[1], np.nan)
    else:  # wcs
        order = np.argsort(-np.abs(q))
        up_set = set(query_genes.iloc[order[:UP_N]])
        down_set = set(query_genes.iloc[order[-DOWN_N:]])
        # 参考谱按强度降序的相对位次，缩放到 [0,1]（1=最强）。
        ranks_desc = block.shape[0] - 1 - np.argsort(np.argsort(-block, axis=0), axis=0)
        scaled = ranks_desc / max(block.shape[0] - 1, 1)
        up_mask = np.isin(common.astype(str), list(up_set))
        down_mask = np.isin(common.astype(str), list(down_set))
        a_up = scaled[up_mask].mean(axis=0) if up_mask.any() else np.zeros(block.shape[1])
        a_down = scaled[down_mask].mean(axis=0) if down_mask.any() else np.zeros(block.shape[1])
        scores = a_down - a_up

    frame = pd.DataFrame({
        "signature_id": chunk.columns.astype(str),
        "score": scores,
        "n_common": len(common),
        "a_up": a_up,
        "a_down": a_down,
    })
    results.append(frame)
    n_total += len(frame)

out = pd.concat(results, ignore_index=True).sort_values("score", ascending=False)
out.to_csv(OUT_TSV, sep="\t", index=False)
with open(OUT_LOG, "w", encoding="utf-8") as handle:
    handle.write(
        f"metric={METRIC}\tquery_genes={len(query)}\treference_signatures={n_total}\n"
        f"rule: full ranking, no condition filtering here (plan §18.2)\n"
    )
