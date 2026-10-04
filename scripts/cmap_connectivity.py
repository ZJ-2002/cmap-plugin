# 中文注释：LINCS/CMap 连通性评分。
# 输入构念：查询签名（gene,value）+ 参考谱矩阵。参考矩阵按魔数自动识别
# 三种物理形态：GCT 文本（#1.2 头）、GCTx（HDF5，/0/DATA/0）、以及两者
# 各自的 .gz 压缩形态（gzipped GCTx 先解压到 WORKDIR 再读，磁盘成本入 log）。
# 规则来源：方案 §18.2——方向相容比较；输出完整排名，条件筛选（细胞类型/
# 剂量/时间）在下游注册表完成，不在此挑最强逆转。
# 分块方向：按签名【列】分块（完整基因集 × 列块）——保证每条参考谱的
# 得分/位次建立在全部共同基因上。按基因行分块会把同一签名拆成多个
# 子集相关的行（历史缺陷，v2 修正）。
# 统计口径：spearman/pearson 全基因相关；wcs 按查询带符号值取 up_n
# （最正）/down_n（最负）基因，参考谱位次缩放 [0,1]，
# score = a_down - a_up（高 = 逆转）。
# 失败处理：格式无法识别/GCT 头错/gctx 形状与元数据不符/共同基因 <10 → 退出 2。
# 解释边界：表达逆转≠药效；LINCS gctx 行标识符是 Entrez ID，符号查询需
# 先用同目录 gene_info 映射（geo_suppl 可取）；非肌肉细胞系为跨模型外推（§18.2）。
import gzip
import os
import shutil
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
WORKDIR = os.environ.get("AUTONOMICS_WORKDIR", "/work")
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
query_indexed = query.set_index(GENE_COL)[VALUE_COL]
query_indexed.index = query_indexed.index.astype(str)

# --- 参考矩阵物理形态嗅探：GCT / GCTx / 两者的 gzip 形态。---
gz = False
with open(GCT, "rb") as handle:
    head = handle.read(8)
if head[:2] == b"\x1f\x8b":
    gz = True
    with gzip.open(GCT, "rb") as handle:
        inner = handle.read(8)
    kind = "gct" if inner[:3] == b"#1." else "gctx"
    source = f"{kind}.gz"
    if kind == "gctx":
        # h5py 需要真实文件：流式解压到 WORKDIR（LINCS Level5 解压后 ~20GB 级）。
        target = os.path.join(WORKDIR, "reference.gctx")
        with gzip.open(GCT, "rb") as src, open(target, "wb") as dst:
            shutil.copyfileobj(src, dst, 1 << 20)
        GCT = target
        source = "gctx.gz→解压到WORKDIR"
elif head[:3] == b"#1.":
    kind, source = "gct", "gct"
elif head[:3] == b"\x89HD":
    kind, source = "gctx", "gctx"
else:
    fail(f"unrecognized reference format (magic={head[:4]!r}); expected GCT (#1.2) / GCTx (HDF5), plain or gzip")

if kind == "gct":
    # GCT 头：#1.2 / nrow ncol / NAME [Description] sample...；数据列 = 表头末 n_cols 列。
    opener = gzip.open if gz else open
    with opener(GCT, "rt", encoding="utf-8") as handle:
        first = handle.readline()
        if not first.startswith("#1.2"):
            fail("GCT does not start with #1.2 header line")
        dims = handle.readline().split()
        n_rows, n_cols = int(dims[0]), int(dims[1])
        names = handle.readline().rstrip("\n").split("\t")
    id_name = names[0]
    if len(names) - 1 < n_cols:
        fail(f"GCT header has {len(names) - 1} data columns, declared n_cols={n_cols}")
    data_names = names[-n_cols:]
    compression = "gzip" if gz else None

    def block_of(j0, j1):
        """全基因 × 列块：按行流读、只保留目标列，拼回完整基因集。"""
        wanted = data_names[j0:j1]
        keep = {id_name} | set(wanted)
        parts = []
        for chunk in pd.read_csv(
            GCT, sep="\t", skiprows=2, index_col=0, chunksize=1000,
            compression=compression, usecols=lambda c: c in keep,
        ):
            parts.append(chunk)
        frame = pd.concat(parts)
        return frame[wanted]
else:
    import h5py

    handle = h5py.File(GCT, "r")
    dataset = handle["/0/DATA/0"]

    def decode(array):
        return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in array[...]]

    rows = decode(handle["/0/META/ROW/id"])
    cols = decode(handle["/0/META/COL/id"])
    if dataset.shape != (len(rows), len(cols)):
        fail(f"gctx /0/DATA/0 shape {dataset.shape} != ({len(rows)} genes, {len(cols)} signatures)")
    n_rows, n_cols = len(rows), len(cols)

    def block_of(j0, j1):
        matrix = np.asarray(dataset[:, j0:j1], dtype=np.float64)
        return pd.DataFrame(matrix, index=rows, columns=cols[j0:j1])

results = []
n_total = 0
n_common_seen = None
warned_overlap = False
for j0 in range(0, n_cols, CHUNK):
    j1 = min(j0 + CHUNK, n_cols)
    block = block_of(j0, j1)
    block.index = block.index.astype(str)
    common = query_indexed.index.intersection(block.index)
    if n_common_seen is None:
        if len(common) < 10:
            fail(f"fewer than 10 common genes between query and reference rows ({len(common)}); "
                 "LINCS gctx rows are Entrez ids — map query symbols via gene_info first")
        n_common_seen = len(common)
    q = query_indexed.loc[common].to_numpy(dtype=float)
    data = block.loc[common].to_numpy(dtype=float)

    if METRIC == "spearman":
        qc = rankdata(q)
        qc = qc - qc.mean()
        ranks = rankdata(data, axis=0).astype(np.float64)
        ranks = ranks - ranks.mean(axis=0, keepdims=True)
        denom = np.sqrt(float((qc ** 2).sum()) * (ranks ** 2).sum(axis=0))
        scores = (qc[:, None] * ranks).sum(axis=0) / np.where(denom == 0, np.nan, denom)
        a_up = np.full(data.shape[1], np.nan)
        a_down = np.full(data.shape[1], np.nan)
    elif METRIC == "pearson":
        qc = q - q.mean()
        centered = data - data.mean(axis=0, keepdims=True)
        denom = np.sqrt((qc ** 2).sum() * (centered ** 2).sum(axis=0))
        scores = (qc[:, None] * centered).sum(axis=0) / np.where(denom == 0, np.nan, denom)
        a_up = np.full(data.shape[1], np.nan)
        a_down = np.full(data.shape[1], np.nan)
    else:  # wcs
        # up/down 集按带符号值取（CMap 惯例）：up = 值最正的前 up_n，
        # down = 值最负的后 down_n。v1 按 |value| 排（强幅度 vs 弱幅度）
        # 是错误口径，已修；up_n+down_n ≥ 共同基因数时两集重叠，WARN 记录。
        order = np.argsort(-q)
        common_list = np.array(common.astype(str).to_list())
        up_set = set(common_list[order[: min(UP_N, len(order))]])
        down_set = set(common_list[order[len(order) - min(DOWN_N, len(order)):]]) if DOWN_N else set()
        overlap = up_set & down_set
        if overlap and not warned_overlap:
            warned_overlap = True
            print(
                f"WARN wcs up/down sets overlap ({len(overlap)} genes); "
                "increase the reference gene overlap or decrease up_n/down_n",
                file=sys.stderr,
            )
        # 参考谱按强度降序的相对位次，缩放到 [0,1]（1=最强）。
        ranks_desc = data.shape[0] - 1 - np.argsort(np.argsort(-data, axis=0), axis=0)
        scaled = ranks_desc / max(data.shape[0] - 1, 1)
        up_mask = np.isin(common_list, list(up_set))
        down_mask = np.isin(common_list, list(down_set))
        a_up = scaled[up_mask].mean(axis=0) if up_mask.any() else np.zeros(data.shape[1])
        a_down = scaled[down_mask].mean(axis=0) if down_mask.any() else np.zeros(data.shape[1])
        scores = a_down - a_up

    results.append(pd.DataFrame({
        "signature_id": block.columns.astype(str),
        "score": scores,
        "n_common": len(common),
        "a_up": a_up,
        "a_down": a_down,
    }))
    n_total += j1 - j0
    print(f"block {j0}-{j1}/{n_cols} done", file=sys.stderr)

if kind == "gctx":
    handle.close()

out = pd.concat(results, ignore_index=True).sort_values("score", ascending=False)
out.to_csv(OUT_TSV, sep="\t", index=False)
with open(OUT_LOG, "w", encoding="utf-8") as handle:
    handle.write(
        f"metric={METRIC}\tquery_genes={len(query)}\treference_signatures={n_total}\n"
        f"source={source}\trows={n_rows}\tcols={n_cols}\tcommon_genes={n_common_seen}\n"
        f"rule: full ranking, no condition filtering here (plan §18.2)\n"
    )
