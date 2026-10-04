# 中文注释：LINCS/CMap 连通性评分。
# 输入构念：查询签名（gene,value）+ 参考谱矩阵。参考矩阵按魔数自动识别
# 三种物理形态：GCT 文本（#1.2 头）、GCTx（HDF5）、以及两者各自的 .gz
# 压缩形态（gzipped GCTx 先解压到 WORKDIR 再读，磁盘成本入 log）。
# 规则来源：方案 §18.2——方向相容比较；输出完整排名，条件筛选（细胞类型/
# 剂量/时间）在下游注册表完成，不在此挑最强逆转。
# 分块方向：按签名【列】分块（完整基因集 × 列块）——保证每条参考谱的
# 得分/位次建立在全部共同基因上。按基因行分块会把同一签名拆成多个
# 子集相关的行（历史缺陷，v2 修正）。
# v3 GCTx 口径修正（审计 F14，对齐 cmapPy parse_gctx 官方语义）：
# - 矩阵在 /0/DATA/0/matrix（/0/DATA/0 是 group 不是 dataset——v2 直接
#   读 group，真实官方结构上 .shape 直接 AttributeError）。
# - 官方存储方向是 (n_signature, n_gene)：COL（签名）在轴 0、ROW（基因）
#   在轴 1，cmapPy 读出后 transpose。v2 假设存储即 gene×signature，方向
#   相反。本版按官方方向读入并转置成 gene×signature，非方阵校验锁定。
# 方向语义（v3 修正）：
# - spearman/pearson：正=同向、负=逆转；逆转候选按【升序】取最前。
# - wcs（探索性位次逆转分，见下）：正=逆转方向；按【降序】取最前。
# - v2 对所有 metric 统一降序是错的（降序≠统一逆转）；本版 direction 列
#   逐行标注 same_direction/reversed，排序随 metric 语义走。
# wcs 诚实命名（v3）：输出列/日志/文档统一称 exploratory rank-reversal
# score（探索性位次逆转分）——up/down 平均相对位次之差，【非】官方
# WTCS/NCS/tau，不可与 CMap 官方分数混用。metric 标识 "wcs" 仅为向后
# 兼容保留。
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

# 逆转方向的排序：相关类最负=最强逆转（升序在前）；rank-reversal 分最正=
# 最强逆转（降序在前）。v2 的统一降序已废止。
REVERSAL_ASCENDING = {"spearman": True, "pearson": True, "wcs": False}

METRIC_LABEL = {
    "spearman": "spearman rho (positive=same direction, negative=reversal)",
    "pearson": "pearson r (positive=same direction, negative=reversal)",
    "wcs": "wcs — exploratory rank-reversal score (a_down - a_up); "
           "NOT the official CMap WTCS/NCS/tau",
}


def fail(message):
    print(message, file=sys.stderr)
    sys.exit(2)


def direction_label(metric, score):
    """逐行方向标注。相关类：符号=同向/逆转；wcs：正=逆转方向。

    不可估计 → "not_estimable"（不用 "NA"：pandas 默认 NA 值表会把
    "NA" 读回成 NaN，下游拿不到可区分的标签）。"""
    if score is None or np.isnan(score):
        return "not_estimable"
    if metric in ("spearman", "pearson"):
        if score > 0:
            return "same_direction"
        if score < 0:
            return "reversed"
        return "flat"
    # wcs：score = a_down - a_up，正 = 查询上调基因在参考谱位次垫底 = 逆转。
    if score > 0:
        return "reversed"
    if score < 0:
        return "concordant"
    return "flat"


def correlation_scores(metric, query_vector, data):
    """全基因相关：spearman 先秩变换；零方差分母 → NaN。"""
    if metric == "spearman":
        qc = rankdata(query_vector)
        qc = qc - qc.mean()
        ranks = rankdata(data, axis=0).astype(np.float64)
        ranks = ranks - ranks.mean(axis=0, keepdims=True)
        left = qc
    else:
        left = query_vector - query_vector.mean()
        ranks = data - data.mean(axis=0, keepdims=True)
    denom = np.sqrt(float((left ** 2).sum()) * (ranks ** 2).sum(axis=0))
    return (left[:, None] * ranks).sum(axis=0) / np.where(denom == 0, np.nan, denom)


def wcs_scores(query_vector, common_list, data, up_n, down_n):
    """探索性位次逆转分：a_down − a_up（up/down 集按查询带符号值取）。

    返回 (scores, a_up, a_down, n_overlap)。位次缩放 [0,1]，1=参考谱中
    最强。v2 修正沿用：up = 值最正前 up_n、down = 值最负后 down_n。
    """
    order = np.argsort(-query_vector)
    up_set = set(common_list[order[: min(up_n, len(order))]])
    down_set = set(common_list[order[len(order) - min(down_n, len(order)):]]) if down_n else set()
    ranks_desc = data.shape[0] - 1 - np.argsort(np.argsort(-data, axis=0), axis=0)
    scaled = ranks_desc / max(data.shape[0] - 1, 1)
    up_mask = np.isin(common_list, list(up_set))
    down_mask = np.isin(common_list, list(down_set))
    a_up = scaled[up_mask].mean(axis=0) if up_mask.any() else np.zeros(data.shape[1])
    a_down = scaled[down_mask].mean(axis=0) if down_mask.any() else np.zeros(data.shape[1])
    return a_down - a_up, a_up, a_down, len(up_set & down_set)


def sniff_reference(path, workdir):
    """返回 (kind, source, gz, path)。gzipped GCTx 解压到 WORKDIR 后改路径。"""
    with open(path, "rb") as handle:
        head = handle.read(8)
    if head[:2] == b"\x1f\x8b":
        with gzip.open(path, "rb") as handle:
            inner = handle.read(8)
        kind = "gct" if inner[:3] == b"#1." else "gctx"
        source = f"{kind}.gz"
        if kind == "gctx":
            # h5py 需要真实文件：流式解压到 WORKDIR（LINCS Level5 解压后 ~20GB 级）。
            target = os.path.join(workdir, "reference.gctx")
            with gzip.open(path, "rb") as src, open(target, "wb") as dst:
                shutil.copyfileobj(src, dst, 1 << 20)
            path = target
            source = "gctx.gz→解压到WORKDIR"
        return kind, source, True, path
    if head[:3] == b"#1.":
        return "gct", "gct", False, path
    if head[:3] == b"\x89HD":
        return "gctx", "gctx", False, path
    fail(f"unrecognized reference format (magic={head[:4]!r}); expected GCT (#1.2) / GCTx (HDF5), plain or gzip")


def read_gct_header(path, gz):
    """GCT 头：#1.2 / nrow ncol / NAME [Description] sample...；数据列 = 表头末 n_cols 列。"""
    opener = gzip.open if gz else open
    with opener(path, "rt", encoding="utf-8") as handle:
        first = handle.readline()
        if not first.startswith("#1.2"):
            fail("GCT does not start with #1.2 header line")
        dims = handle.readline().split()
        n_rows, n_cols = int(dims[0]), int(dims[1])
        names = handle.readline().rstrip("\n").split("\t")
    if len(names) - 1 < n_cols:
        fail(f"GCT header has {len(names) - 1} data columns, declared n_cols={n_cols}")
    return names[0], names[-n_cols:], n_rows, n_cols


def open_gctx(path):
    """官方 GCTx：/0/DATA/0/matrix（cmapPy parse_gctx 语义），存储方向
    (n_signature, n_gene) → 读出转置 gene×signature。"""
    import h5py

    handle = h5py.File(path, "r")
    if "/0/DATA/0/matrix" not in handle:
        if "/0/DATA/0" in handle and isinstance(handle["/0/DATA/0"], h5py.Dataset):
            fail("gctx has a dataset at /0/DATA/0 (legacy non-official layout); "
                 "official GCTx stores the matrix at /0/DATA/0/matrix")
        fail("gctx missing /0/DATA/0/matrix (official cmapPy layout)")

    def decode(array):
        return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in array[...]]

    rows = decode(handle["/0/META/ROW/id"])   # 基因（GCT 行）
    cols = decode(handle["/0/META/COL/id"])   # 签名（GCT 列）
    dataset = handle["/0/DATA/0/matrix"]
    expected = (len(cols), len(rows))  # 官方存储：轴0=COL(签名)、轴1=ROW(基因)
    if dataset.shape != expected:
        fail(
            f"gctx /0/DATA/0/matrix shape {dataset.shape} != official storage "
            f"orientation (n_signatures, n_genes)={expected}; cmapPy writes COL "
            "on axis 0 / ROW on axis 1 and transposes on load (parse_gctx)"
        )
    return handle, dataset, rows, cols


def main():
    query_path = os.environ["AUTONOMICS_INPUT0"]
    reference_path = os.environ["AUTONOMICS_INPUT1"]
    out_tsv = os.environ["AUTONOMICS_OUTPUT0"]
    out_log = os.environ["AUTONOMICS_OUTPUT1"]

    if METRIC not in ("spearman", "pearson", "wcs"):
        fail(f"metric must be spearman|pearson|wcs, got {METRIC}")

    query = pd.read_csv(query_path, sep="\t")
    for column in (GENE_COL, VALUE_COL):
        if column not in query.columns:
            fail(f"query TSV missing column: {column}")
    query = query[[GENE_COL, VALUE_COL]].dropna()
    query = query.groupby(GENE_COL, as_index=False)[VALUE_COL].mean()
    query_indexed = query.set_index(GENE_COL)[VALUE_COL]
    query_indexed.index = query_indexed.index.astype(str)

    kind, source, gz, reference_path = sniff_reference(reference_path, WORKDIR)

    handle = None
    if kind == "gct":
        id_name, data_names, n_rows, n_cols = read_gct_header(reference_path, gz)
        compression = "gzip" if gz else None

        def block_of(j0, j1):
            """全基因 × 列块：按行流读、只保留目标列，拼回完整基因集。"""
            wanted = data_names[j0:j1]
            keep = {id_name} | set(wanted)
            parts = []
            for chunk in pd.read_csv(
                reference_path, sep="\t", skiprows=2, index_col=0, chunksize=1000,
                compression=compression, usecols=lambda c: c in keep,
            ):
                parts.append(chunk)
            frame = pd.concat(parts)
            return frame[wanted]
    else:
        handle, dataset, rows, cols = open_gctx(reference_path)
        n_rows, n_cols = len(rows), len(cols)

        def block_of(j0, j1):
            # 存储方向 (签名, 基因)：签名块切轴 0，读出转置成 gene×signature。
            matrix = np.asarray(dataset[j0:j1, :], dtype=np.float64).T
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
        common_list = np.array(common.astype(str).to_list())

        if METRIC in ("spearman", "pearson"):
            scores = correlation_scores(METRIC, q, data)
            a_up = np.full(data.shape[1], np.nan)
            a_down = np.full(data.shape[1], np.nan)
        else:  # wcs — exploratory rank-reversal score
            scores, a_up, a_down, n_overlap = wcs_scores(q, common_list, data, UP_N, DOWN_N)
            if n_overlap and not warned_overlap:
                warned_overlap = True
                print(
                    f"WARN wcs up/down sets overlap ({n_overlap} genes); "
                    "increase the reference gene overlap or decrease up_n/down_n",
                    file=sys.stderr,
                )

        directions = [direction_label(METRIC, float(value)) for value in scores]
        results.append(pd.DataFrame({
            "signature_id": block.columns.astype(str),
            "metric": METRIC,
            "score": scores,
            "direction": directions,
            "n_common": len(common),
            "a_up": a_up,
            "a_down": a_down,
        }))
        n_total += j1 - j0
        print(f"block {j0}-{j1}/{n_cols} done", file=sys.stderr)

    if handle is not None:
        handle.close()

    # 排序随 metric 逆转语义：相关类升序（最负=最强逆转在前），
    # rank-reversal 分降序（最正=最强逆转在前）。完整排名照旧输出。
    ascending = REVERSAL_ASCENDING[METRIC]
    out = pd.concat(results, ignore_index=True).sort_values("score", ascending=ascending)
    out.to_csv(out_tsv, sep="\t", index=False)
    with open(out_log, "w", encoding="utf-8") as handle:
        handle.write(
            f"metric={METRIC}\tlabel: {METRIC_LABEL[METRIC]}\n"
            f"sort: reversal-first ordering ({'ascending' if ascending else 'descending'} "
            f"for {METRIC}); direction column labels each signature "
            "same_direction/reversed per metric convention\n"
            f"query_genes={len(query)}\treference_signatures={n_total}\n"
            f"source={source}\trows={n_rows}\tcols={n_cols}\tcommon_genes={n_common_seen}\n"
            "gctx: /0/DATA/0/matrix, official (n_signature, n_gene) storage, "
            "transposed to gene x signature on load (cmapPy parse_gctx semantics)\n"
            f"rule: full ranking, no condition filtering here (plan §18.2)\n"
        )


if __name__ == "__main__":
    main()
