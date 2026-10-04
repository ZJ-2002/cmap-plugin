# 离线测试（审计 F14）：官方结构 GCTx fixture（/0/DATA/0/matrix，存储方向
# (n_signature, n_gene)）+ 已知答案 + 方向/排序语义 + 契约锁定。不打网络。
import gzip
import importlib.util
import os
import subprocess
import sys
import tomllib
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
from scipy.stats import pearsonr, spearmanr

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "cmap_connectivity.py"
MANIFEST = REPO / "manifest.toml"
PYTHON = sys.executable

# 23 基因 × 7 谱（非方阵——方向装错的口子会被 7≠23 卡死）。
N_GENES, N_SIGS = 23, 7
RNG_SEED = 7
rng = np.random.default_rng(RNG_SEED)
QUERY = np.round(rng.normal(size=N_GENES), 4)
GENES = [f"GENE{i:02d}" for i in range(N_GENES)]
SIGS = [f"sig{i}" for i in range(N_SIGS)]


def build_matrix():
    """sig0=查询同向量、sig1=精确负向量、sig2=查询×2、sig3=2×查询+噪声、
    sig4=常向量（零方差→NaN）、sig5/sig6=独立噪声。"""
    noise = rng.normal(size=N_GENES)
    matrix = np.column_stack([
        QUERY.copy(),
        -QUERY.copy(),
        2.0 * QUERY,
        2.0 * QUERY + 0.05 * noise,
        np.full(N_GENES, 3.25),
        rng.normal(size=N_GENES),
        rng.normal(size=N_GENES),
    ])
    return np.round(matrix, 4)


MATRIX = build_matrix()


def write_query(path):
    pd.DataFrame({"gene": GENES, "value": QUERY}).to_csv(path, sep="\t", index=False)


def write_gct(path):
    lines = ["#1.2", f"{N_GENES}\t{N_SIGS}", "NAME\tDescription\t" + "\t".join(SIGS)]
    for gene, row in zip(GENES, MATRIX):
        lines.append(gene + "\tna\t" + "\t".join(f"{v:.4f}" for v in row))
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_gctx(path, official=True):
    """官方 cmapPy 布局：/0/META/ROW|COL/id + /0/DATA/0/matrix，
    存储方向 (n_signature, n_gene)——COL 轴 0、ROW 轴 1（读出需转置）。"""
    with h5py.File(path, "w") as handle:
        row = handle.create_group("/0/META/ROW")
        col = handle.create_group("/0/META/COL")
        row.create_dataset("id", data=np.array(GENES, dtype=h5py.string_dtype("utf-8")))
        col.create_dataset("id", data=np.array(SIGS, dtype=h5py.string_dtype("utf-8")))
        data = handle.create_group("/0/DATA/0")
        data.create_dataset("matrix", data=MATRIX.T, dtype="f8")  # (7, 23)


def run_node(tmp_path, reference, extra_env=None, expect_ok=True):
    tmp_path.mkdir(parents=True, exist_ok=True)
    query = tmp_path / "query.tsv"
    write_query(query)
    out_tsv = tmp_path / "cmap_scores.tsv"
    out_log = tmp_path / "cmap.log"
    env = dict(os.environ)
    env.update({
        "AUTONOMICS_INPUT0": str(query),
        "AUTONOMICS_INPUT1": str(reference),
        "AUTONOMICS_OUTPUT0": str(out_tsv),
        "AUTONOMICS_OUTPUT1": str(out_log),
        "AUTONOMICS_WORKDIR": str(tmp_path),
    })
    env.update(extra_env or {})
    result = subprocess.run([PYTHON, str(SCRIPT)], env=env,
                            capture_output=True, text=True)
    if expect_ok:
        assert result.returncode == 0, result.stderr
    table = pd.read_csv(out_tsv, sep="\t") if result.returncode == 0 else None
    log = out_log.read_text(encoding="utf-8") if out_log.exists() else ""
    return table, log, result


def load_module():
    spec = importlib.util.spec_from_file_location("cmap_connectivity", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules["cmap_connectivity"] = module
    spec.loader.exec_module(module)
    return module


# --- 官方结构：v2 在其上必 AttributeError（审计缺陷复现锁定）。---

def test_v2_dataset_path_attributeerror_on_official_layout(tmp_path):
    gctx = tmp_path / "ref.gctx"
    write_gctx(gctx)
    with h5py.File(gctx, "r") as handle:
        node = handle["/0/DATA/0"]  # v2 直接把这个 group 当 dataset
        with pytest.raises(AttributeError):
            node.shape  # noqa: B018 — 正是 v2 在官方结构上的崩溃点
    # v3 读官方布局无碍，且转置方向正确（见下）。
    mod = load_module()
    handle, dataset, rows, cols = mod.open_gctx(str(gctx))
    assert dataset.shape == (N_SIGS, N_GENES)
    assert rows == GENES and cols == SIGS
    handle.close()


def test_v3_rejects_legacy_dataset_at_data0(tmp_path):
    # 旧的非官方布局（数据直接挂在 /0/DATA/0 dataset 上，v2 冒烟 fixture 形态）
    # → 明确退出 2，指出官方路径，不再静默错读。
    legacy = tmp_path / "legacy.gctx"
    with h5py.File(legacy, "w") as handle:
        handle.create_dataset("/0/DATA/0", data=MATRIX)
        handle.create_dataset("/0/META/ROW/id",
                              data=np.array(GENES, dtype=h5py.string_dtype("utf-8")))
        handle.create_dataset("/0/META/COL/id",
                              data=np.array(SIGS, dtype=h5py.string_dtype("utf-8")))
    _, _, result = run_node(tmp_path / "run", legacy, expect_ok=False)
    assert result.returncode == 2
    assert "/0/DATA/0/matrix" in result.stderr


def test_gctx_orientation_mismatch_exits_2(tmp_path):
    # 存储方向写成 (gene, signature)（v2 的错误假设）→ 官方校验退出 2。
    wrong = tmp_path / "wrong.gctx"
    with h5py.File(wrong, "w") as handle:
        handle.create_dataset("/0/META/ROW/id",
                              data=np.array(GENES, dtype=h5py.string_dtype("utf-8")))
        handle.create_dataset("/0/META/COL/id",
                              data=np.array(SIGS, dtype=h5py.string_dtype("utf-8")))
        handle.create_dataset("/0/DATA/0/matrix", data=MATRIX)  # (23, 7) ≠ (7, 23)
    _, _, result = run_node(tmp_path / "run", wrong, expect_ok=False)
    assert result.returncode == 2
    assert "(n_signatures, n_genes)" in result.stderr


# --- 已知答案：gct 与官方结构 gctx 逐分一致，且对拍 scipy。---

@pytest.mark.parametrize("metric,scipy_fn", [("spearman", spearmanr), ("pearson", pearsonr)])
def test_known_answers_gct_and_gctx_agree(tmp_path, metric, scipy_fn):
    gct = tmp_path / "ref.gct"
    gctx = tmp_path / "ref.gctx"
    write_gct(gct)
    write_gctx(gctx)
    env = {"CMAP_METRIC": metric}
    table_gct, _, _ = run_node(tmp_path / "gct", gct, env)
    table_gctx, log, _ = run_node(tmp_path / "gctx", gctx, env)
    # 两种物理形态读的是同一矩阵 → 分数逐位一致（转置方向正确的强证据）。
    pd.testing.assert_frame_equal(table_gct, table_gctx)
    assert "source=gctx" in log
    # 已知答案：sig0=+1 同向；sig1=−1 逆转；sig4 零方差 → NaN。
    scores = dict(zip(table_gct["signature_id"], table_gct["score"]))
    assert abs(scores["sig0"] - 1.0) < 1e-12
    assert abs(scores["sig1"] + 1.0) < 1e-12
    assert np.isnan(scores["sig4"])
    # 全谱对拍 scipy（NaN 谱除外）。
    for j, sig in enumerate(SIGS):
        expected = scipy_fn(QUERY, MATRIX[:, j]).statistic
        got = scores[sig]
        if np.isnan(expected):
            assert np.isnan(got)
        else:
            assert abs(got - expected) < 1e-12


def test_direction_labels_and_metric_ascending_sort(tmp_path):
    gctx = tmp_path / "ref.gctx"
    write_gctx(gctx)
    table, log, _ = run_node(tmp_path / "run", gctx, {"CMAP_METRIC": "spearman"})
    by_id = table.set_index("signature_id")
    assert by_id.loc["sig0", "direction"] == "same_direction"
    assert by_id.loc["sig1", "direction"] == "reversed"
    assert by_id.loc["sig4", "direction"] == "not_estimable"  # NaN 分（不用 "NA"，pandas 会读回 NaN）
    # 相关类按升序：最强逆转（sig1=−1）在最前，v2 统一降序会把 sig0 排最前；
    # NaN 分（sig4 零方差）pandas 默认垫尾。
    assert table.iloc[0]["signature_id"] == "sig1"
    assert table.iloc[-2]["signature_id"] == "sig0"
    assert table.iloc[-1]["signature_id"] == "sig4"
    assert (table["metric"] == "spearman").all()
    assert "ascending" in log


# --- wcs：诚实命名 + 手算金标准 + 降序方向。---

def test_wcs_hand_computed_golden_and_descending(tmp_path):
    up_n = down_n = 5
    gctx = tmp_path / "ref.gctx"
    write_gctx(gctx)
    table, log, _ = run_node(tmp_path / "run", gctx, {
        "CMAP_METRIC": "wcs", "CMAP_UP_N": str(up_n), "CMAP_DOWN_N": str(down_n),
    })
    # 手算：up=值最正前 5、down=最负后 5；位次缩放 [0,1]（1=参考谱最强）。
    order = np.argsort(-QUERY)
    up_set = set(np.array(GENES)[order[:up_n]])
    down_set = set(np.array(GENES)[order[len(order) - down_n:]])
    by_id = table.set_index("signature_id")
    for j, sig in enumerate(SIGS):
        column = MATRIX[:, j]
        ranks_desc = N_GENES - 1 - np.argsort(np.argsort(-column))
        scaled = ranks_desc / (N_GENES - 1)
        up_mean = scaled[[g in up_set for g in GENES]].mean()
        down_mean = scaled[[g in down_set for g in GENES]].mean()
        expected = down_mean - up_mean
        assert abs(by_id.loc[sig, "score"] - expected) < 1e-12
        assert abs(by_id.loc[sig, "a_up"] - up_mean) < 1e-12
        assert abs(by_id.loc[sig, "a_down"] - down_mean) < 1e-12
    # 完全逆转谱 sig1：查询 up 基因在 −query 谱中位次垫底 → 正分居首（降序）。
    assert by_id.loc["sig1", "direction"] == "reversed"
    assert table.iloc[0]["signature_id"] == "sig1"
    assert by_id.loc["sig1", "score"] > 0 and by_id.loc["sig0", "score"] < 0
    # 诚实命名：log 必须声明非官方 WTCS/NCS/tau。
    assert "NOT the official CMap WTCS/NCS/tau" in log
    assert "exploratory rank-reversal" in log


# --- gz 形态、分块一致性、契约锁定。---

def test_gz_forms_sniff_and_chunk_invariance(tmp_path):
    gct_gz = tmp_path / "ref.gct.gz"
    raw = tmp_path / "plain.gct"
    write_gct(raw)
    with open(raw, "rb") as src, gzip.open(gct_gz, "wb") as dst:
        dst.writelines(src)
    gctx = tmp_path / "ref.gctx"
    write_gctx(gctx)
    gctx_gz = tmp_path / "ref.gctx.gz"
    with open(gctx, "rb") as src, gzip.open(gctx_gz, "wb") as dst:
        dst.writelines(src)

    full, _, _ = run_node(tmp_path / "full", gctx, {"CMAP_CHUNK": "2000"})
    chunked, _, _ = run_node(tmp_path / "chunked", gctx, {"CMAP_CHUNK": "3"})
    from_gct_gz, _, _ = run_node(tmp_path / "gctgz", gct_gz)
    from_gctx_gz, log, _ = run_node(tmp_path / "gctxgz", gctx_gz)
    pd.testing.assert_frame_equal(full, chunked)      # 分块不改变分数
    pd.testing.assert_frame_equal(full, from_gct_gz)  # gct.gz 同数据同分
    pd.testing.assert_frame_equal(full, from_gctx_gz)  # gctx.gz 解压后同分
    assert "gctx.gz" in log


def test_output_columns_match_manifest_declaration(tmp_path):
    gctx = tmp_path / "ref.gctx"
    write_gctx(gctx)
    table, _, _ = run_node(tmp_path / "run", gctx)
    doc = tomllib.loads(MANIFEST.read_text(encoding="utf-8"))["nodes"][0]["doc"]
    declared = doc.split("输出：cmap_scores.tsv ——", 1)[1].split("+ cmap.log", 1)[0]
    declared_columns = [c.strip() for c in declared.strip().split("/")]
    assert declared_columns == list(table.columns)
    assert "direction" in declared_columns and "metric" in declared_columns


def test_import_is_side_effect_free():
    load_module()  # main 守卫：import 不触发 env 读取/退出


def test_direction_label_unit():
    mod = load_module()
    assert mod.direction_label("spearman", 0.5) == "same_direction"
    assert mod.direction_label("spearman", -0.5) == "reversed"
    assert mod.direction_label("pearson", 0.0) == "flat"
    assert mod.direction_label("pearson", float("nan")) == "not_estimable"
    assert mod.direction_label("wcs", 0.3) == "reversed"
    assert mod.direction_label("wcs", -0.3) == "concordant"
    assert mod.direction_label("wcs", 0.0) == "flat"
