"""
Protein List (GO 解析用) 作成アプリ
==================================

Streamlit アプリ。ProteinList 用の Excel を読み込み、"T-half" シートに
"UniProtGO" シートから Gene Ontology 情報を転記して、GO 解析用の
Protein list を作成します。

処理内容
--------
1. "T-half" の Gene Name (A列) と "UniProtGO" の From (A列) を照合。
   同じ Gene Name に複数行がある場合の採用ルール:
     (1) Reviewed (C列) が "reviewed" の行を最優先
     (2) その中で Length (H列) が最大の行を採用
2. 採用行から下記を転記:
     C列 <- UniProtGO J列  Gene Ontology (cellular component)
     D列 <- UniProtGO I列  Gene Ontology (biological process)
     E列 <- UniProtGO L列  Gene Ontology (molecular function)
     F/G/H列 <- 採用行の Entry (B列)   (F=C用, G=D用, H=E用)
3. C/D/E 列のセミコロン区切り文字列を分解し、重複を除いた一覧を
   I/J/K 列に出力 (セミコロンは削除)。
4. 各文字列が該当列に出現する回数を L/M/N 列に出力
   (I->L, J->M, K->N)。
5. cellular component / biological process / molecular function ごとに、
   各 GO 用語を持つタンパク質の T1/2 (HalfLife) の箱ひげ図を作成 (対数軸)。
   用語は出現頻度でフィルタ (>=20 タンパク質 かつ 全タンパク質の <=90%)。
   PNG でダウンロード可能。
6. 補足解析 (2026-10-27 追加): フィルタを通過した各 GO 用語について
   - 注釈あり vs なしの T1/2 を両側 Mann–Whitney U 検定 (順位の正規近似、同順位補正) で比較し、
     Cliff's δ (正 = 長寿命) と Benjamini–Hochberg の q 値を算出
   - Median 上位/下位の用語どうしの重なり (Jaccard 係数 = 共通タンパク質数 / 和集合、
     overlap 係数 = 共通タンパク質数 / 小さい方の用語のタンパク質数)
   - 条件付き比較: ある用語 (例: myelin sheath) から別の用語 (例: mitochondrion) の
     注釈を持つタンパク質を除いて検定
   結果は Excel (GO_term_statistics.xlsx) でダウンロード可能。

必要ライブラリ:
    streamlit, pandas, openpyxl, matplotlib

実行方法:
    streamlit run protein_list_go_app.py
"""

import io
import math
import random
import re
import statistics
import textwrap
from collections import Counter

import matplotlib
matplotlib.use("Agg")  # サーバー環境向けの非対話バックエンド
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
import pandas as pd  # noqa: E402
import streamlit as st  # noqa: E402

# 箱ひげ図の配色 (添付の例に合わせた薄紫)
BOX_FILL = "#9a9ae6"
BOX_EDGE = "#5b5bd6"
BOX_MEDIAN = "#3b3ba8"

# Median 上位/下位の強調色 (上位=赤, 下位=青)
HIGH_FILL = "#e06666"   # 赤
HIGH_EDGE = "#cc0000"
LOW_FILL = "#5b9bd5"    # 青
LOW_EDGE = "#1155cc"

# 箱ひげ図を作成する 3 カテゴリ (ラベル, out_df の列名)
BOX_CATEGORIES = [
    ("cellular component", "GO cellular component"),
    ("biological process", "GO biological process"),
    ("molecular function", "GO molecular function"),
]

# ---------------------------------------------------------------------------
# UniProtGO シートの元カラム名
# ---------------------------------------------------------------------------
COL_FROM = "From"
COL_ENTRY = "Entry"
COL_REVIEWED = "Reviewed"
COL_LENGTH = "Length"
COL_BP = "Gene Ontology (biological process)"   # I列
COL_CC = "Gene Ontology (cellular component)"    # J列
COL_MF = "Gene Ontology (molecular function)"    # L列

# 出力する T-half シートのヘッダー (A..N)
OUT_HEADERS = [
    "Gene Name",                       # A
    "HalfLife",                        # B
    "GO cellular component",           # C  <- UniProtGO J
    "GO biological process",           # D  <- UniProtGO I
    "GO molecular function",           # E  <- UniProtGO L
    "Entry (cellular component)",      # F  <- 採用行 Entry (C用)
    "Entry (biological process)",      # G  <- 採用行 Entry (D用)
    "Entry (molecular function)",      # H  <- 採用行 Entry (E用)
    "Unique GO cellular component",    # I  <- C の重複なし一覧
    "Unique GO biological process",    # J  <- D の重複なし一覧
    "Unique GO molecular function",    # K  <- E の重複なし一覧
    "Count (cellular component)",      # L  <- I の出現回数
    "Count (biological process)",      # M  <- J の出現回数
    "Count (molecular function)",      # N  <- K の出現回数
]


# ---------------------------------------------------------------------------
# 処理ロジック
# ---------------------------------------------------------------------------
def _split_terms(cell):
    """GO セルを ';' で分割し、前後空白を除いた非空の語リストを返す。"""
    if cell is None or (isinstance(cell, float) and pd.isna(cell)):
        return []
    return [t.strip() for t in str(cell).split(";") if t.strip()]


def select_best_row(candidates):
    """1 遺伝子に対する UniProtGO 候補行から採用行を 1 つ選ぶ。
    (1) Reviewed == 'reviewed' を優先  (2) その中で Length 最大。"""
    is_reviewed = (
        candidates[COL_REVIEWED].astype(str).str.strip().str.lower() == "reviewed"
    )
    reviewed = candidates[is_reviewed]
    pool = reviewed if len(reviewed) else candidates
    lengths = pd.to_numeric(pool[COL_LENGTH], errors="coerce").fillna(-1)
    return pool.loc[lengths.idxmax()]


def _val(row, col):
    """欠損は空文字に、それ以外は元の値を返す。"""
    if col not in row or pd.isna(row[col]):
        return ""
    return row[col]


def process(th_df, go_df):
    """完成した T-half データフレーム (A..N) と、遺伝子単位の中間表を返す。"""
    # 必須カラムの存在チェック
    missing = [c for c in (COL_FROM, COL_ENTRY, COL_REVIEWED, COL_LENGTH,
                           COL_BP, COL_CC, COL_MF) if c not in go_df.columns]
    if missing:
        raise ValueError(
            "UniProtGO シートに次の列が見つかりません: " + ", ".join(missing)
        )
    if th_df.shape[1] < 2:
        raise ValueError("T-half シートには Gene Name と HalfLife の 2 列が必要です。")

    gene_col = th_df.columns[0]      # A: Gene Name
    half_col = th_df.columns[1]      # B: HalfLife

    # From ごとに候補行をまとめておく
    go_by_gene = {k: v for k, v in go_df.groupby(go_df[COL_FROM].astype(str))}

    rows = []
    n_matched = 0
    for _, r in th_df.iterrows():
        gene = r[gene_col]
        cc = bp = mf = entry = ""
        cand = go_by_gene.get(str(gene))
        if cand is not None and len(cand):
            best = select_best_row(cand)
            cc = _val(best, COL_CC)
            bp = _val(best, COL_BP)
            mf = _val(best, COL_MF)
            entry = _val(best, COL_ENTRY)
            n_matched += 1
        # 採用行は 1 行なので、その Entry が C/D/E 共通の出所 -> F=G=H
        rows.append({
            "Gene Name": gene,
            "HalfLife": r[half_col],
            "GO cellular component": cc,
            "GO biological process": bp,
            "GO molecular function": mf,
            "Entry (cellular component)": entry,
            "Entry (biological process)": entry,
            "Entry (molecular function)": entry,
        })
    per_gene = pd.DataFrame(rows)

    # 重複なし一覧 + 出現回数 (列ごと)
    def unique_counts(series):
        counter = Counter()
        for cell in series:
            for term in _split_terms(cell):
                counter[term] += 1
        ordered = list(counter.keys())          # 初出順を保持
        return ordered, [counter[t] for t in ordered]

    cc_terms, cc_counts = unique_counts(per_gene["GO cellular component"])
    bp_terms, bp_counts = unique_counts(per_gene["GO biological process"])
    mf_terms, mf_counts = unique_counts(per_gene["GO molecular function"])

    n = max(len(per_gene), len(cc_terms), len(bp_terms), len(mf_terms))

    def pad(lst):
        lst = list(lst)
        return lst + [""] * (n - len(lst))

    out = pd.DataFrame({
        "Gene Name": pad(per_gene["Gene Name"]),
        "HalfLife": pad(per_gene["HalfLife"]),
        "GO cellular component": pad(per_gene["GO cellular component"]),
        "GO biological process": pad(per_gene["GO biological process"]),
        "GO molecular function": pad(per_gene["GO molecular function"]),
        "Entry (cellular component)": pad(per_gene["Entry (cellular component)"]),
        "Entry (biological process)": pad(per_gene["Entry (biological process)"]),
        "Entry (molecular function)": pad(per_gene["Entry (molecular function)"]),
        "Unique GO cellular component": pad(cc_terms),
        "Unique GO biological process": pad(bp_terms),
        "Unique GO molecular function": pad(mf_terms),
        "Count (cellular component)": pad(cc_counts),
        "Count (biological process)": pad(bp_counts),
        "Count (molecular function)": pad(mf_counts),
    })[OUT_HEADERS]

    stats = {
        "total_genes": len(per_gene),
        "matched_genes": n_matched,
        "unique_cc": len(cc_terms),
        "unique_bp": len(bp_terms),
        "unique_mf": len(mf_terms),
    }
    return out, stats


def to_excel_bytes(out_df, go_df):
    """完成 T-half と元の UniProtGO を含む Excel をバイト列で返す。"""
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        out_df.to_excel(writer, sheet_name="T-half", index=False)
        go_df.to_excel(writer, sheet_name="UniProtGO", index=False)
    buffer.seek(0)
    return buffer.getvalue()


# ---------------------------------------------------------------------------
# 箱ひげ図
# ---------------------------------------------------------------------------
_JP_FONT_CHECKED = False
_JP_FONT = None


def _jp_fontprops():
    """日本語対応フォントがあれば FontProperties を、無ければ None を返す。
    (GO 用語ラベルは英語なので既定フォントで問題なく描画される。)"""
    global _JP_FONT_CHECKED, _JP_FONT
    if _JP_FONT_CHECKED:
        return _JP_FONT
    _JP_FONT_CHECKED = True
    try:
        import matplotlib.font_manager as fm
        available = {f.name for f in fm.fontManager.ttflist}
        for name in ("Noto Sans CJK JP", "Noto Sans JP", "IPAexGothic",
                     "IPAGothic", "Yu Gothic", "Hiragino Sans", "MS Gothic",
                     "TakaoGothic", "VL Gothic", "Meiryo"):
            if name in available:
                _JP_FONT = fm.FontProperties(family=name)
                break
    except Exception:  # noqa: BLE001
        _JP_FONT = None
    return _JP_FONT


def collect_term_halflives(per_gene, cat_col, half_col="HalfLife"):
    """GO 用語 -> その用語を持つタンパク質の T1/2 値リスト、を返す。
    T1/2 が欠損 / 0 以下 の行は除外 (対数軸のため)。"""
    term_vals = {}
    hl = pd.to_numeric(per_gene[half_col], errors="coerce")
    for cell, v in zip(per_gene[cat_col], hl):
        if pd.isna(v) or v <= 0:
            continue
        for term in _split_terms(cell):
            term_vals.setdefault(term, []).append(float(v))
    return term_vals


def filter_terms(term_vals, total_proteins, min_count=20, max_frac=0.90,
                 top_n=0):
    """出現頻度でフィルタし、頻度の高い順に並べた (用語, 値リスト) を返す。
    条件: min_count <= (その用語を持つタンパク質数) <= max_frac * 全タンパク質。"""
    hi = max_frac * total_proteins
    items = [(t, v) for t, v in term_vals.items() if min_count <= len(v) <= hi]
    items.sort(key=lambda x: len(x[1]), reverse=True)
    if top_n and top_n > 0:
        items = items[:top_n]
    return items


# ---------------------------------------------------------------------------
# 補足解析: 用語ごとの検定 (Mann–Whitney + Cliff's δ + BH-FDR) と重なり (Jaccard)
# ---------------------------------------------------------------------------
_GO_ID_RE = re.compile(r"\[(GO:\d+)\]")


def _go_id(term):
    m = _GO_ID_RE.search(str(term))
    return m.group(1) if m else str(term)


def _membership(per_gene, cat_col, half_col="HalfLife"):
    """T1/2 > 0 の遺伝子について、(遺伝子名, T1/2, 用語集合) の DataFrame を返す。"""
    hl = pd.to_numeric(per_gene[half_col], errors="coerce")
    ok = hl > 0
    df = pd.DataFrame({
        "gene": per_gene.loc[ok, "Gene Name"].astype(str).values,
        "hl": hl[ok].astype(float).values,
        "terms": [set(_split_terms(c)) for c in per_gene.loc[ok, cat_col]],
    })
    return df


def mann_whitney(x, y):
    """両側 Mann–Whitney U 検定 (順位の正規近似、同順位補正、連続修正なし)。
    戻り値: (U_x, P, Cliff's δ = 2U/(n_x n_y) − 1)"""
    x = pd.Series(x, dtype=float); y = pd.Series(y, dtype=float)
    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        return float("nan"), float("nan"), float("nan")
    allv = pd.concat([x, y], ignore_index=True)
    ranks = allv.rank(method="average")
    u = float(ranks.iloc[:nx].sum() - nx * (nx + 1) / 2.0)
    n = nx + ny
    t = allv.value_counts().values
    tie = float(((t ** 3) - t).sum())
    var = nx * ny / 12.0 * ((n + 1) - tie / (n * (n - 1)))
    if var <= 0:
        return u, float("nan"), 2 * u / (nx * ny) - 1
    z = (u - nx * ny / 2.0) / math.sqrt(var)
    p = math.erfc(abs(z) / math.sqrt(2.0))
    return u, p, 2 * u / (nx * ny) - 1


def bh_qvalues(p):
    """Benjamini–Hochberg の q 値 (欠損は欠損のまま)。"""
    p = pd.Series(p, dtype=float)
    ok = p.notna()
    q = pd.Series(float("nan"), index=p.index)
    pv = p[ok].sort_values()
    m = len(pv)
    if m == 0:
        return q
    adj = (pv * m / pd.Series(range(1, m + 1), index=pv.index)).iloc[::-1].cummin().iloc[::-1].clip(upper=1.0)
    q[adj.index] = adj
    return q


def term_tests(per_gene, cat_col, items):
    """フィルタ通過済みの用語 (items) について、注釈あり vs なしを検定した表を返す。"""
    mem = _membership(per_gene, cat_col)
    rows = []
    for term, _ in items:
        m = mem["terms"].map(lambda s, t=term: t in s)
        a, b = mem.loc[m, "hl"], mem.loc[~m, "hl"]
        _, p, d = mann_whitney(a, b)
        rows.append({"GO term": term, "GO ID": _go_id(term), "n (annotated)": int(m.sum()),
                     "median T1/2 annotated (days)": float(a.median()) if len(a) else float("nan"),
                     "median T1/2 not annotated (days)": float(b.median()) if len(b) else float("nan"),
                     "Cliff's delta": d, "P (Mann-Whitney, two-sided)": p})
    res = pd.DataFrame(rows)
    if len(res):
        res["q (BH)"] = bh_qvalues(res["P (Mann-Whitney, two-sided)"])
        res = res.sort_values("Cliff's delta", ascending=False).reset_index(drop=True)
    return res


def jaccard_matrix(per_gene, cat_col, terms):
    """指定した用語どうしの Jaccard 係数 (共通タンパク質数 / 和集合) の正方行列。"""
    mem = _membership(per_gene, cat_col)
    sets = {t: set(mem.loc[mem["terms"].map(lambda s, t=t: t in s), "gene"]) for t in terms}
    mat = pd.DataFrame(index=list(terms), columns=list(terms), dtype=float)
    for a in terms:
        for b in terms:
            u = sets[a] | sets[b]
            mat.loc[a, b] = len(sets[a] & sets[b]) / len(u) if u else float("nan")
    return mat


def jaccard_pairs(per_gene, cat_col, terms):
    """用語の組ごとの重なり: 各用語のタンパク質数、共通数、Jaccard 係数、overlap 係数
    (共通数 / 小さい方の用語のタンパク質数; 入れ子の関係では 1 に近づく)。Jaccard 降順。"""
    mem = _membership(per_gene, cat_col)
    sets = {t: set(mem.loc[mem["terms"].map(lambda s, t=t: t in s), "gene"]) for t in terms}
    rows = []
    terms = list(terms)
    for i, a in enumerate(terms):
        for b in terms[i + 1:]:
            inter = len(sets[a] & sets[b]); union = len(sets[a] | sets[b])
            small = min(len(sets[a]), len(sets[b]))
            rows.append({"GO term A": a, "GO term B": b, "n A": len(sets[a]), "n B": len(sets[b]),
                         "n both": inter, "Jaccard": inter / union if union else float("nan"),
                         "Overlap coefficient": inter / small if small else float("nan")})
    return pd.DataFrame(rows).sort_values("Jaccard", ascending=False).reset_index(drop=True)


def conditional_test(per_gene, cat_col, term_id, exclude_id):
    """term_id を持つタンパク質と持たないタンパク質を、exclude_id を持つタンパク質を除いて比較。"""
    mem = _membership(per_gene, cat_col)
    has = lambda gid: mem["terms"].map(lambda s: any(_go_id(t) == gid for t in s))  # noqa: E731
    keep = ~has(exclude_id)
    m = has(term_id)
    a, b = mem.loc[m & keep, "hl"], mem.loc[~m & keep, "hl"]
    n_all = int(m.sum()); n_overlap = int((m & ~keep).sum())
    _, p, d = mann_whitney(a, b)
    return {"term": term_id, "excluded": exclude_id, "n annotated (all)": n_all,
            "n also annotated with excluded term": n_overlap,
            "fraction overlapping": n_overlap / n_all if n_all else float("nan"),
            "n annotated after exclusion": int(len(a)),
            "median T1/2 annotated (days)": float(a.median()) if len(a) else float("nan"),
            "median T1/2 other (days)": float(b.median()) if len(b) else float("nan"),
            "Cliff's delta": d, "P (Mann-Whitney, two-sided)": p}


def stats_excel_bytes(tests, jac, pairs, cond, settings):
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as w:
        tests.to_excel(w, sheet_name="Term tests", index=False)
        jac.to_excel(w, sheet_name="Jaccard matrix")
        pairs.to_excel(w, sheet_name="Jaccard pairs", index=False)
        pd.DataFrame([cond]).to_excel(w, sheet_name="Conditional test", index=False)
        settings.to_excel(w, sheet_name="Settings", index=False)
    return buf.getvalue()


def _high_low_indices(data, highlight_n):
    """各箱の median を基に、median が高い上位 N 個と低い下位 N 個の
    インデックス集合を返す。箱数が少なく重複する場合は上位を優先。"""
    n = len(data)
    k = max(0, min(int(highlight_n), n))
    if k == 0:
        return set(), set(), 0
    medians = [statistics.median(v) for v in data]
    order = sorted(range(n), key=lambda i: medians[i])  # median 昇順
    low_idx = set(order[:k])
    high_idx = set(order[-k:])
    low_idx -= high_idx                                  # 重複時は上位(赤)を優先
    return high_idx, low_idx, k


def select_high_low_subset(items, n):
    """median 上位 n 個 + 下位 n 個だけを抜き出し、median 降順に並べて返す。
    戻り値: (部分 items, 上位(赤) の個数)。"""
    high_idx, low_idx, _ = _high_low_indices([v for _, v in items], n)
    med = lambda i: statistics.median(items[i][1])  # noqa: E731
    high_sorted = sorted(high_idx, key=med, reverse=True)
    low_sorted = sorted(low_idx, key=med, reverse=True)
    subset = [items[i] for i in high_sorted + low_sorted]
    return subset, len(high_sorted)


def _wrap_label(text, max_len=20):
    """長い軸ラベルを 2 行に折り返す (単語境界で分割)。"""
    if len(text) <= max_len:
        return text
    width = max(max_len, (len(text) + 1) // 2)
    lines = textwrap.wrap(text, width=width, break_long_words=False,
                          break_on_hyphens=False)
    if len(lines) > 2:                       # 3 行以上は 2 行にまとめる
        lines = [lines[0], " ".join(lines[1:])]
    return "\n".join(lines)


def _format_label(term, hide_go_id):
    """表示用ラベルを作成。hide_go_id=True なら ' [GO:xxxx]' を除去。"""
    if hide_go_id:
        term = re.sub(r"\s*\[GO:\d+\]", "", term).strip()
    return _wrap_label(term)


def build_boxplot(items, category_label, highlight_n=10, only_high_low=False,
                  hide_go_id=False, show_all_points=False):
    """フィルタ済みの (用語, 値リスト) から箱ひげ図 Figure を作成して返す。
    median が高い上位 highlight_n 個を赤、低い下位 highlight_n 個を青で強調。
    only_high_low=True の場合は、その上位/下位のみ (最大 2×highlight_n 個) を描画。
    hide_go_id=True で軸ラベルの [GO:xxxx] を非表示、
    show_all_points=True で全データ点をプロット。"""
    if only_high_low and highlight_n:
        items, n_high = select_high_low_subset(items, highlight_n)
        high_idx = set(range(n_high))
        low_idx = set(range(n_high, len(items)))
    else:
        high_idx, low_idx, _ = _high_low_indices([v for _, v in items], highlight_n)
    labels = [_format_label(t, hide_go_id) for t, _ in items]
    data = [v for _, v in items]

    width = max(12.0, len(labels) * 0.42)
    fig, ax = plt.subplots(figsize=(width, 6.0))

    bp = ax.boxplot(
        data, patch_artist=True, widths=0.6,
        showfliers=not show_all_points,      # 全点表示時は外れ値マーカーを重複させない
        flierprops=dict(marker="o", markersize=2.5, markerfacecolor="#555555",
                        markeredgecolor="none", alpha=0.5),
    )
    # 箱ごとに塗り分け (上位=赤 / 下位=青 / その他=薄紫)。フチ・線は黒。
    for i, box in enumerate(bp["boxes"]):
        if i in high_idx:
            fill = HIGH_FILL
        elif i in low_idx:
            fill = LOW_FILL
        else:
            fill = BOX_FILL
        box.set(facecolor=fill, alpha=0.8, edgecolor="black", linewidth=1)
    # ひげ・キャップ・中央線は黒
    for whisker in bp["whiskers"]:
        whisker.set(color="black", linewidth=1)
    for cap in bp["caps"]:
        cap.set(color="black", linewidth=1)
    for median in bp["medians"]:
        median.set(color="black", linewidth=1.4)

    # 全データを点で表示 (横方向にジッターを付与、乱数は固定シードで再現可能)
    if show_all_points:
        rng = random.Random(0)
        xs, ys = [], []
        for i, vals in enumerate(data):
            for v in vals:
                xs.append(i + 1 + rng.uniform(-0.18, 0.18))
                ys.append(v)
        ax.scatter(xs, ys, s=5, color="#333333", alpha=0.3, linewidths=0,
                   zorder=3)

    ax.set_yscale("log")
    ax.set_ylabel("T1/2 (days)")
    jp = _jp_fontprops()
    if jp is not None:                       # 日本語フォントがある環境のみラベル表示
        ax.set_xlabel("項目", fontproperties=jp)
    ax.set_title(f"Gene Ontology ({category_label})")
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    # 強調した箱の項目名も同色に
    for i, lbl in enumerate(ax.get_xticklabels()):
        if i in high_idx:
            lbl.set_color(HIGH_EDGE)
        elif i in low_idx:
            lbl.set_color(LOW_EDGE)
    ax.grid(axis="y", which="both", color="#eeeeee")
    ax.set_axisbelow(True)

    handles = []                              # 凡例 (英語表記でフォント差を回避)
    if high_idx:
        handles.append(Patch(facecolor=HIGH_FILL, edgecolor="black",
                             label=f"Top {len(high_idx)} highest median"))
    if low_idx:
        handles.append(Patch(facecolor=LOW_FILL, edgecolor="black",
                             label=f"Top {len(low_idx)} lowest median"))
    if len(data) > len(high_idx) + len(low_idx):
        handles.append(Patch(facecolor=BOX_FILL, edgecolor="black",
                             label="Others"))
    if handles:
        ax.legend(handles=handles, loc="upper right", fontsize=8,
                  framealpha=0.9)

    fig.tight_layout()
    return fig


def figure_to_png(fig, dpi=150):
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    buffer.seek(0)
    return buffer.getvalue()


# ---------------------------------------------------------------------------
# Streamlit UI
# ---------------------------------------------------------------------------
def main():
    st.set_page_config(page_title="Protein List (GO解析用) 作成", layout="wide")
    st.title("🧬 Protein List (GO 解析用) 作成ツール")
    st.caption(
        "T-half シートの Gene Name と UniProtGO シートの From を照合し、"
        "Gene Ontology (cellular component / biological process / molecular "
        "function) と Entry を転記します。"
    )

    with st.expander("処理ルールの詳細", expanded=False):
        st.markdown(
            """
- **照合**: `T-half` の Gene Name (A列) ↔ `UniProtGO` の From (A列)
- **採用行の選択** (同名が複数ある場合):
  1. `Reviewed` (C列) が **reviewed** の行を最優先
  2. その中で **Length** (H列) が最大の行を採用
- **転記**:
  - C列 ← UniProtGO J列 (cellular component)
  - D列 ← UniProtGO I列 (biological process)
  - E列 ← UniProtGO L列 (molecular function)
  - F / G / H列 ← 採用行の Entry (B列) — それぞれ C / D / E に対応
- **重複なし一覧** (セミコロン区切りを分解・重複除去):
  - I列 ← C列の一覧 / J列 ← D列の一覧 / K列 ← E列の一覧
- **出現回数**: L列 ← I、M列 ← J、N列 ← K
            """
        )

    st.sidebar.header("設定")
    th_name = st.sidebar.text_input("T-half シート名", value="T-half")
    go_name = st.sidebar.text_input("UniProtGO シート名", value="UniProtGO")

    st.sidebar.markdown("---")
    st.sidebar.subheader("箱ひげ図フィルタ")
    min_count = st.sidebar.number_input(
        "最小タンパク質数 (この数以上)", min_value=1, value=20, step=1,
    )
    max_pct = st.sidebar.slider(
        "最大割合 (全タンパク質に対する %、この割合以下)",
        min_value=10, max_value=100, value=90, step=1,
    )
    top_n = st.sidebar.number_input(
        "表示上限 (頻度の高い上位 N 件、0 = すべて)", min_value=0, value=0, step=10,
    )
    highlight_n = st.sidebar.number_input(
        "Median 上位/下位を強調する件数 (赤=上位, 青=下位)",
        min_value=0, value=10, step=1,
    )
    only_high_low = st.sidebar.checkbox(
        "Median 上位/下位のみをグラフ化 (2×上記件数のみ表示)", value=False,
    )
    hide_go_id = st.sidebar.checkbox(
        "GO ID を非表示 (用語名のみ表示)", value=False,
    )
    show_all_points = st.sidebar.checkbox(
        "全データを点で表示 (外れ値以外も表示)", value=False,
    )
    max_frac = max_pct / 100.0

    uploaded = st.file_uploader(
        "Excel ファイルをアップロード (.xlsx)", type=["xlsx", "xlsm"]
    )

    if uploaded is None:
        st.info("上のボックスから ProteinList の Excel をアップロードしてください。")
        return

    # 新しいファイルがアップロードされたら前回の結果を破棄
    file_sig = (uploaded.name, getattr(uploaded, "size", None))
    if st.session_state.get("file_sig") != file_sig:
        st.session_state.pop("result", None)
        st.session_state.pop("png_cache", None)
        st.session_state["file_sig"] = file_sig

    # シート一覧を確認
    try:
        xl = pd.ExcelFile(uploaded)
    except Exception as exc:  # noqa: BLE001
        st.error(f"Excel を読み込めませんでした: {exc}")
        return

    if th_name not in xl.sheet_names or go_name not in xl.sheet_names:
        st.error(
            f"シートが見つかりません。ファイル内のシート: {xl.sheet_names}。"
            "サイドバーでシート名を修正してください。"
        )
        return

    th_df = xl.parse(th_name)
    go_df = xl.parse(go_name)
    st.success(
        f"読み込み完了 — T-half: {len(th_df):,} 行 / UniProtGO: {len(go_df):,} 行"
    )

    if st.button("▶ 実行 (Protein list を作成)", type="primary"):
        with st.spinner("照合・転記を実行中..."):
            try:
                out_df, stats = process(th_df, go_df)
            except Exception as exc:  # noqa: BLE001
                st.error(f"処理中にエラーが発生しました: {exc}")
                return
        # 箱ひげ図用に「用語 -> T1/2 値」を一度だけ計算して保持
        term_vals = {
            col: collect_term_halflives(out_df, col)
            for _, col in BOX_CATEGORIES
        }
        total_proteins = int(
            (pd.to_numeric(out_df["HalfLife"], errors="coerce") > 0).sum()
        )
        st.session_state["result"] = {
            "out_df": out_df, "go_df": go_df, "stats": stats,
            "term_vals": term_vals, "total": total_proteins,
        }
        st.session_state.pop("png_cache", None)

    res = st.session_state.get("result")
    if not res:
        return

    out_df = res["out_df"]
    stats = res["stats"]

    # サマリー
    st.subheader("サマリー")
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("対象遺伝子数", f"{stats['total_genes']:,}")
    c2.metric(
        "照合できた遺伝子", f"{stats['matched_genes']:,}",
        delta=f"未照合 {stats['total_genes'] - stats['matched_genes']:,}",
    )
    c3.metric("CC ユニーク語数", f"{stats['unique_cc']:,}")
    c4.metric(
        "BP / MF ユニーク語数",
        f"{stats['unique_bp']:,} / {stats['unique_mf']:,}",
    )

    # プレビュー
    st.subheader("プレビュー — T-half (A〜H: 遺伝子単位)")
    st.dataframe(out_df.iloc[:, :8].head(50), use_container_width=True)

    st.subheader("プレビュー — 重複なし一覧と出現回数 (I〜N)")
    preview_cols = OUT_HEADERS[8:]
    st.dataframe(
        out_df[preview_cols][out_df["Unique GO cellular component"] != ""].head(50),
        use_container_width=True,
    )

    # Excel ダウンロード
    with st.spinner("ダウンロード用 Excel を作成中..."):
        excel_bytes = to_excel_bytes(out_df, res["go_df"])
    st.download_button(
        label="⬇ 結果 Excel をダウンロード (ProteinList_GO.xlsx)",
        data=excel_bytes,
        file_name="ProteinList_GO.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        type="primary",
    )
    st.caption(
        "ダウンロードファイルには、完成した T-half シートと元の UniProtGO "
        "シートが含まれます。"
    )

    # 箱ひげ図
    st.markdown("---")
    st.header("📦 箱ひげ図 — 各 GO 用語の T1/2 分布")
    filter_note = (
        f"フィルタ: 出現頻度 **{min_count} 以上** かつ "
        f"**全タンパク質 ({res['total']:,}) の {max_pct}% 以下** "
        f"(= {max_frac * res['total']:.0f} 以下)、頻度の高い順"
    )
    if top_n:
        filter_note += f" / 上位 **{top_n}** 件のみ"
    st.caption(filter_note + "。縦軸 T1/2 は対数目盛です。")
    if highlight_n:
        if only_high_low:
            st.caption(
                f"🔴 median が高い上位 **{highlight_n}** 個 と "
                f"🔵 median が低い下位 **{highlight_n}** 個 **のみ** を "
                "median 降順で表示しています。"
            )
        else:
            st.caption(
                f"🔴 median が高い上位 **{highlight_n}** 個 / "
                f"🔵 median が低い下位 **{highlight_n}** 個 を強調表示しています。"
            )

    render_boxplots(res, min_count, max_frac, top_n, highlight_n, only_high_low,
                    hide_go_id, show_all_points)

    render_term_statistics(res, min_count, max_frac, top_n, highlight_n)


def render_term_statistics(res, min_count, max_frac, top_n, highlight_n):
    """補足解析: 用語ごとの検定 (BH-FDR)、重なり (Jaccard)、条件付き比較。"""
    st.markdown("---")
    st.header("📊 GO 用語の検定と重なりの確認（補足）")
    st.caption(
        "箱ひげ図と同じフィルタを通過した用語について、注釈ありと注釈なしのタンパク質の T1/2 を "
        "両側 Mann–Whitney U 検定（正規近似・同順位補正）で比較し、Cliff's δ（正 = 長寿命）と "
        "Benjamini–Hochberg の q 値を求めます。"
    )
    out_df = res["out_df"]
    labels = {lab: col for lab, col in BOX_CATEGORIES}
    cat_label = st.selectbox("カテゴリ", list(labels), index=0, key="stat_cat")
    cat_col = labels[cat_label]
    items = filter_terms(res["term_vals"][cat_col], res["total"], min_count, max_frac, top_n)
    if not items:
        st.warning("この条件に合う GO 用語がありません。")
        return
    tests = term_tests(out_df, cat_col, items)
    q = tests["q (BH)"]
    c1, c2, c3 = st.columns(3)
    c1.metric("検定した用語", f"{len(tests):,}")
    d = tests["Cliff's delta"]
    c2.metric("q < 0.05（長寿命）", f"{int(((q < 0.05) & (d > 0)).sum()):,}")
    c3.metric("q < 0.05（短寿命）", f"{int(((q < 0.05) & (d < 0)).sum()):,}")
    st.dataframe(tests, use_container_width=True, height=360)

    st.subheader("重なり（Jaccard 係数）— Median 上位/下位の用語")
    k = max(1, int(highlight_n) if highlight_n else 10)
    sub, _ = select_high_low_subset(items, k)
    terms = [t for t, _ in sub]
    jac = jaccard_matrix(out_df, cat_col, terms)
    pairs = jaccard_pairs(out_df, cat_col, terms)
    st.caption("Jaccard = 両方の用語に注釈されたタンパク質数 / どちらかに注釈されたタンパク質数。"
               "Overlap 係数 = 共通数 / 小さい方の用語のタンパク質数（入れ子の関係では 1 に近づく）。")
    st.dataframe(pairs.head(30), use_container_width=True)

    st.subheader("条件付き比較（ある用語の注釈を持つタンパク質を除外）")
    cc1, cc2 = st.columns(2)
    term_id = cc1.text_input("検定する用語の GO ID", value="GO:0043209", help="既定: myelin sheath")
    excl_id = cc2.text_input("除外する用語の GO ID", value="GO:0005739", help="既定: mitochondrion")
    cond = conditional_test(out_df, cat_col, term_id.strip(), excl_id.strip())
    st.dataframe(pd.DataFrame([cond]).T.rename(columns={0: "Value"}).astype(str), use_container_width=True)

    settings = pd.DataFrame({
        "Setting": ["Category", "Min proteins", "Max fraction", "Top N", "Highlight N (Jaccard terms)",
                    "Test", "Multiple testing", "Total proteins (T1/2 > 0)"],
        "Value": [cat_label, min_count, max_frac, top_n, k,
                  "two-sided Mann–Whitney U (normal approximation, tie-corrected); Cliff's delta",
                  "Benjamini–Hochberg across tested terms", res["total"]],
    })
    st.download_button(
        "⬇ 検定・重なりの結果をダウンロード (GO_term_statistics.xlsx)",
        data=stats_excel_bytes(tests, jac, pairs, cond, settings),
        file_name=f"GO_term_statistics_{cat_label.replace(' ', '_')}.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        key="dl_term_stats",
    )


def render_boxplots(res, min_count, max_frac, top_n, highlight_n,
                    only_high_low=False, hide_go_id=False,
                    show_all_points=False):
    """3 カテゴリの箱ひげ図を描画し、PNG ダウンロードボタンを表示する。"""
    png_cache = st.session_state.setdefault("png_cache", {})
    subset_mode = bool(only_high_low and highlight_n)
    for cat_label, cat_col in BOX_CATEGORIES:
        items = filter_terms(
            res["term_vals"][cat_col], res["total"], min_count, max_frac, top_n
        )
        if subset_mode:
            n_plotted = min(2 * int(highlight_n), len(items))
            st.subheader(f"{cat_label} — median 上位/下位 {n_plotted:,} 用語")
        else:
            st.subheader(f"{cat_label} — {len(items):,} 用語")
        if not items:
            st.warning(
                "この条件に合う GO 用語がありません。フィルタを調整してください。"
            )
            continue

        cache_key = (cat_col, int(min_count), round(max_frac, 4),
                     int(top_n), int(highlight_n), subset_mode,
                     bool(hide_go_id), bool(show_all_points))
        cached = png_cache.get(cat_col)
        if cached and cached[0] == cache_key:
            png = cached[1]
        else:
            n_boxes = min(2 * int(highlight_n), len(items)) if subset_mode \
                else len(items)
            with st.spinner(f"{cat_label} の箱ひげ図を作成中... ({n_boxes:,} 箱)"):
                png = figure_to_png(
                    build_boxplot(items, cat_label, highlight_n, subset_mode,
                                  hide_go_id, show_all_points)
                )
            png_cache[cat_col] = (cache_key, png)

        _show_image(png)
        st.download_button(
            label=f"⬇ {cat_label} の箱ひげ図をダウンロード (PNG)",
            data=png,
            file_name=f"boxplot_{cat_label.replace(' ', '_')}.png",
            mime="image/png",
            key=f"dl_box_{cat_col}",
        )


def _show_image(png_bytes):
    """Streamlit のバージョン差を吸収して画像を横幅いっぱいに表示。"""
    try:
        st.image(png_bytes, use_container_width=True)
    except TypeError:
        st.image(png_bytes, use_column_width=True)


if __name__ == "__main__":
    main()
