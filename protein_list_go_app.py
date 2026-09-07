"""
Protein List Builder for Gene Ontology (GO) Analysis
====================================================

A Streamlit application that reads a ProteinList Excel workbook and transfers
Gene Ontology annotations from the "UniProtGO" sheet onto the "T-half" sheet,
producing a protein list ready for GO analysis.

Workflow
--------
1. Match Gene Name (column A) of "T-half" against From (column A) of
   "UniProtGO". When several rows share the same Gene Name, the representative
   row is selected as follows:
     (1) rows whose Reviewed (column C) equals "reviewed" take priority
     (2) among those, the row with the largest Length (column H) is used
2. Transfer the following fields from the selected row:
     column C <- UniProtGO column J  Gene Ontology (cellular component)
     column D <- UniProtGO column I  Gene Ontology (biological process)
     column E <- UniProtGO column L  Gene Ontology (molecular function)
     columns F/G/H <- Entry (column B) of the selected row
                      (F for C, G for D, H for E)
3. Split the semicolon-separated strings in columns C/D/E, remove duplicates,
   and write the resulting term lists to columns I/J/K (semicolons dropped).
4. Write the number of occurrences of each term in the corresponding column to
   columns L/M/N (I->L, J->M, K->N).
5. For each GO namespace (cellular component / biological process /
   molecular function), draw a box plot of the T1/2 (HalfLife) values of the
   proteins carrying each GO term, on a logarithmic axis. Terms are filtered by
   frequency (>= 20 proteins and <= 90% of all proteins). The plots can be
   downloaded as PNG files.

Requirements:
    streamlit, pandas, openpyxl, matplotlib

Usage:
    streamlit run protein_list_go_app.py

SPDX-License-Identifier: BSD-3-Clause
"""

import io
import random
import re
import statistics
import textwrap
from collections import Counter

import matplotlib
matplotlib.use("Agg")  # non-interactive backend for headless/server environments
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
import pandas as pd  # noqa: E402
import streamlit as st  # noqa: E402


# Default box plot colours (light purple)
BOX_FILL = "#9a9ae6"
BOX_EDGE = "#5b5bd6"
BOX_MEDIAN = "#3b3ba8"

# Highlight colours for the highest/lowest medians (high = red, low = blue)
HIGH_FILL = "#e06666"   # red
HIGH_EDGE = "#cc0000"
LOW_FILL = "#5b9bd5"    # blue
LOW_EDGE = "#1155cc"

# The three GO namespaces plotted as box plots (label, column name in out_df)
BOX_CATEGORIES = [
    ("cellular component", "GO cellular component"),
    ("biological process", "GO biological process"),
    ("molecular function", "GO molecular function"),
]

# ---------------------------------------------------------------------------
# Source column names of the UniProtGO sheet
# ---------------------------------------------------------------------------
COL_FROM = "From"
COL_ENTRY = "Entry"
COL_REVIEWED = "Reviewed"
COL_LENGTH = "Length"
COL_BP = "Gene Ontology (biological process)"   # column I
COL_CC = "Gene Ontology (cellular component)"    # column J
COL_MF = "Gene Ontology (molecular function)"    # column L

# Headers of the generated T-half sheet (A..N)
OUT_HEADERS = [
    "Gene Name",                       # A
    "HalfLife",                        # B
    "GO cellular component",           # C  <- UniProtGO J
    "GO biological process",           # D  <- UniProtGO I
    "GO molecular function",           # E  <- UniProtGO L
    "Entry (cellular component)",      # F  <- Entry of selected row (for C)
    "Entry (biological process)",      # G  <- Entry of selected row (for D)
    "Entry (molecular function)",      # H  <- Entry of selected row (for E)
    "Unique GO cellular component",    # I  <- deduplicated terms of C
    "Unique GO biological process",    # J  <- deduplicated terms of D
    "Unique GO molecular function",    # K  <- deduplicated terms of E
    "Count (cellular component)",      # L  <- occurrence count of I
    "Count (biological process)",      # M  <- occurrence count of J
    "Count (molecular function)",      # N  <- occurrence count of K
]


# ---------------------------------------------------------------------------
# Processing logic
# ---------------------------------------------------------------------------
def _split_terms(cell):
    """Split a GO cell on ';' and return the non-empty, stripped terms."""
    if cell is None or (isinstance(cell, float) and pd.isna(cell)):
        return []
    return [t.strip() for t in str(cell).split(";") if t.strip()]


def select_best_row(candidates):
    """Pick a single representative UniProtGO row for one gene.

    (1) prefer Reviewed == 'reviewed'  (2) then take the largest Length.
    """
    is_reviewed = (
        candidates[COL_REVIEWED].astype(str).str.strip().str.lower() == "reviewed"
    )
    reviewed = candidates[is_reviewed]
    pool = reviewed if len(reviewed) else candidates
    lengths = pd.to_numeric(pool[COL_LENGTH], errors="coerce").fillna(-1)
    return pool.loc[lengths.idxmax()]


def _val(row, col):
    """Return the cell value, or an empty string when it is missing."""
    if col not in row or pd.isna(row[col]):
        return ""
    return row[col]


def process(th_df, go_df):
    """Return the completed T-half data frame (A..N) plus summary statistics."""
    # Verify that the required columns are present
    missing = [c for c in (COL_FROM, COL_ENTRY, COL_REVIEWED, COL_LENGTH,
                           COL_BP, COL_CC, COL_MF) if c not in go_df.columns]
    if missing:
        raise ValueError(
            "The following columns were not found in the UniProtGO sheet: "
            + ", ".join(missing)
        )
    if th_df.shape[1] < 2:
        raise ValueError(
            "The T-half sheet must contain two columns: Gene Name and HalfLife."
        )

    gene_col = th_df.columns[0]      # A: Gene Name
    half_col = th_df.columns[1]      # B: HalfLife

    # Group the candidate rows by their From value
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
        # A single row is selected, so its Entry is the source of C/D/E -> F=G=H
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

    # Deduplicated term list + occurrence counts (per column)
    def unique_counts(series):
        counter = Counter()
        for cell in series:
            for term in _split_terms(cell):
                counter[term] += 1
        ordered = list(counter.keys())          # keep first-appearance order
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
    """Return an Excel workbook (as bytes) with the T-half and UniProtGO sheets."""
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        out_df.to_excel(writer, sheet_name="T-half", index=False)
        go_df.to_excel(writer, sheet_name="UniProtGO", index=False)
    buffer.seek(0)
    return buffer.getvalue()


# ---------------------------------------------------------------------------
# Box plots
# ---------------------------------------------------------------------------
def collect_term_halflives(per_gene, cat_col, half_col="HalfLife"):
    """Map each GO term to the T1/2 values of the proteins carrying it.

    Rows whose T1/2 is missing or <= 0 are dropped (the axis is logarithmic).
    """
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
    """Filter terms by frequency and return (term, values) sorted by frequency.

    Condition: min_count <= (number of proteins with the term)
               <= max_frac * total number of proteins.
    """
    hi = max_frac * total_proteins
    items = [(t, v) for t, v in term_vals.items() if min_count <= len(v) <= hi]
    items.sort(key=lambda x: len(x[1]), reverse=True)
    if top_n and top_n > 0:
        items = items[:top_n]
    return items


def _high_low_indices(data, highlight_n):
    """Return the index sets of the N highest and N lowest medians.

    When there are too few boxes and the two sets overlap, the high set wins.
    """
    n = len(data)
    k = max(0, min(int(highlight_n), n))
    if k == 0:
        return set(), set(), 0
    medians = [statistics.median(v) for v in data]
    order = sorted(range(n), key=lambda i: medians[i])  # ascending median
    low_idx = set(order[:k])
    high_idx = set(order[-k:])
    low_idx -= high_idx                                  # on overlap, high (red) wins
    return high_idx, low_idx, k


def select_high_low_subset(items, n):
    """Keep only the n highest and n lowest medians, sorted by median (desc).

    Returns: (subset of items, number of high (red) entries).
    """
    high_idx, low_idx, _ = _high_low_indices([v for _, v in items], n)
    med = lambda i: statistics.median(items[i][1])  # noqa: E731
    high_sorted = sorted(high_idx, key=med, reverse=True)
    low_sorted = sorted(low_idx, key=med, reverse=True)
    subset = [items[i] for i in high_sorted + low_sorted]
    return subset, len(high_sorted)


def _wrap_label(text, max_len=20):
    """Wrap a long axis label onto two lines (breaking at word boundaries)."""
    if len(text) <= max_len:
        return text
    width = max(max_len, (len(text) + 1) // 2)
    lines = textwrap.wrap(text, width=width, break_long_words=False,
                          break_on_hyphens=False)
    if len(lines) > 2:                       # collapse three or more lines to two
        lines = [lines[0], " ".join(lines[1:])]
    return "\n".join(lines)


def _format_label(term, hide_go_id):
    """Build the display label; hide_go_id=True strips ' [GO:xxxx]'."""
    if hide_go_id:
        term = re.sub(r"\s*\[GO:\d+\]", "", term).strip()
    return _wrap_label(term)


def build_boxplot(items, category_label, highlight_n=10, only_high_low=False,
                  hide_go_id=False, show_all_points=False):
    """Build a box plot Figure from the filtered (term, values) pairs.

    The highlight_n terms with the highest medians are drawn in red and the
    highlight_n terms with the lowest medians in blue.
    only_high_low=True plots only those terms (at most 2 x highlight_n).
    hide_go_id=True removes [GO:xxxx] from the axis labels and
    show_all_points=True overlays every individual data point.
    """
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
        showfliers=not show_all_points,      # avoid drawing outliers twice
        flierprops=dict(marker="o", markersize=2.5, markerfacecolor="#555555",
                        markeredgecolor="none", alpha=0.5),
    )
    # Fill each box (high = red / low = blue / other = light purple); outlines black.
    for i, box in enumerate(bp["boxes"]):
        if i in high_idx:
            fill = HIGH_FILL
        elif i in low_idx:
            fill = LOW_FILL
        else:
            fill = BOX_FILL
        box.set(facecolor=fill, alpha=0.8, edgecolor="black", linewidth=1)
    # Whiskers, caps and medians in black
    for whisker in bp["whiskers"]:
        whisker.set(color="black", linewidth=1)
    for cap in bp["caps"]:
        cap.set(color="black", linewidth=1)
    for median in bp["medians"]:
        median.set(color="black", linewidth=1.4)

    # Overlay every data point with horizontal jitter (fixed seed = reproducible)
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
    ax.set_title(f"Gene Ontology ({category_label})")
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    # Colour the tick labels of the highlighted boxes to match
    for i, lbl in enumerate(ax.get_xticklabels()):
        if i in high_idx:
            lbl.set_color(HIGH_EDGE)
        elif i in low_idx:
            lbl.set_color(LOW_EDGE)
    ax.grid(axis="y", which="both", color="#eeeeee")
    ax.set_axisbelow(True)

    handles = []
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
    st.set_page_config(page_title="Protein List Builder for GO Analysis",
                       layout="wide")
    st.title("🧬 Protein List Builder for GO Analysis")
    st.caption(
        "Matches Gene Name in the T-half sheet against From in the UniProtGO "
        "sheet and transfers the Gene Ontology annotations (cellular component "
        "/ biological process / molecular function) together with the Entry "
        "accession."
    )

    with st.expander("Processing rules", expanded=False):
        st.markdown(
            """
- **Matching**: Gene Name (column A) of `T-half` ↔ From (column A) of `UniProtGO`
- **Row selection** (when several rows share a gene name):
  1. rows whose `Reviewed` (column C) is **reviewed** take priority
  2. among those, the row with the largest **Length** (column H) is used
- **Transfer**:
  - column C ← UniProtGO column J (cellular component)
  - column D ← UniProtGO column I (biological process)
  - column E ← UniProtGO column L (molecular function)
  - columns F / G / H ← Entry (column B) of the selected row — matching C / D / E
- **Deduplicated term lists** (semicolon-separated values split and deduplicated):
  - column I ← terms of column C / column J ← terms of column D / column K ← terms of column E
- **Occurrence counts**: column L ← I, column M ← J, column N ← K
            """
        )

    st.sidebar.header("Settings")
    th_name = st.sidebar.text_input("T-half sheet name", value="T-half")
    go_name = st.sidebar.text_input("UniProtGO sheet name", value="UniProtGO")

    st.sidebar.markdown("---")
    st.sidebar.subheader("Box plot filters")
    min_count = st.sidebar.number_input(
        "Minimum number of proteins (at least)", min_value=1, value=20, step=1,
    )
    max_pct = st.sidebar.slider(
        "Maximum share (% of all proteins, at most)",
        min_value=10, max_value=100, value=90, step=1,
    )
    top_n = st.sidebar.number_input(
        "Display limit (top N most frequent terms, 0 = all)",
        min_value=0, value=0, step=10,
    )
    highlight_n = st.sidebar.number_input(
        "Number of highest/lowest medians to highlight (red = high, blue = low)",
        min_value=0, value=10, step=1,
    )
    only_high_low = st.sidebar.checkbox(
        "Plot only the highest/lowest medians (2 x the number above)",
        value=False,
    )
    hide_go_id = st.sidebar.checkbox(
        "Hide GO IDs (show term names only)", value=False,
    )
    show_all_points = st.sidebar.checkbox(
        "Show all data points (not only outliers)", value=False,
    )
    max_frac = max_pct / 100.0

    uploaded = st.file_uploader(
        "Upload an Excel file (.xlsx)", type=["xlsx", "xlsm"]
    )

    if uploaded is None:
        st.info("Upload your ProteinList Excel workbook using the box above.")
        return

    # Discard previous results when a new file is uploaded
    file_sig = (uploaded.name, getattr(uploaded, "size", None))
    if st.session_state.get("file_sig") != file_sig:
        st.session_state.pop("result", None)
        st.session_state.pop("png_cache", None)
        st.session_state["file_sig"] = file_sig

    # Check the available sheets
    try:
        xl = pd.ExcelFile(uploaded)
    except Exception as exc:  # noqa: BLE001
        st.error(f"The Excel file could not be read: {exc}")
        return

    if th_name not in xl.sheet_names or go_name not in xl.sheet_names:
        st.error(
            f"Sheet not found. Sheets in this file: {xl.sheet_names}. "
            "Please correct the sheet names in the sidebar."
        )
        return

    th_df = xl.parse(th_name)
    go_df = xl.parse(go_name)
    st.success(
        f"Loaded — T-half: {len(th_df):,} rows / UniProtGO: {len(go_df):,} rows"
    )

    if st.button("▶ Run (build protein list)", type="primary"):
        with st.spinner("Matching and transferring annotations..."):
            try:
                out_df, stats = process(th_df, go_df)
            except Exception as exc:  # noqa: BLE001
                st.error(f"An error occurred during processing: {exc}")
                return
        # Compute "term -> T1/2 values" once and keep it for the box plots
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

    # Summary
    st.subheader("Summary")
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Genes processed", f"{stats['total_genes']:,}")
    c2.metric(
        "Genes matched", f"{stats['matched_genes']:,}",
        delta=f"unmatched {stats['total_genes'] - stats['matched_genes']:,}",
    )
    c3.metric("Unique CC terms", f"{stats['unique_cc']:,}")
    c4.metric(
        "Unique BP / MF terms",
        f"{stats['unique_bp']:,} / {stats['unique_mf']:,}",
    )

    # Preview
    st.subheader("Preview — T-half (A–H: per gene)")
    st.dataframe(out_df.iloc[:, :8].head(50), use_container_width=True)

    st.subheader("Preview — deduplicated terms and counts (I–N)")
    preview_cols = OUT_HEADERS[8:]
    st.dataframe(
        out_df[preview_cols][out_df["Unique GO cellular component"] != ""].head(50),
        use_container_width=True,
    )

    # Excel download
    with st.spinner("Preparing the Excel file for download..."):
        excel_bytes = to_excel_bytes(out_df, res["go_df"])
    st.download_button(
        label="⬇ Download result Excel (ProteinList_GO.xlsx)",
        data=excel_bytes,
        file_name="ProteinList_GO.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        type="primary",
    )
    st.caption(
        "The downloaded file contains the completed T-half sheet and the "
        "original UniProtGO sheet."
    )

    # Box plots
    st.markdown("---")
    st.header("📦 Box plots — T1/2 distribution per GO term")
    filter_note = (
        f"Filter: frequency **>= {min_count}** and "
        f"**<= {max_pct}% of all proteins ({res['total']:,})** "
        f"(= {max_frac * res['total']:.0f} or fewer), sorted by frequency"
    )
    if top_n:
        filter_note += f" / top **{top_n}** terms only"
    st.caption(filter_note + ". The T1/2 axis uses a logarithmic scale.")
    if highlight_n:
        if only_high_low:
            st.caption(
                f"Showing **only** the 🔴 **{highlight_n}** terms with the "
                f"highest medians and the 🔵 **{highlight_n}** terms with the "
                "lowest medians, in descending order of median."
            )
        else:
            st.caption(
                f"Highlighting the 🔴 **{highlight_n}** terms with the highest "
                f"medians and the 🔵 **{highlight_n}** terms with the lowest "
                "medians."
            )

    render_boxplots(res, min_count, max_frac, top_n, highlight_n, only_high_low,
                    hide_go_id, show_all_points)


def render_boxplots(res, min_count, max_frac, top_n, highlight_n,
                    only_high_low=False, hide_go_id=False,
                    show_all_points=False):
    """Draw the box plots for the three namespaces with PNG download buttons."""
    png_cache = st.session_state.setdefault("png_cache", {})
    subset_mode = bool(only_high_low and highlight_n)
    for cat_label, cat_col in BOX_CATEGORIES:
        items = filter_terms(
            res["term_vals"][cat_col], res["total"], min_count, max_frac, top_n
        )
        if subset_mode:
            n_plotted = min(2 * int(highlight_n), len(items))
            st.subheader(
                f"{cat_label} — {n_plotted:,} terms with highest/lowest medians"
            )
        else:
            st.subheader(f"{cat_label} — {len(items):,} terms")
        if not items:
            st.warning(
                "No GO term matches these criteria. Please adjust the filters."
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
            with st.spinner(
                f"Building the {cat_label} box plot... ({n_boxes:,} boxes)"
            ):
                png = figure_to_png(
                    build_boxplot(items, cat_label, highlight_n, subset_mode,
                                  hide_go_id, show_all_points)
                )
            png_cache[cat_col] = (cache_key, png)

        _show_image(png)
        st.download_button(
            label=f"⬇ Download the {cat_label} box plot (PNG)",
            data=png,
            file_name=f"boxplot_{cat_label.replace(' ', '_')}.png",
            mime="image/png",
            key=f"dl_box_{cat_col}",
        )


def _show_image(png_bytes):
    """Display an image full width, tolerating Streamlit API differences."""
    try:
        st.image(png_bytes, use_container_width=True)
    except TypeError:
        st.image(png_bytes, use_column_width=True)


if __name__ == "__main__":
    main()
