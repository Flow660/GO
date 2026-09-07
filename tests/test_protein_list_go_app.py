"""Unit tests for the processing logic of protein_list_go_app.

Streamlit is stubbed out so that the tests can run without a running server.

Run with:  pytest
"""

import importlib.util
import sys
import types
from pathlib import Path

import pandas as pd
import pytest

APP_PATH = Path(__file__).resolve().parents[1] / "protein_list_go_app.py"


@pytest.fixture(scope="module")
def app():
    """Import the application module with a stubbed streamlit dependency."""
    if "streamlit" not in sys.modules:
        stub = types.ModuleType("streamlit")
        stub.__getattr__ = lambda name: (lambda *a, **k: None)
        sys.modules["streamlit"] = stub
    spec = importlib.util.spec_from_file_location("protein_list_go_app", APP_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def th_df():
    return pd.DataFrame(
        {"Gene Name": ["AAA", "BBB", "CCC"], "HalfLife": [1.0, 10.0, 100.0]}
    )


@pytest.fixture
def go_df():
    return pd.DataFrame(
        {
            "From": ["AAA", "AAA", "AAA", "BBB"],
            "Entry": ["P00001", "P00002", "P00003", "Q00001"],
            "Reviewed": ["unreviewed", "reviewed", "reviewed", "reviewed"],
            "Entry Name": ["a1", "a2", "a3", "b1"],
            "Length": [999, 100, 300, 200],
            "Gene Ontology (biological process)": [
                "bp1 [GO:0000001]",
                "bp1 [GO:0000001]; bp2 [GO:0000002]",
                "bp2 [GO:0000002]",
                "bp1 [GO:0000001]",
            ],
            "Gene Ontology (cellular component)": [
                "cc1 [GO:0000011]",
                "cc1 [GO:0000011]",
                "cc1 [GO:0000011]; cc2 [GO:0000012]",
                "cc2 [GO:0000012]",
            ],
            "Gene Ontology (molecular function)": [
                "mf1 [GO:0000021]",
                "mf1 [GO:0000021]",
                "mf2 [GO:0000022]",
                "mf1 [GO:0000021]",
            ],
        }
    )


def test_split_terms(app):
    assert app._split_terms("a [GO:1]; b [GO:2] ;") == ["a [GO:1]", "b [GO:2]"]
    assert app._split_terms(None) == []
    assert app._split_terms(float("nan")) == []


def test_select_best_row_prefers_reviewed_then_longest(app, go_df):
    best = app.select_best_row(go_df[go_df["From"] == "AAA"])
    # P00001 is longer but unreviewed, so the longest reviewed record wins.
    assert best["Entry"] == "P00003"


def test_process_transfers_annotations(app, th_df, go_df):
    out, stats = app.process(th_df, go_df)

    assert list(out.columns) == app.OUT_HEADERS
    assert stats["total_genes"] == 3
    assert stats["matched_genes"] == 2          # CCC has no UniProt record

    aaa = out.iloc[0]
    assert aaa["Entry (cellular component)"] == "P00003"
    assert aaa["GO cellular component"] == "cc1 [GO:0000011]; cc2 [GO:0000012]"
    # F, G and H all carry the accession of the single selected record.
    assert (
        aaa["Entry (cellular component)"]
        == aaa["Entry (biological process)"]
        == aaa["Entry (molecular function)"]
    )

    # Unmatched genes keep their half-life but have empty annotations.
    ccc = out.iloc[2]
    assert ccc["Gene Name"] == "CCC"
    assert ccc["GO biological process"] == ""


def test_process_counts_unique_terms(app, th_df, go_df):
    out, stats = app.process(th_df, go_df)

    terms = list(out["Unique GO cellular component"])
    counts = list(out["Count (cellular component)"])
    mapping = {t: c for t, c in zip(terms, counts) if t != ""}

    assert mapping == {"cc1 [GO:0000011]": 1, "cc2 [GO:0000012]": 2}
    assert stats["unique_cc"] == 2


def test_process_requires_expected_columns(app, th_df, go_df):
    with pytest.raises(ValueError, match="UniProtGO"):
        app.process(th_df, go_df.drop(columns=["Reviewed"]))


def test_collect_term_halflives_excludes_non_positive(app):
    per_gene = pd.DataFrame(
        {
            "HalfLife": [2.0, 0.0, None, 8.0],
            "GO cellular component": ["cc1", "cc1", "cc1", "cc1; cc2"],
        }
    )
    term_vals = app.collect_term_halflives(per_gene, "GO cellular component")
    assert term_vals == {"cc1": [2.0, 8.0], "cc2": [8.0]}


def test_filter_terms_applies_both_bounds(app):
    term_vals = {
        "rare": [1.0] * 2,
        "good": [1.0] * 5,
        "ubiquitous": [1.0] * 10,
    }
    items = app.filter_terms(term_vals, total_proteins=10, min_count=3,
                             max_frac=0.9)
    assert [t for t, _ in items] == ["good"]


def test_filter_terms_top_n_keeps_most_frequent(app):
    term_vals = {"a": [1.0] * 3, "b": [1.0] * 8, "c": [1.0] * 5}
    items = app.filter_terms(term_vals, total_proteins=100, min_count=1,
                             max_frac=1.0, top_n=2)
    assert [t for t, _ in items] == ["b", "c"]


def test_high_low_indices(app):
    data = [[1.0], [5.0], [10.0], [50.0]]
    high, low, k = app._high_low_indices(data, 1)
    assert k == 1 and high == {3} and low == {0}


def test_select_high_low_subset_orders_by_median(app):
    items = [("a", [1.0]), ("b", [100.0]), ("c", [10.0]), ("d", [50.0])]
    subset, n_high = app.select_high_low_subset(items, 1)
    assert n_high == 1
    assert [t for t, _ in subset] == ["b", "a"]


def test_format_label_hides_go_id(app):
    assert app._format_label("nucleus [GO:0005634]", True) == "nucleus"
    # Labels of 20 characters or fewer are left on a single line.
    assert app._format_label("nucleus [GO:0005634]", False) == \
        "nucleus [GO:0005634]"


def test_wrap_label_splits_long_labels_over_two_lines(app):
    wrapped = app._wrap_label(
        "positive regulation of transcription by RNA polymerase II"
    )
    assert wrapped.count("\n") == 1
    assert wrapped.replace("\n", " ") == \
        "positive regulation of transcription by RNA polymerase II"


def test_to_excel_bytes_roundtrip(app, th_df, go_df):
    import io

    out, _ = app.process(th_df, go_df)
    data = app.to_excel_bytes(out, go_df)
    sheets = pd.ExcelFile(io.BytesIO(data)).sheet_names
    assert sheets == ["T-half", "UniProtGO"]


def test_build_boxplot_returns_figure(app):
    items = [("cc1", [1.0, 2.0, 3.0]), ("cc2", [10.0, 20.0, 30.0])]
    fig = app.build_boxplot(items, "cellular component", highlight_n=1)
    png = app.figure_to_png(fig)
    assert png.startswith(b"\x89PNG")
