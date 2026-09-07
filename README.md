# Protein List Builder for Gene Ontology (GO) Analysis

A single-file [Streamlit](https://streamlit.io/) application that builds a
Gene Ontology (GO) annotated protein list from a protein half-life dataset and a
UniProt GO export, and visualises the protein half-life (T1/2) distribution of
each GO term as a box plot.

The application was developed to support the analysis reported in the associated
publication (see [Citation](#citation)).

---

## Features

- **Annotation transfer** — matches each gene in the half-life table against a
  UniProt GO export and copies the cellular component, biological process and
  molecular function annotations, together with the UniProt accession (`Entry`).
- **Deterministic row selection** — when several UniProt records map to the same
  gene name, reviewed (Swiss-Prot) records are preferred and, among those, the
  record with the greatest sequence length is used, so results are reproducible.
- **Term inventory** — splits the semicolon-separated GO strings, deduplicates
  them and reports how many proteins carry each term.
- **Box plots** — for each GO namespace, plots the T1/2 distribution of the
  proteins carrying each term on a logarithmic axis, with the terms of highest
  and lowest median highlighted in red and blue.
- **Exports** — the annotated workbook downloads as `ProteinList_GO.xlsx` and
  each plot downloads as a PNG image.

---

## Installation

Python 3.9 or newer is required.

```bash
git clone <repository-url>
cd <repository-directory>

python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate

pip install -r requirements.txt
```

## Usage

```bash
streamlit run protein_list_go_app.py
```

Streamlit opens a browser tab (by default <http://localhost:8501>). Upload your
Excel workbook, adjust the filters in the sidebar if required, and press
**▶ Run (build protein list)**.

An example workbook is provided at
[`data/ProteinList_monkey_common.xlsx`](data/ProteinList_monkey_common.xlsx).

---

## Input format

The uploaded `.xlsx` (or `.xlsm`) workbook must contain two sheets. The sheet
names are configurable in the sidebar; the defaults are shown below.

### Sheet `T-half`

| Column | Header      | Description                                        |
| ------ | ----------- | -------------------------------------------------- |
| A      | `Gene Name` | Gene symbol used as the matching key                |
| B      | `HalfLife`  | Protein half-life T1/2 (days); non-positive and missing values are excluded from the plots |

### Sheet `UniProtGO`

A UniProt ID-mapping / retrieval export. The following columns are required
(additional columns are carried through untouched):

| Header                                | Description                                    |
| ------------------------------------- | ---------------------------------------------- |
| `From`                                | Gene symbol submitted to UniProt (matching key) |
| `Entry`                               | UniProt accession                              |
| `Reviewed`                            | `reviewed` / `unreviewed`                      |
| `Length`                              | Sequence length (amino acids)                  |
| `Gene Ontology (biological process)`  | Semicolon-separated GO terms                   |
| `Gene Ontology (cellular component)`  | Semicolon-separated GO terms                   |
| `Gene Ontology (molecular function)`  | Semicolon-separated GO terms                   |

---

## Output format

The generated `T-half` sheet contains columns A–N:

| Column | Header                         | Content                                             |
| ------ | ------------------------------ | --------------------------------------------------- |
| A      | `Gene Name`                    | copied from the input                               |
| B      | `HalfLife`                     | copied from the input                               |
| C      | `GO cellular component`        | from the selected UniProt record                    |
| D      | `GO biological process`        | from the selected UniProt record                    |
| E      | `GO molecular function`        | from the selected UniProt record                    |
| F–H    | `Entry (…)`                    | UniProt accession of the selected record for C/D/E  |
| I–K    | `Unique GO …`                  | deduplicated term list of C/D/E, in order of first appearance |
| L–N    | `Count (…)`                    | number of proteins carrying the term in I/J/K       |

The downloaded workbook also retains the original `UniProtGO` sheet.

---

## Method summary

1. **Matching** — `T-half!Gene Name` is matched against `UniProtGO!From` as text.
2. **Record selection** — for genes with multiple UniProt records:
   1. records with `Reviewed == "reviewed"` are retained if any exist;
   2. among the retained records, the one with the largest `Length` is selected.
3. **Term statistics** — GO strings are split on `;`, whitespace-trimmed and
   counted across all genes.
4. **Box plots** — proteins with a missing or non-positive T1/2 are excluded
   (the y-axis is logarithmic). Terms are retained when they are carried by at
   least `min_count` proteins (default 20) and by no more than `max_frac` of all
   proteins (default 90%), which removes both rare and near-universal terms.
   The terms with the highest and lowest medians are highlighted. Point jitter
   uses a fixed random seed, so figures are reproducible.

### Sidebar parameters

| Parameter | Default | Description |
| --------- | ------- | ----------- |
| T-half / UniProtGO sheet name | `T-half` / `UniProtGO` | names of the input sheets |
| Minimum number of proteins | 20 | lower frequency cut-off for a GO term |
| Maximum share | 90% | upper frequency cut-off, as a share of all proteins |
| Display limit | 0 (all) | keep only the N most frequent terms |
| Highlight count | 10 | number of highest/lowest medians to highlight |
| Plot only highest/lowest medians | off | restrict the plot to the highlighted terms |
| Hide GO IDs | off | strip `[GO:xxxxxxx]` from the axis labels |
| Show all data points | off | overlay every individual protein as a jittered point |

---

## Requirements

```
streamlit
pandas
openpyxl
matplotlib
```

Pinned lower bounds are listed in [`requirements.txt`](requirements.txt).

---

## Repository layout

```
.
├── protein_list_go_app.py                 # the Streamlit application
├── data/
│   └── ProteinList_monkey_common.xlsx     # example input workbook
├── tests/
│   └── test_protein_list_go_app.py        # unit tests for the processing logic
├── requirements.txt
├── CITATION.cff
├── CONTRIBUTING.md
├── LICENSE
└── README.md
```

---

## Tests

```bash
pip install pytest
pytest
```

The tests exercise the annotation-transfer and term-counting logic on small
synthetic inputs and do not require Streamlit to be running.

---

## Citation

If you use this software, please cite the associated publication and this
repository. Citation metadata is provided in
[`CITATION.cff`](CITATION.cff); GitHub renders it as a
*"Cite this repository"* button.

> Publication details are to be completed once the article is published; please
> update the `preferred-citation` block in `CITATION.cff` accordingly.

---

## License

Distributed under the BSD 3-Clause License. See [`LICENSE`](LICENSE) for the
full text.

Copyright (c) 2026, Takeda Pharmaceutical Company Limited. All rights reserved.

Gene Ontology annotations are obtained from [UniProt](https://www.uniprot.org/)
and the [Gene Ontology Consortium](https://geneontology.org/); please observe
their respective terms of use when redistributing derived data.
