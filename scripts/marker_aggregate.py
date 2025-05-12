"""
Aggregate metatranscriptome CPM counts by functional marker genes.

Required input files
====================
1. **Per‑marker HMMER domtblout** – one file per entry in `MARKER_DOMTBL`, produced by
   `hmmsearch --domtblout <file> <profile.hmm> <proteins.faa>`.
2. **Gene‑level CPM matrix** – tab‑separated table where the first column is
   `OMRGC_ID` (or any unique gene identifier) and the remaining columns are
   samples.  Path configured in `COUNTS_TSV` below.

Output files generated
======================
* **all_sum_CPM.csv** – tidy summary matrix (rows = samples, columns =
  marker genes) containing summed CPM counts.  Overwritten if it already exists.

--------------------------------------------------------------------------
This script does **one thing**:
    • Reads each domtblout, extracts every gene (ORF / OM‑RGC ID) that matches a
      given marker profile.
    • Pulls CPM values for those genes from the huge counts table.
    • Sums CPMs per sample for every marker.

Usage
-----
Adjust the dictionaries right below (``MARKER_DOMTBL``, ``SYNONYM_MAP``) and the
``COUNTS_TSV`` path, then run:

    python marker_aggregate.py

Dependencies: pandas ≥ 1.0
"""

from __future__ import annotations
import re
import pandas as pd
from pathlib import Path
from typing import Dict, List

###############################################################################
# 1. User-editable configuration
###############################################################################

# Mapping   marker_name  →  domtblout file produced by `hmmsearch --domtblout`
# Add or remove entries as needed.
MARKER_DOMTBL: Dict[str, str] = {
    "nod":   "nod.tbl",
    # --- aerobic markers ------------------------------------------------------
    "cox1":  "cox1.domtbl",
    "cox3":  "cox3.domtbl",
    "cytbd": "cytbd.domtbl",
    "sod":   "sod.domtbl",
    "pmoA":  "pmoA.domtbl",
    # --- anaerobic markers ----------------------------------------------------
    "nap":   "nap.domtbl",
    "narG":  "narG.domtbl",
    "nir":   "nir.domtbl",
    "nor":   "nor.domtbl",
    "nos":   "nos.domtbl",
}

# HMM names that should be collapsed onto a single “marker” label
# (use *lower-case* keys for case-insensitive matching).
SYNONYM_MAP: Dict[str, str] = {
    "ccoN".lower(): "cox1",
    "ctaD".lower(): "cox1",
    "cyoB".lower(): "cox1",
    "cydA".lower(): "cytbd",
    "cox1".lower(): "cox1",
    "cox3".lower(): "cox3",
    "Sod_Fe_N".lower(): "sod",
    "sodA".lower(): "sod",
    "dsrab": "dsr",   # example
    "napA": "nap",
    "narG": "narG",
    "nirS": "nir",
    "nirK": "nir",
    "norB": "nor",
    "nosZ": "nos",
    # … extend as needed …
}

# Path to the huge CPM matrix (tab-separated, one row per OM-RGC gene, one column per sample)
COUNTS_TSV = "counts_cpm.tsv"  # <-- change me if needed

###############################################################################
# 2. Helper functions
###############################################################################

def parse_domtbl(path: str | Path, marker: str) -> pd.DataFrame:
    """Return a 1-column DataFrame of OM-RGC IDs that hit *marker*.

    Assumes the **target/sequence** ID is in the **first** whitespace-separated
    column of each non-comment line (standard `--domtblout` layout).
    """
    ids: List[str] = []
    with open(path, "r") as handle:
        for ln in handle:
            if ln.startswith("#"):  # skip comments
                continue
            seq_id = ln.split(maxsplit=1)[0]
            # “clean” the ID if it contains trailing domain coords (e.g. “…/1-250”)
            seq_id = re.sub(r"/\d+.*", "", seq_id)
            ids.append(seq_id)

    df = pd.DataFrame({"OMRGC_ID": ids})
    df["marker"] = marker  # annotate
    return df.drop_duplicates()

###############################################################################
# 3. Build the gene-to-marker lookup table
###############################################################################

all_hits: List[pd.DataFrame] = []
for marker, file in MARKER_DOMTBL.items():
    fp = Path(file)
    if not fp.exists():
        raise FileNotFoundError(f"Expected domtblout for marker '{marker}' at {file}")
    df_marker = parse_domtbl(fp, marker)
    all_hits.append(df_marker)

dedup = pd.concat(all_hits, ignore_index=True).drop_duplicates("OMRGC_ID")

# Apply synonym mapping (case-insensitive)
idx = dedup["marker"].str.lower().map(lambda x: SYNONYM_MAP.get(x, x))
dedup["marker"] = idx

###############################################################################
# 4. Stream-read the massive CPM matrix, filtering to relevant IDs
###############################################################################

chunks: List[pd.DataFrame] = []
CHUNKSIZE = 200_000  # adjust based on RAM
cols_to_keep: List[str] | None = None  # will discover header lazily

for chunk in pd.read_csv(COUNTS_TSV, sep="\t", chunksize=CHUNKSIZE):
    # Lazily determine which columns we actually want (all samples + OMRGC_ID)
    if cols_to_keep is None:
        samples = [c for c in chunk.columns if c != "OMRGC_ID"]
        cols_to_keep = ["OMRGC_ID"] + samples
    chunk = chunk[cols_to_keep]
    chunk = chunk[chunk["OMRGC_ID"].isin(dedup["OMRGC_ID"])].copy()
    if not chunk.empty:
        chunks.append(chunk)

if not chunks:
    raise RuntimeError("No overlapping OMRGC_IDs between CPM table and HMMER hits!")

counts_df = pd.concat(chunks, ignore_index=True)

###############################################################################
# 5. Aggregate CPMs per sample × marker
###############################################################################

# Map each gene to its marker label
counts_df["marker"] = counts_df["OMRGC_ID"].map(dedup.set_index("OMRGC_ID")["marker"])

# Melt → long, aggregate, pivot → wide
long_df = counts_df.melt(id_vars=["OMRGC_ID", "marker"], var_name="sample", value_name="CPM")
long_df = long_df[long_df["CPM"] > 0]  # ignore zero rows to speed up groupby

final = (
    long_df
    .groupby(["sample", "marker"], sort=False)["CPM"].sum()
    .unstack(fill_value=0)
    .sort_index(axis=1)
)

###############################################################################
# 6. Write result
###############################################################################

final.to_csv("all_sum_CPM.csv")
print("✔ Aggregated CPM matrix written to all_sum_CPM.csv")
