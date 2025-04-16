# 440finalproject

# Overview

This repository provides a complete pipeline for the multi‑omic analysis of dark oxygen production in anoxic waters via nitrate/oxygen dismutation (NOD). It includes:

- **Exploratory notebooks**: Jupyter notebooks demonstrating data QC, visualization, and initial hypothesis testing.
- **Data processing scripts**: Python modules for parsing metagenomic and metatranscriptomic read tables, presence/absence calls, and normalization.
- **HMMER pipeline**: Shell and Python wrappers to build and apply a consensus HMM profile for nod genes using HMMER3 (Eddy, 2017). 
- **Phylogenetic analysis files**: MEGA6 input files and exported Newick trees for nod gene family reconstruction (Tamura et al., 2013).
- **Plotting scripts**: Standalone `.py` files (e.g., `plot_nod_vs_oxygen.py`) that generate figures used in the manuscript.

> **Citations:**
>
> - Eddy, S. R. (2017). HMMER: biosequence analysis using profile hidden Markov models. [*http://hmmer.org/*](http://hmmer.org/).
> - Tamura, K., Stecher, G., Peterson, D., Filipski, A., & Kumar, S. (2013). MEGA6: Molecular Evolutionary Genetics Analysis version 6.0. *Molecular Biology and Evolution*.

---

# Data

All analyses use three large environmental metagenomic/metatranscriptomic datasets:

1. **TARA Oceans** (Salazar et al., 2019): global ocean transect, samples filtered for O₂ < 20 μM.
2. **Cariaco Basin** (Geller‑McGrath et al., 2023): anoxic particle‐associated and free‐living fractions.
3. **North American Rivers** (Borton et al., 2024): fertilizer‐runoff freshwater samples.

Raw read count tables and accompanying metadata are too large to host here. You can download them from:

- TARA Oceans project portal: [*https://doi.org/10.1016/j.cell.2019.10.014*](https://doi.org/10.1016/j.cell.2019.10.014)
- Cariaco Basin data repository: [*https://doi.org/10.1038/s41467-023-36026-w*](https://doi.org/10.1038/s41467-023-36026-w)
- GROWdb North American Rivers: [*https://doi.org/10.1038/s41586-024-08240-z*](https://doi.org/10.1038/s41586-024-08240-z)

Processed presence/absence matrices and filtered datasets are included under `data/processed/`.

---

# Folder structure

```plaintext
├── data/ # Filtered & formatted input tables (CSV)│
│               
├── notebooks/          # Jupyter notebooks for exploration & QC
│
├── src/                # Core Python modules
│ 
├── scripts/            # Generated figure outputs (PNG, PDF)
│   └── plot_nod_vs_oxygen.py
├── figures/            
│
├── results/            # Analysis results & tables used in manuscript
│
├── docs/               # Supplementary methods & protocol descriptions
│
└── README.md           # This file
```

Each subdirectory contains a brief `README.md` explaining its contents.

---
# Installation:

To run the code make sure that the .csv file is in the same folder as the .py file
This code was developed using python 3.13.3.

   Dependencies required:

   - Python ≥ 3.8
   - numpy, pandas, matplotlib, scikit-learn, scipy, biopython

 **How to Install HMMER3**

   ```bash
   # On macOS with Homebrew:
   brew install hmmer

   # Or from source:
   wget http://eddylab.org/software/hmmer/hmmer-3.3.2.tar.gz
   tar xzf hmmer-3.3.2.tar.gz && cd hmmer-3.3.2
   ./configure && make && make install
   ```

5. **How to Install MEGA6**

   Download the MEGA6 desktop application from [*https://www.megasoftware.net/*](https://www.megasoftware.net/) to perform phylogenetic reconstructions.

## Usage

- Place raw data files in `data/raw/` (follow naming conventions in `data/README.md`).
- Run the HMMER pipeline:
  ```bash
  python src/hmm_pipeline.py --input data/raw/ --output data/processed/
  ```
- Generate figures:
  ```bash
  python scripts/plot_nod_vs_oxygen.py
  ```

For detailed instructions, see the README files in each subfolder.

