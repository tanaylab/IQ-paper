# IceQream (IQ) paper companion repository

Welcome to the companion repository for the **“IceQream: Quantitative chromosome accessibility analysis using physical TF models”** manuscript. 

This repo contains everything required to:
1. Reproduce **figures** in the paper (R / Python notebooks).
2. Download the raw data and set-up genome databases.
3. Train and evaluate the deep-learning **benchmark models** compared in Figure 5.

---

## Directory layout

| Path | Contents |
|------|----------|
| `analysis/` | Jupyter notebooks that generate Figures 1–5. Each notebook is named after the figure it produces (e.g. `Figure3.ipynb`). |
| `benchmarks/` | Training & evaluation code for the deep-learning benchmarks (Figure 5). See the dedicated `benchmarks/README.md` for a quick-start guide. |
| `code/` | Lightweight helper scripts (mostly R) for data download & preprocessing. In particular, `code/download_data.R` fetches all public tracks used throughout the paper. |

---

## Getting started

```bash
# clone the repo 
git clone https://github.com/tanaylab/IQ-paper.git

cd IQ-paper

# 1. Create a fresh conda / venv environment (optional but recommended)

# 2. Install Python & R dependencies
pip install -r benchmarks/requirements.txt  # deep-learning deps
# R packages are listed at the top of each notebook / script

# 3. Download data & build genome DBs (creates ./data/)
Rscript code/download_data.R
```

The download script downloads approximately **15 GB** of data required for all analyses.

---

## Reproducing the figures

All primary figure notebooks live in `analysis/`. Launch them via JupyterLab or VS Code:

```bash
jupyter lab  # then open analysis/Figure1.ipynb …
```

Each notebook is fully self-contained and will write any intermediate results to `./output/`.

---

## Running the benchmarks (Figure 5)

Deep-learning baselines and the IQ ensemble are defined in `benchmarks/`.

```bash
# Example: train every model & collect R² scores
cd benchmarks
bash commands.sh                  # runs >10 jobs sequentially
python scripts/collect_r2.py --runs output/* > all_models_r2.csv
```

Hardware tips, expected runtimes and per-model instructions are documented in `benchmarks/README.md`.

---

## Citing IceQream

If you use IceQream or the resources in this repository, please cite:

```
Bercovich A, Lifshitz A *et al.*
"IceQream: Quantitative chromosome accessibility analysis using physical TF models" (2025)
```

---

For questions, bug reports or feature requests **please open a GitHub issue**.
