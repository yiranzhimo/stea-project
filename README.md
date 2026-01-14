# STEA

STEA is a deconvolution tool for spatial transcriptomics datasets.

## Requirements
- Python 3.9+
- [Scanpy](https://scanpy.readthedocs.io) ≥ 1.9
- pandas, numpy, seaborn, matplotlib
- scipy, scikit-learn, statsmodels

Install the dependencies with pip (feel free to pin versions to your stack):
```bash
pip install scanpy pandas numpy seaborn matplotlib scipy scikit-learn statsmodels
```

### Conda Environment Setup
Prefer Conda? Create an isolated env, install core scientific stacks via `conda-forge`, then fall back to pip for anything missing:
```bash
# 1) Create and activate the environment
conda create -n stea python=3.10 -c conda-forge
conda activate stea

# 2) Install common deps with conda
conda install -c conda-forge scanpy seaborn matplotlib scipy scikit-learn statsmodels

# 3) Fill the gaps (if any) with pip inside the same env
pip install pandas numpy
```
If you already have a Scanpy-ready env, simply `conda activate <env>` and skip to the installation step below.

## Installation

### Install from GitHub
```bash
# Install directly from GitHub
pip install git+https://github.com/yiranzhimo/stea-project.git

# Or clone and install in development mode
git clone https://github.com/yiranzhimo/stea-project.git
cd stea-project
pip install -e .
```



## Tips & Troubleshooting
- Ensure `adata.X` still contains raw counts before calling `qcplot`; it will log-transform data if necessary and store a copy in `adata.layers['counts']` when absent.
- For KL-based functions, large `num_permutations` values improve stability but increase runtime; start with 500–1000 and scale up as needed.
- The cosine similarity helper materializes dense matrices; for very large datasets consider subsetting genes or converting to dense only for the genes of interest.
