# Posture Toward Diversity Determines Bistability in Higher-Order Contagion

Spring 2026 MOCS project. Extends the constant-$\omega$ wAME framework of Lamata-Otín et al. (2026) with three state-dependent group composition turnover kernels (diversity-tension, tension-shifted, diversity-preferred) and characterizes the bistability response across the $(\lambda, \nu)$ plane on three networks.

## Repository structure

```
.
├── main.tex                    Paper source (LaTeX)
├── higher_order_contagion.bib  Bibliography
├── vacc_sweep.py               Cluster sweep script (one α, full (λ,ν) grid)
├── submit_sweep.sh             SLURM submission for the VACC cluster
├── analysis.ipynb              Post-processing: heatmaps, ridgelines, summary table
├── replicate_sweep.ipynb       Reproduce a single-α heatmap locally
├── network_data_eda.ipynb      Network statistics + paper Table 1 source
├── discussion.md               Discussion outline
├── Data/
│   └── group_statistics.txt    p(n) and g(k) for every network
├── Files/
│   └── vacc/                   .npz outputs from the sweep (one file per α)
└── figures/                    Generated PNGs and CSVs (paper-ready)
```

## Reproducing results

### A single heatmap (local, ~2 min)

`replicate_sweep.ipynb` runs one $(\lambda, \nu)$ sweep at a chosen `(NETWORK, KERNEL, ALPHA)` and plots the resulting bistability heatmap. Edit the constants at the top:

```python
NETWORK = "Thiers13"             # Thiers13, Synthetic_poisson_k2, Synthetic_poisson_k3, ...
KERNEL  = "diversity_tension"    # diversity_tension, tension_shifted, mixture
ALPHA   = 0.1
```

Defaults are 20×15 grid resolution with `T_MAX=300` to stay fast. Match the cluster runs with 50×40 and `T_MAX=1000` for full resolution.

### Full sweeps (cluster)

`vacc_sweep.py` runs one α across the full $(\lambda, \nu)$ grid. `submit_sweep.sh` parallelizes across α values on UVM's VACC. Outputs land in `Files/vacc/` as `.npz`.

### Post-processing

- `analysis.ipynb` reads `Files/vacc/`, produces the figures in `figures/`, and writes `peak_bistability_summary.csv`.
- `network_data_eda.ipynb` produces `network_statistics.csv` (the structural summary in the Methods).

## Networks studied

| Network | Type | $\langle k\rangle$ | Notes |
|---|---|---|---|
| Synthetic Poisson, $\langle k\rangle = 2$ | Synthetic | 2.00 | Shares Thiers13's $p(n)$, Poisson $g(k)$ |
| Synthetic Poisson, $\langle k\rangle = 3$ | Synthetic | 3.00 | Shares Thiers13's $p(n)$, Poisson $g(k)$ |
| Thiers13 | Empirical | 1.04 | French high-school contact network |

The synthetics share Thiers13's group-size distribution by design and differ only in $g(k)$, isolating the effect of $\langle k\rangle$ on bistability.

## Requirements

Python 3 with `numpy`, `scipy`, `matplotlib`, `pandas`. Notebooks expect to be run from the project root.
