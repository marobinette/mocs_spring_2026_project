# Parameter Sweep Documentation

## Overview

Three distinct parameter sweeps exist across this codebase. They form a coherent stack: the analytical sweep (`find_tricritical_points_joint`) maps the phase boundary for the constant-ω model; the fixed-point tracer (`collect_fixed_points_by_lam`) produces full branch portraits at a given ν; and the trajectory sweep (`run_sweep`) is the brute-force approach that works even where analytics fail — which is exactly where the kernel model lives.

The analytical tricritical condition (Eq. 16 of the paper) does not transfer to the diversity-tension kernel because ω is now `w(n, i)` — a function of group state — rather than a scalar that factors out of the derivative calculations. The trajectory sweep is therefore the primary investigation tool, not a fallback.

---

## Sweep 1 — `run_sweep` (`kerenel.ipynb`, cells 6–9)

**Scope:** 2D grid over (λ, ν). The main sweep for this project.

**Locations:**
| Role | File | Cell / Line |
|---|---|---|
| Grid and IC setup | `kerenel.ipynb` | Cell 6 (id `99416687`) |
| `run_sweep` definition | `kerenel.ipynb` | Cell 7 (id `f8a69430`) |
| Run A — kernel execution | `kerenel.ipynb` | Cell 8 (id `1f89f0a4`) |
| Run B — baseline execution | `kerenel.ipynb` | Cell 9 (id `197c2234`) |

### Grid

| Parameter | Range | Spacing | Local N | Notes |
|---|---|---|---|---|
| λ (transmission rate) | 1e-4 to 1e-1 | log | 15 | Expand downward if transition sits below grid |
| ν (synergy exponent) | 1.0 to 12.0 | linear | 15 | 1 = linear contagion, 12 = strong synergy |

### Initial Conditions

Two initial conditions per cell (Laurent's two-IC strategy):

- `I0_low = 1e-3` — invasion attempt; system climbs to the first stable branch from below
- `I0_high = 0.99` — persistence attempt; system relaxes to the first stable branch from above

### Mechanism

For each (λ, ν) cell, integrates `integrate_I_traj_kernel` twice via LSODA to stationarity (`T_MAX=20`, `TRAJ_POINTS=500`). Records `I_traj[-1]` (the final infected fraction) from each run. Returns three `(N_LAM, N_NU)` arrays:

- `I_low` — stationary I\* reached from I0_low
- `I_high` — stationary I\* reached from I0_high
- `delta = I_high - I_low` — **bistability indicator**: large where two stable branches coexist, collapses to zero where only one stable state exists

### Two Runs

**Run A — Diversity-tension kernel:**
```python
w_func = w_diversity_tension  # w(n, i, alpha) = alpha * (i/n) * (1 - i/n)
w_args = (alpha,)             # alpha = 30.0
```
Switching rate peaks at φ = i/n = 0.5 (maximally mixed groups) and is zero at homogeneous endpoints (φ = 0 or 1). Saved to `Files/sweep_kernel_only.npz`.

**Run B — Constant-ω baseline:**
```python
w_func = w_constant           # w(n, i, omega) = omega
w_args = (omega_scalar,)      # omega_scalar = 5.0
```
Reproduces the paper's scalar-ω model on the same network. Saved alongside Run A in `Files/sweep_local.npz`.

### Output

Three-panel heatmap (`figures/sweep_bistability_comparison.png`):
1. Kernel δI
2. Baseline δI
3. Difference (kernel − baseline)

The heatmap of δI is the computational substitute for the analytical tricritical line — it maps where the interesting regimes live without requiring the derivative calculation.

### Integration Settings

```python
MU          = 1.0
TRAJ_POINTS = 500
T_MAX       = 200.0   # increased from 20.0 — see note below
```

These are lightweight local settings. The VACC sweep should use larger grids, denser I0 sampling in cells that show bistability signal, and multiple α values.

**Note on T_MAX:** With `mu=1` the natural relaxation timescale is 1 time unit. T=20 is only 20 characteristic times, which is too short for cells near the transition threshold where critical slowing down makes convergence arbitrarily slow. An unconverged trajectory from `I0_high` records a spuriously high `I*(high)`, producing a false δI signal. T=200 is safe for all cells in this grid. `TRAJ_POINTS` controls output resolution only (LSODA takes adaptive internal steps), so increasing T_MAX does not affect memory — the 500 output snapshots are simply spread over a longer window.

---

## Sweep 2 — `collect_fixed_points_by_lam` (`fixed_points.py`, line 99)

**Scope:** λ at fixed ν — a single 1D trace, not a 2D grid. Produces Figure 3-style branch portraits.

**Locations:**
| Role | File | Line |
|---|---|---|
| `collect_fixed_points_by_lam` | `wAMEs-main/src/wAMEs/fixed_points.py` | 99 |
| `find_fixed_points_for_lambda` (called per λ) | `wAMEs-main/src/wAMEs/fixed_points.py` | 3 |
| `pack_per_lam_pts` (pack results for plotting) | `wAMEs-main/src/wAMEs/fixed_points.py` | 195 |
| `plot_rank_tracked_branches` (visualize branches) | `wAMEs-main/src/wAMEs/fixed_points.py` | 375 |

### Mechanism

For each λ in `lam_vals`, evaluates the self-consistent mean-field map M(r) on a log-grid of 4000 points (`r_min` to `r_max`), brackets sign changes of G(r) = M(r) − r, and refines each root with `brentq` (xtol=1e-12). Near-duplicate roots are removed. Each root is classified by slope M′(r\*): stable if slope < 1, unstable otherwise.

From each fixed point r\*, reconstructs the full stationary state and computes:
- `I*` — stationary infected fraction via `I_from_r`
- `P = 1 / Ỹ₄` — localization observable (inverse normalized fourth-order moment of group infection levels); large P means activity is spread broadly (delocalized), small P means concentrated in the largest groups (localized)

Fixed points are sorted by decreasing I\* and packed into dense arrays by `pack_per_lam_pts`. Branches are then tracked by rank index and plotted by `plot_rank_tracked_branches`, which draws stable segments as solid lines, unstable segments as dashed lines, and optionally marks the endpoints of unstable segments.

### Key Difference from `run_sweep`

This uses the **scalar-ω mean-field map** (`mf_map_w` in `core.py`), not the kernel vector field. It finds equilibria exactly rather than by simulating trajectories. It applies to the constant-ω wAMEs model as implemented in the library — not to the diversity-tension kernel.

The unstable branch (the dashed middle line in Figure 3) is directly observable here via root-finding. In `run_sweep`, the unstable branch is invisible by construction: `I*(high) ≠ I*(low)` at a given λ is the only signal that it exists.

---

## Sweep 3 — `find_tricritical_points_joint` (`thresholds.py`, line 362)

**Scope:** ν × w (synergy exponent × rewiring rate), searching for tricritical points in (ν, w, λ) space. Applies to the **constant-ω model only**.

**Locations:**
| Role | File | Line |
|---|---|---|
| `find_tricritical_points_joint` | `wAMEs-main/src/wAMEs/thresholds.py` | 362 |
| `tricritical_condition` (evaluates d²M/dr² at threshold) | `wAMEs-main/src/wAMEs/thresholds.py` | 232 |
| `invasion_threshold_w` (computes λc at each candidate point) | `wAMEs-main/src/wAMEs/thresholds.py` | 161 |

### Mechanism

Two complementary scans, combined and deduplicated:

**Scan A — Fix ν, scan w:**
- Outer grid: `n_nu_outer=41` points, linear from `nu_min` to `nu_max`
- Inner grid: `n_w_inner=400` points, log-spaced from `w_min` to `w_max`
- For each fixed ν, evaluates `tricritical_condition(ν, w)` = d²M/dr²|_{r=0, λ=λc} across the inner w grid, brackets sign changes, refines with `brentq`

**Scan B — Fix w, scan ν:**
- Outer grid: `n_w_outer=101` points, log-spaced
- Inner grid: `n_nu_inner=400` points, linear from `nu_inner_min` to `nu_inner_max`
- Symmetric to Scan A with ν and w swapped

At each candidate tricritical point (ν\*, w\*), also computes λc via `invasion_threshold_w`. Results from both scans are merged, deduplicated at `rtol_dup=1e-4`, and sorted by ν to yield arrays `(nu_vals, w_vals, lambda_vals)`.

### Why This Does Not Extend to the Kernel

The tricritical condition is derived by setting the second derivative of M to zero. In the constant-ω model, ω is a scalar and can be pulled outside the sums in the derivative calculation. In the kernel model, ω is `w(n, i)` — it has a different value for every (n, i) pair and cannot be factored out. Every sum that previously simplified now carries an extra (n, i) dependency. Extending this computation to the kernel is a significant physics project (see `param_sweep_guidance.md`, §4).

---

## Alignment Summary

| Guidance Recommendation | Code Status |
|---|---|
| Two ICs per cell (I0_low, I0_high) | ✅ `run_sweep` Cells 6–7 |
| λ log-spaced 1e-4 to 1e-1 | ✅ `np.logspace(-4, -1, 15)` |
| ν linear 1 to 12 | ✅ `np.linspace(1.0, 12.0, 15)` |
| δI as bistability proxy for tricritical line | ✅ `delta = I_high - I_low` |
| Kernel vs constant-ω baseline comparison | ✅ Run A / Run B |
| Tricritical line not computed for kernel | ✅ `find_tricritical_points_joint` has no kernel version |
| Extend λ downward if transition sits below grid | Pending — depends on sweep results |
| Larger VACC sweep, denser I0 near bistability | Pending — VACC grid not yet designed |

---

## Next Steps (from `param_sweep_guidance.md`)

1. Read the local heatmap — where is δI nonzero in the kernel panel vs baseline?
2. Validate the extremes — does ν ≈ 1 give a continuous transition? Does ν ≈ 10 show bistability?
3. Calibrate the λ range — extend downward if the transition boundary sits below 1e-4
4. Design the VACC sweep — larger grid, multiple α values, denser I0 sampling in bistable cells
5. Draft the abstract once at least one concrete finding about bistability is in hand
