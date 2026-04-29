# Localization Analysis: Allegiance Kernel vs Scalar Baseline
**Date:** 2026-04-21  
**Network:** CNS (Corporate Network Structure)  
**Notebook:** `kernel.ipynb`, cell `586e43f9`

---

## 1. Model Overview

The wAME (weighted Approximate Master Equations) framework tracks the joint distribution of group states `f_{n,i}` — the fraction of groups of size `n` containing `i` infected members — alongside the per-membership susceptibility `s_m`. Infection spreads within groups; nodes can also switch groups at a rate controlled by `w`.

**Core state variables:**
- `f_{n,i}`: fraction of groups of size `n` with `i` infected members, shape `(nmax+1, nmax+1)`
- `s_m`: probability that a node with `m` group memberships is susceptible, shape `(mmax+1,)`
- `I = Σ_m (1 - s_m) g_m`: global infected fraction (scalar)

---

## 2. Network Parameters (CNS)

| Parameter | Symbol | Value | Description |
|---|---|---|---|
| Max group size | `nmax` | 20 | Largest group size in CNS network |
| Max membership | `mmax` | 233 | Maximum number of groups a node belongs to |
| Mean membership | `<m>` | ~2.33 | Average groups per node (multi-group overlap) |
| Group-size dist. | `p_n` | empirical | Probability a group has size `n` |
| Membership dist. | `g_m` | empirical | Probability a node has `m` memberships |

The CNS network was chosen because it has significant multi-group membership (`<m> > 1`), which is required for the wAME group-switching dynamics to be non-trivial.

---

## 3. Dynamical Parameters

| Parameter | Symbol | Value | Description |
|---|---|---|---|
| Recovery rate | `mu` (μ) | 1.0 | Rate at which infected nodes recover |
| Synergy exponent | `nu` (ν) | 9.5 | Controls nonlinearity of infection rate |
| Allegiance scale | `scale` | 5.0 | Amplitude of the allegiance switching kernel |
| Scalar rewiring rate | `w` | 5.0 | Baseline group-switching rate (non-kernel) |
| Lambda range | `lam_vals` | logspace(-6, -3, 50) | Transmission prefactor sweep |

### 3.1 Infection rate

```
β(n, i) = λ · i^ν
```

With `ν = 9.5`, the rate is strongly superlinear in the number of infected group members — a synergistic contagion model. Infection is negligible when few members are infected and explosive when the group is mostly infected.

### 3.2 Mean-field infection pressure

```
r = [Σ_{n,i} β(n,i) · (n-i) · f_{n,i} · p_n] / [Σ_{n,i} (n-i) · f_{n,i} · p_n]
```

`r` is the per-contact rate at which a susceptible node acquires infection from a random group partner.

---

## 4. Switching Rate Functions

### 4.1 Scalar baseline

```python
w(n, i) = w   # constant for all (n, i)
```

The switching rate does not depend on group composition. Used as the analytical reference case.

### 4.2 Allegiance kernel (`w_allegiance`)

```python
def w_allegiance(n, i, scale):
    phi = i / n
    return scale * 4.0 * phi * (1.0 - phi)
```

- Peaks at `phi = 0.5` (maximally mixed groups): nodes flee mixed groups
- Zero at `phi = 0` (all susceptible) and `phi = 1` (all infected): nodes stay in homogeneous groups
- `scale` controls the amplitude relative to `w = scale = 5.0`

This models **group allegiance**: nodes are loyal to groups that share their health state and leave groups where they are in the minority.

---

## 5. Core Functions

### `core.py`

| Function | Signature | Description |
|---|---|---|
| `build_inf_mat` | `(lam, nu, nmax)` | Builds infection matrix `β(n,i) = λ·i^ν` for all valid `(n,i)` pairs |
| `switching_matrix` | `(w_func, nmax, args)` | Builds switching matrix by evaluating `w_func(n,i,*args)` for all `(n,i)` |
| `compute_In_from_fni` | `(fni, nmax)` | Computes `I_n = Σ_i (i/n) · f_{n,i}` — mean infected fraction within size-`n` groups |
| `y4_tilde` | `(In, pn_use)` | Returns `Ỹ₄ = Σ_n p_n I_n⁴ / (Σ_n p_n I_n²)²` |
| `normalize_group_distribution` | `(pn, nmin=2)` | Zeros out `p_n` for `n < nmin` and renormalizes |
| `state_from_mf_w` | `(r, mu, w, inf_mat, state_meta)` | Analytically reconstructs stationary `(s_m, f_{n,i})` from scalar `r`, `w` |
| `mf_map_w` | `(r, mu, w, inf_mat, state_meta)` | Evaluates the self-consistent map `M(r)` — used for scalar fixed-point search |
| `infected_fraction` | `(sm, gm)` | Returns `I = Σ_m (1-s_m) g_m` |

### `fixed_points.py`

| Function | Signature | Description |
|---|---|---|
| `find_fixed_points_for_lambda` | `(inf_mat, mu, w, state_meta, r_min, r_max, n_grid)` | Scans log grid for sign changes of `M(r)-r`, refines roots with Brent's method, classifies stability from `M'(r*)` |
| `collect_fixed_points_by_lam` | `(lam_vals, nu, nmax, mu, w, state_meta, m_arr, gm, pn_filtered, r_min, r_max, n_grid_root)` | Loops over `lam_vals`, calls fixed-point finder, reconstructs state, computes `I*` and `P*` at each root |
| `pack_per_lam_pts` | `(per_lam_pts, lam_vals)` | Converts list-of-dicts to dense arrays of shape `(Kmax, N)` for plotting |

### `temporal_dynamics.py`

| Function | Signature | Description |
|---|---|---|
| `integrate_I_traj_kernel` | `(lam, w_func, state_meta, nmax, mmax, gm, mu, nu, w_args, I0, traj_points, t_max)` | Integrates the ODE forward in time with kernel `w_func`; returns `(t, I_traj, fni_traj)` |
| `initialize` | `(state_meta, initial_density)` | Sets `s_m = 1 - I0` uniformly and draws `f_{n,i}` from Binomial(n, I0) |
| `vector_field_w_kernel` | `(v, t, inf_mat, w_mat, state_meta, mu)` | JIT-compiled RHS of the ODE system using a precomputed switching matrix |

---

## 6. Localization Measure

### 6.1 Effective participation ratio P

```
P = (Σ_n p_n I_n²)² / (Σ_n p_n I_n⁴)   ∈  [p_{n_max}, 1]
```

Implemented as `P = 1.0 / y4_tilde(In, pn_filtered)`.

- **P = 1**: infection is uniformly distributed across all group sizes (delocalized)
- **P = p_{n_max}**: infection is concentrated in groups of the largest size only (maximally localized)

### 6.2 Epidemic threshold λ_c

The minimum `λ` at which infection can invade from near-zero. Detected as the first `λ` where `I_ker[I0=1e-15]` exceeds `I_cut = 1e-6` (kernel) or where the scalar low-IC stable branch is non-trivial (scalar).

### 6.3 Localization transition λ*

The value of `λ` at which `dP / d(log λ)` is maximized on the endemic branch (tracked via `I0 = 0.99`). Computed as:

```python
lam_mid = sqrt(lam_vals[:-1] * lam_vals[1:])   # geometric midpoints
dP      = diff(P_endemic) / diff(log(lam_vals))
lam_star = lam_mid[argmax(dP)]
```

### 6.4 Localization criterion (professor's definition)

> **Localized ↔ λ* > λ_c**

The system is localized if the sharpest response in the participation ratio occurs *above* (i.e., after) the invasion threshold. If λ* < λ_c, the endemic branch's sharpest transition lies in the bistable region (below the invasion threshold), and the system is not localized in this sense.

---

## 7. Sweep Parameters

```python
scale         = 5.0                        # fixed after cell-25097866 side-effect
lam_vals      = np.logspace(-6, -3, 50)    # 50 log-spaced points
r_max         = lam_vals[-1] * nmax**nu * 10  # ~3.3e10, required for ν=9.5
traj_points   = 500
t_max         = 200.0                      # long enough to reach steady state
I_cut         = 1e-6                       # mask P where I* too small
I0_vals       = [1e-15, 1e-5, 0.99]       # low IC (invasion), medium, high IC (persistence)
```

**Note on `r_max`:** With `ν = 9.5`, the endemic fixed point has `r* ~ λ · nmax^ν`. At `lam = 1e-3` and `nmax = 20`, `r_max ~ 20^9.5 · 1e-3 ≈ 3.3e10`. Earlier runs used `r_max = 10.0`, which missed all fixed points.

---

## 8. Initial Results (scale = 5.0)

```
Scalar  lam_c = 1.000e-06   lam* = 1.073e-06   localized = True
Kernel  lam_c = 3.556e-06   lam* = 1.073e-06   localized = False
Allegiance shifts lam* by factor 1.00x relative to scalar
```

### 8.1 Interpretation

| Quantity | Scalar | Kernel (scale=5) |
|---|---|---|
| λ_c | 1.000e-06 | 3.556e-06 (+3.6×) |
| λ* | 1.073e-06 | 1.073e-06 (unchanged) |
| Localized? | **Yes** (barely) | **No** |

**The allegiance kernel raises λ_c by ~3.5× without shifting λ*.** The mechanism:

- Allegiance drives nodes out of mixed groups (high φ(1-φ)), disrupting the early-spread dynamics that allow invasion from near-zero. This pushes λ_c upward.
- On the endemic (high-IC) branch, the sharpest change in group composition — and thus the sharpest change in P — occurs at roughly the same λ for both models, because the transition in group-size infection profiles is dominated by the synergy exponent ν, not the switching rate.
- Because λ* is unchanged but λ_c has moved up, the kernel's λ* now falls *below* λ_c — in the bistable window — so the system is classified as not localized.

**The two effects of the allegiance kernel pull in opposite directions:**
1. Segregation into all-infected groups → concentrates infection → *increases* localization
2. Disruption of mixed-group spreading → raises invasion barrier → *decreases* localization by the λ* > λ_c criterion

At scale = 5, effect (2) dominates.

### 8.2 Bistability

With ν = 9.5, the system exhibits a wide bistable window for both scalar and kernel models. The three initial conditions (`I0 = 1e-15`, `1e-5`, `0.99`) reveal:

- **Low IC (I0 = 1e-15):** decays to disease-free for λ < λ_c; invades only above λ_c
- **High IC (I0 = 0.99):** the endemic state persists down to much smaller λ (bistable region)
- **Medium IC (I0 = 1e-5):** behavior interpolates between the two

---

## 9. Open Questions

1. **Does stronger allegiance (scale >> 5) eventually make λ* > λ_c_kernel?**  
   At some scale the segregation effect may win: if all-infected clusters are self-sustaining, λ_c may stop rising while λ* shifts up with increasing group homogeneity.

2. **Quenched limit (scale = 0, w = 0):**  
   With no switching at all, does the localization pattern change fundamentally? This serves as a second baseline.

3. **Dependence on ν:**  
   With a lower synergy exponent (ν ~ 1–2), the fixed-point landscape is simpler (no bistability). Does the localization signal change qualitatively?

4. **nmax limitation:**  
   CNS has nmax = 20. Larger groups (higher nmax) would provide more resolution in P across group sizes and potentially a stronger localization signal.
