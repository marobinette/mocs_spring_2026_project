# Findings: Diversity-Tension Kernel and Bistability

**Date:** 2026-05-03  
**Networks analyzed:** `Synthetic_poisson_k5`, `Thiers13`  
**Sweep:** α ∈ {0.1, 0.18, 0.32, 0.56, 1.0, 1.78, 3.16, 5.62, 10, 17.8, 31.6, 56.2, 100} (log-spaced), λ ∈ [10⁻⁴, 10⁻¹] (50 pts), ν ∈ [1, 12] (40 pts)  
**Bistability indicator:** δI = I\*(high I₀) − I\*(low I₀); a cell is "bistable" if δI > 0.05  
**Ridgeline check:** Re-integrated peak-δI cells to T_MAX=300 (15× longer) to check stationarity

---

## Critical caveat: much of the sweep signal may be pseudo-bistability

Before discussing the α-sweep results, the ridgeline analysis raises a red flag that must be addressed first.

We re-ran the ODE at the peak-δI parameter point for each network (using T_MAX=300 vs. the sweep's T_MAX=20) and compared where low-I₀ and high-I₀ trajectories end up. In both cases, **both initial conditions converged to the same I\***:

| Network              | α (bistable) | I₀=10⁻³ → I\* | I₀=0.99 → I\* |
|----------------------|-------------|---------------|---------------|
| Thiers13             | 3.16        | 0.248         | 0.248         |
| Synthetic_poisson_k5 | 0.1         | 0.578         | 0.578         |

These peak-δI cells have only **one attractor**. The large δI values in the sweep (up to 0.569 and 0.190 respectively) were artifacts of the sweep's short integration horizon: at T=20, the low-I₀ run had not yet reached the endemic fixed point while the high-I₀ run had. The two initial conditions were at different transient stages of approaching the same state — **not in different basins of attraction**.

This means the δI signal in the sweep is at least partially, and possibly predominantly, a measure of **differential convergence speed** rather than the existence of multiple stable states. The "bistable cells" counted by the threshold δI > 0.05 mix genuinely bistable points with slow-converging monostable ones.

This does not make the sweep results meaningless, but it demands a different interpretation. The α-dependence of the bistable cell count tells us how α modulates the convergence timescale and the character of the approach to stationarity — which is interesting in its own right — but claims about multiple attractors cannot be made from the current data without longer integration times.

A further hint: for the "flat" α=100 run in Thiers13, the ridgeline shows a small residual difference (I\*(low)=0.078 vs I\*(high)=0.094 at T_MAX=300). Whether this reflects true bistability with shallow basins or a very slow transient is unclear and warrants closer investigation at that parameter point.

---

## 1. What the α-sweep does tell us

With the pseudo-bistability caveat in mind, the sweep tables still show a consistent and interpretable pattern.

**Synthetic_poisson_k5** — δI signal decreases *monotonically* with α:

| α      | bistable cells | max δI |
|--------|---------------|--------|
| 0.1    | 140           | 0.569  |
| 0.32   | 66            | 0.214  |
| 1.0    | 8             | 0.069  |
| 1.78+  | 0             | <0.04  |
| baseline (ω=5 const.) | **301** | **0.882** |

**Thiers13** — δI signal is *non-monotonic*, peaking at intermediate α:

| α      | bistable cells | max δI |
|--------|---------------|--------|
| 0.1    | 413           | 0.188  |
| 1.0    | 480           | 0.197  |
| 3.16   | **527**       | 0.190  |
| 10     | 515           | 0.181  |
| 31.6   | 383           | 0.175  |
| 100    | 230           | 0.137  |
| baseline (ω=5 const.) | 175 | 0.540 |

Even if these numbers measure convergence speed rather than true bistability, the direction of the effect is real: the diversity-tension kernel accelerates convergence from low I₀ in Synthetic_poisson_k5 (shrinking the apparent bistable region) while it slows convergence from both I₀ conditions in Thiers13 in a way that expands it. The question of what this means for the true attractor structure requires the longer sweeps described in Section 5.

---

## 2. The structural key: membership, not group size

Both networks have the same `nmax = 6` and nearly identical mean group size (≈ 2.31). Their group-size distributions are essentially the same. The dramatic difference lies in the **membership distribution**:

| Network              | nmax | mean\_n | mmax | mean\_k |
|----------------------|------|---------|------|---------|
| Synthetic_poisson_k5 | 6    | 2.31    | 25   | 5.00    |
| Thiers13             | 6    | 2.31    | 5    | 1.04    |

In Synthetic_poisson_k5, each node simultaneously belongs to ~5 groups. In Thiers13, the average node belongs to barely more than one group.

**High membership (Synthetic_poisson_k5, mean_k = 5):** A node's infection state is shaped by a superposition of five simultaneous group memberships. Switching out of one mixed group provides no escape from infection pressure from the other four. The diversity-tension kernel's selective switching in mixed groups is therefore diluted — no single group's composition strongly determines a node's trajectory. The constant-ω baseline at ω=5 already generates substantial turnover; the kernel makes switching faster in exactly the mixed groups that would otherwise mediate slow convergence from low I₀, which accelerates the approach to the endemic state and shrinks the apparent bistable region.

**Low membership (Thiers13, mean_k ≈ 1):** A node is tightly coupled to a single group's composition. The kernel's sorting pressure at moderate α is felt acutely, pushing groups toward homogeneous compositions (all-S or all-I). This slows convergence from *both* initial conditions (because homogeneous groups are near fixed points where infection cannot spread or die out quickly), which makes the transient window longer and expands the apparent bistable region relative to the constant-ω baseline.

---

## 3. The genuinely new finding: size-stratified group structure

The most mechanistically revealing result from the ridgeline analysis is **not** related to bistability at all. At α=100 in Synthetic_poisson_k5, the stationary fₙᵢ distribution shows a striking **size-stratification** absent at α=0.1:

- **Small groups (n=2, 3):** fₙᵢ concentrated near φ=0 (almost entirely susceptible)  
- **Large groups (n=5, 6):** fₙᵢ concentrated near φ=1 (almost entirely infected)  
- Overall I\* ≈ 0.157 — a relatively low endemic state with a highly heterogeneous group structure

At α=0.1, all group sizes have fₙᵢ concentrated near φ=1, and I\* ≈ 0.578 — high infection with no size-dependent differentiation.

This stratification is a direct mechanical consequence of the kernel. The switching rate w(n, i, α) = 4α·(i/n)·(1−i/n) vanishes at φ=0 and φ=1 and peaks at φ=0.5. At high α:

- A group of size 2 with one infected node faces switching rate w(2, 1, 100) = 100. Mixed size-2 groups are eliminated almost instantly, driving them to the all-susceptible absorbing state (φ=0), where w=0 and they remain.
- A group of size 6 with 3 infected nodes faces w(6, 3, 100) = 100 as well, but the all-infected state (φ=1) is also absorbing (w=0). Large groups are "trapped" in their all-infected configuration once they reach it.

The kernel therefore creates a **size-dependent trapping mechanism**: large groups accumulate infection and lock in, small groups purge it and lock out. This produces a lower global I\* than at α=0.1 (the infection is concentrated in a smaller share of the population, the members of large groups), but it is a *single stable state* with heterogeneous internal structure — not bistability.

For Thiers13, this size-stratification is weaker. All group sizes remain concentrated near φ=1 at both α values. With mean_k ≈ 1, nodes rarely switch at all, so the kernel's trapping effect on small groups cannot operate — there is not enough switching flux to drive small groups to all-susceptible states.

---

## 4. Revised picture

| What we thought | What the evidence now suggests |
|---|---|
| Higher α → more bistability | Higher α changes convergence speed and group structure, but direction depends on network |
| Kernel promotes bistability in Thiers13 | Kernel expands the *apparent* bistable region but may be extending transients, not creating new attractors |
| Group-size distributions differ between networks | Group-size distributions are essentially identical; the difference is entirely in **membership** (mean_k) |
| The bistable cells in the sweep reflect multiple stable states | At least the peak-δI cells are monostable with slow convergence; need longer T_MAX to assess the rest |
| Ridgelines would show different fₙᵢ between I₀ conditions | Ridgelines show the same fₙᵢ (same attractor), but reveal a novel size-stratification at high α |

---

## 5. Next steps

### 5.1 Re-run the sweep with much longer T_MAX (highest priority)

The clearest priority is to repeat the sweep at T_MAX=500–1000 and verify whether any cells in the (λ, ν) grid show true bistability (I\* genuinely different between low and high I₀ at stationarity). This separates the real bistable region from pseudo-bistability. The VACC is the right place to do this given the compute cost. Start with a single α (say the sweep-detected "best" α per network) and a coarser (λ, ν) grid.

### 5.2 Investigate the Thiers13 α=100 residual δI

At α=100 and T_MAX=300, Thiers13 still shows I\*(low)=0.078 vs I\*(high)=0.094 at the tested parameter point. Extend T_MAX to 1000–2000 at this specific (λ, ν) point to determine whether the difference persists (true bistability) or vanishes (very slow transient). If it persists, map out where in (λ, ν) space this true bistability exists.

### 5.3 Characterize the size-stratification systematically

The size-stratified attractor seen in Synthetic_poisson_k5 at high α is interesting independent of bistability. Run a sweep over α at fixed (λ, ν) and measure the *spread* in mean φ across group sizes (a size-stratification index). Map when and where stratification appears, and whether it correlates with I\* suppression. This may be the most novel finding from the kernel.

### 5.4 Isolate the membership effect with controlled synthetics

Run the α sweep on synthetic networks where only mean_k is varied (k=1, 2, 3, 5, 10) while holding the group-size distribution fixed. This would test whether the transition between the two behavioral regimes (monotone suppression vs. non-monotone with a peak) is predicted by mean_k alone, and if so, where the crossover occurs.

### 5.5 Match the baseline switching rate to the kernel mean

The constant-ω baseline uses ω=5, but the diversity-tension kernel at a given α produces a mean switching rate that varies with the population state (it is zero in all-homogeneous states). At α=3.16 and the Thiers13 stationary state, the effective mean switching rate may be much less than 5. Computing this effective rate and rerunning the baseline at the matched value would give a cleaner comparison of structure (kernel vs. constant) vs. rate.

### 5.6 Analytical invasion threshold

In the low-I limit, linearizing around I=0 gives an invasion eigenvalue that determines whether a small outbreak can grow. The diversity-tension kernel enters through S_w, which depends on fₙᵢ and the kernel values. Near I=0, fₙᵢ is concentrated at i=0 (all groups mostly susceptible), so S_w → 1 and the kernel vanishes — meaning the kernel has **no effect on the invasion threshold** at leading order. This means any kernel-induced bistability (if real) cannot come from shifting the lower threshold; it must come from modifying the upper endemic attractor or the separatrix between basins. Verifying this analytically would sharpen the mechanistic story considerably.
