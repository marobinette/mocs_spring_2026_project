# VACC Parameter Sweep — Full Summary
**Date:** 2026-04-29

---

## Background

The local 15×15 sweep revealed a clear and interpretable result: the diversity-tension kernel (α=30) completely eliminates bistability relative to the constant-ω baseline. The baseline recovered the expected phase structure from the paper — bistability appearing for ν ≥ ~5, forming a diagonal staircase across (λ, ν) space. The kernel panel was uniformly blank. δI = 0 everywhere.

This is not a failure — it is a finding. The kernel appears to convert the discontinuous (bistable) transition to a continuous one. The mechanism is the self-quenching at the homogeneous endpoints: `w(n,i) = α·φ·(1-φ)` goes to zero when φ=0 (fully susceptible groups) and φ=1 (fully infected groups). This freezes group dynamics precisely at the states the two initial conditions start from, preventing the system from sustaining the coexisting branches that produce bistability.

The immediate question this raised: is this α-dependent? At α=30, max ω=7.5. As α decreases toward zero, the kernel should weaken and eventually the system should recover bistability. There must be a critical α below which bistability re-emerges.

---

## The VACC Sweep Strategy

**Grid:** 50 λ values × 40 ν values = 2,000 cells per run. λ runs log-spaced from 1e-5 to 1e-1 (extended left from the local grid's 1e-4, to catch transitions that may sit at lower transmission rates). ν runs linear from 1.0 to 12.0.

**Two initial conditions per cell (Laurent's two-IC strategy):**
- `I0_low = 1e-3` — invasion attempt. System starts near the absorbing state and climbs to the first stable active branch from below, if one exists.
- `I0_high = 0.99` — persistence attempt. System starts near full infection and relaxes to the first stable active branch from above.

For each cell, both trajectories are integrated to stationarity (T_MAX=200, which is 200 characteristic recovery times with μ=1) and the final infected fraction I* is recorded. The bistability indicator is:

**δI = I\*(high) − I\*(low)**

This is large wherever two stable branches coexist and collapses to zero where there is only one stable fixed point. It is the computational substitute for the analytical tricritical line, which cannot be derived analytically for the kernel model (see `param_sweep_guidance.md`, §4).

**Ten runs total:**
- 9 kernel runs: α ∈ {0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 50.0}
- 1 constant-ω baseline: ω=5.0 (paper-style, no kernel)

The α values were chosen to bracket the regime where bistability should re-emerge. α=30 kills it. α=0.1 is nearly quenched (w≈0 everywhere). The interesting transition should occur somewhere in the low-to-mid range, likely between α=1 and α=10.

---

## What We Are Looking For

**Primary question:** At what α does bistability re-emerge, and how does it re-emerge?

There are three possible outcomes for each kernel panel:

1. **δI = 0 everywhere** — the kernel fully suppresses bistability at this α. The transition is continuous.
2. **δI > 0 in some cells, shifted relative to baseline** — bistability exists but the bistable region has moved in (λ, ν) space. Shifted right (higher λ) means the kernel raises the invasion threshold. Shifted left means it lowers it.
3. **δI > 0 in the same region as baseline** — the kernel has no effect on the phase structure at this α.

By comparing all 9 kernel panels against the baseline, we can map out:
- The **critical α** below which bistability first appears in the kernel model
- Whether bistability re-emergence is sharp (sudden appearance) or gradual (growing δI signal)
- Whether the bistable region shifts in (λ, ν) space as α changes — this would tell us whether the kernel primarily affects the invasion threshold, the persistence threshold, or both

**Secondary question:** Does the kernel produce any qualitatively new behavior not seen in the baseline — for example, bistability in a region where constant-ω shows none, or a different shape to the bistable boundary?

---

## VACC Job Details

| File | Job ID | Status |
|---|---|---|
| `sweep_kernel_alpha_0.5.npz` | 4040435_1 | ✅ Complete |
| `sweep_kernel_alpha_2.0.npz` | 4040435_3 | ✅ Complete |
| `sweep_kernel_alpha_5.0.npz` | 4040435_4 | ✅ Complete |
| `sweep_kernel_alpha_10.0.npz` | 4040435_5 | ✅ Complete |
| `sweep_kernel_alpha_20.0.npz` | 4040435_6 | ✅ Complete |
| `sweep_kernel_alpha_30.0.npz` | 4040435_7 | ✅ Complete |
| `sweep_kernel_alpha_50.0.npz` | 4040435_8 | ✅ Complete |
| `sweep_kernel_alpha_0.1.npz` | 4041880_0 | ⏳ Running |
| `sweep_kernel_alpha_1.0.npz` | 4041880_2 | ⏳ Running |
| `sweep_baseline.npz` | 4041881 | ⏳ Running |

First 7 jobs hit the original 6-hour time limit (~85% complete). Resubmitted with 12-hour limit.

---

## Next Steps

1. **Wait for the 3 remaining files** — currently running on VACC (jobs 4041880, 4041881).

2. **Pull results to local machine** — once all 10 files are in `Files/vacc/`, run from local terminal:
   ```bash
   scp -r mrobine1@vacc-login.uvm.edu:/gpfs1/home/m/r/mrobine1/mocs_spring_2026_project/Files/vacc/ ~/mocs_spring_2026_project/Files/
   ```

3. **Build the analysis figure** — a grid of heatmaps, one panel per α value plus the baseline, all on the same color scale. The evolution of the δI pattern across α values will show exactly where and how bistability re-emerges.

4. **Identify the critical α** — find the lowest α at which δI > 0 appears. This is the threshold below which the diversity-tension kernel is weak enough that bistability survives.

5. **Draft the abstract** — once the critical α is identified and the direction of any shift in the bistable region is clear, you have a concrete finding: the kernel suppresses bistability above a critical amplitude, converting a discontinuous transition to a continuous one. That is the central result.

6. **Decide whether a follow-up sweep is needed** — if the critical α falls between two of the sampled values (e.g., between α=5 and α=10), a finer sweep around that range would sharpen the result. This can be run on VACC with the same script using `--alpha` directly.
