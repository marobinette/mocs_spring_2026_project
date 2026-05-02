# Presentation Outline — 7 Minutes
## Diversity-Tension Kernel in Temporal Higher-Order Contagion

---

## 1. Motivation: What the Original Paper Established (~1.5 min)

- **Higher-order contagion on temporal networks:** Real groups gain and lose members over time. The paper (Lamata-Otín et al. 2026) showed that ignoring this — as static hypergraph models do — qualitatively mischaracterizes both the nature and location of critical transitions.

- **The core finding:** A single constant plasticity rate ω governs everything. Increasing ω can convert a *continuous* (smooth) transition into a *discontinuous* (explosive) one with bistability — two stable states coexisting for the same parameters.

- **Four dynamical phases emerge** when real group-size heterogeneity is added: continuous, discontinuous (absorbing-active bistability), hybrid continuous with active bistability, and a three-state regime (absorbing + two active states). The richer phases only exist at *finite* ω — not in the static or mean-field limits.

- **The gap this creates:** ω is a fixed external parameter. It doesn't respond to what is actually happening inside the groups at any moment. Every group switches at the same rate regardless of whether it is homogeneous or deeply divided. This is socially unrealistic.

- **Our motivating question:** What if the rate at which people leave a group depends on the *diversity* of that group? Groups where everyone agrees (homogeneous) should be stable; groups caught between two factions (mixed) should be volatile. This is the idea of **diversity tension**.

---

## 2. The Kernel — Contagion Dynamics Shape Group Dynamics (~1 min)

- **The causal arrow in the original paper:** Group dynamics (plasticity rate ω) → shape → contagion dynamics (phase portrait). ω is set from outside; the system is driven.

- **Our reversal:** Contagion dynamics (local prevalence φ = i/n inside each group) → shape → group switching rate w(n, i) → feeds back into → contagion dynamics. The system is *coevolutionary*.

- **The diversity-tension kernel:**
  - w(n, i) = α · φ · (1 − φ), where φ = i/n is the infected fraction within a group
  - φ = 0 or 1 (homogeneous groups): **w = 0** — no switching, people stay
  - φ = 0.5 (maximally mixed groups): **w = α/4** — maximum switching, people flee
  - α controls the overall amplitude of the effect

- **The social mechanism (Laurent's framing):** *"If the group is 50-50, people have to do too much code switching. Whereas if you're the only one in the minority, you just use the code of the other people."* Allegiance and homophily drive people out of divided groups and into groups that match their own state.

- **Dynamical consequence — self-quenching at endpoints:** The kernel vanishes wherever groups are homogeneous. This means the system *self-regulates*: once disease dynamics push groups toward all-infected or all-susceptible, the switching rate drops to zero and those states become frozen. The system quenches itself.

---

## 3. Method — Leveraging the wAME Framework (~1 min)

- **Foundation:** We keep the paper's full wAME (approximate master equations) structure intact. The two state variables are the same: s_k (susceptible fraction among nodes in k simultaneous groups) and f_{n,i} (fraction of size-n groups with i infected members).

- **One structural modification — replace ω with w(n, i):** In the paper, every node switches at the same rate ω. In our model, the switching rate depends on which (n, i) state the group is in. This means the pool of people who are actively switching is no longer a random cross-section of the population.

- **A new derived quantity — S_w:** Because switchers come disproportionately from mixed groups, the probability that a randomly chosen switcher is *susceptible* is no longer simply 1 − I (the global susceptible fraction). We must track a flux-weighted susceptible fraction:

  > S_w = Σ (n−i) · w(n,i) · f_{n,i} · p_n  /  Σ n · w(n,i) · f_{n,i} · p_n

  When w is constant, S_w reduces exactly to 1 − I (recovering the paper). When w is the diversity-tension kernel, S_w is biased toward ≈ 0.5 because only mixed groups are contributing switchers.

- **Network — Synthetic_delta_k5:** We use a synthetic network where every individual belongs to exactly 5 groups simultaneously (g_k = δ_{k,5}). This activates the full multi-membership structure of the model (ρ > 0, the shared-membership infection pathway), while keeping the network regular and interpretable — avoiding the artifacts in empirical datasets like CNS.

---

## 4. The Tricritical Line Problem — Why We Sweep (~1 min)

- **What the paper did analytically:** The tricritical line (Eq. 16) separates continuous from discontinuous transitions in (ω, λ) space. It is derived by setting d²M/dr² = 0 at threshold. Because ω is a scalar constant, it factors cleanly out of every sum in the AMEs — the algebra closes.

- **Why this fails with the kernel:** w(n, i) varies independently for every (n, i) pair. It cannot be factored out of the sums. Every derivative picks up an additional (n, i) dependence that prevents the condition from closing analytically. As Laurent put it: *"This is like a big physics project — maybe a summer thing, or for the master's student in Spain."*

- **The numerical alternative — vacc_sweep.py:**
  - Sweep a 50 × 40 grid in (λ, ν) space: λ from 10⁻⁵ to 10⁻¹ (log scale), ν from 1.0 to 12.0
  - For each (λ, ν) cell, integrate the full wAME system to steady state from **two initial conditions**: I₀ = 0.001 (invasion) and I₀ = 0.99 (persistence)
  - Record δI = I*(high) − I*(low) — nonzero δI means two distinct attractors exist (bistability)
  - Run one job per α value on the VACC cluster; baseline (constant ω = 5) run for comparison

- **The δI heatmap is the computational tricritical map:** The boundary of the δI > 0 region in (λ, ν) space *is* the tricritical line evaluated at fixed ω (or fixed α), without requiring the analytic derivation.

---

## 5. Early Results — How They Connect Back to the Tricritical Line (~1.5 min)

- **Baseline confirms the framework:** The constant-ω = 5 baseline produces 366 bistable cells with max δI = 0.93. The δI > 0 region concentrates at low λ and high ν — exactly the "Discontinuous" upper-left regime from Figure 3a of the paper. The method is recovering the right physics.

- **The kernel systematically suppresses bistability:**

  | α | bistable cells (δI > 0.05) | peak kernel rate |
  |---|---------------------------|-----------------|
  | 0.1 | 271 | 0.025 (≪ μ) |
  | 1.0 | 106 | 0.25 (< μ) |
  | 2.0 | 33 | 0.50 (< μ) |
  | **5.0** | **0** | **1.25 (> μ)** |
  | 10.0 | 0 | 2.50 |
  | 20.0 | 0 | 5.00 |

- **A natural critical scale:** Bistability collapses completely between α = 2 and α = 5 — exactly when the peak kernel rate (α/4) crosses the recovery rate μ = 1. Once mixing is faster than recovery, the system cannot sustain the heterogeneous intermediate states needed to support two attractors.

- **Mechanistic interpretation in terms of Figure 3a:** The paper's tricritical line separates "Continuous" (lower-left) from "Discontinuous" (upper-right) in (ω, λ) space. The kernel's self-quenching at the endpoints (w → 0 at φ = 0, 1) is equivalent to pulling the *effective* ω toward zero whenever the system approaches the homogeneous fixed points. Above the critical α, this pull is strong enough that the system can never remain in the "Discontinuous" regime — bistability is structurally eliminated.

- **What the plots show directly:**
  - **Heatmap grid:** The bistable region (bright plasma colors) shrinks and shifts to ever lower λ as α grows, then vanishes entirely above α = 5
  - **Bistable-cell count plot:** A sharp drop between α = 2 and α = 5, well below the baseline's 366-cell reference line
  - **Overlay contour:** Each α's bistable footprint contracts inward; the baseline (gray) sets the outer boundary that the kernel progressively erodes
