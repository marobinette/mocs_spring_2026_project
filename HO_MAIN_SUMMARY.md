# HO_MAIN_SUMMARY: Multistable Active Phases of Complex Contagions from Temporal Higher-Order Interactions

**Source:** Lamata-Otín, Malizia, Keating, St-Onge, Latora, Gómez-Gardeñes, Hébert-Dufresne (March 4, 2026)

---

## PASS 1 — Main Ideas

### Problem Statement

Prior higher-order contagion models showed that group-based reinforcement can produce discontinuous (explosive) transitions from an inactive to an active adoption state. However, those results were derived on **static hypergraphs** or under **mean-field (annealed) approximations** that ignore temporal dynamics — i.e., the fact that real groups gain and lose members over time. This paper asks: what happens when groups are temporally evolving?

### Core Contribution

The authors build an analytically tractable **Approximate Master Equations (AMEs)** model that couples higher-order contagion with temporal group turnover, parameterized by a single **plasticity rate ω**. The model continuously interpolates between:

- **Quenched limit** (ω = 0): fully static groups — no transition is possible because the system always collapses to the absorbing state
- **Annealed limit** (ω → ∞): infinitely fast reshuffling — reduces to standard mean-field predictions

### Three Main Results

**(i) Temporality reshapes phase transitions (homogeneous structures)**

- Increasing ω can convert a continuous transition into a discontinuous (first-order) one accompanied by bistability
- The invasion threshold λ_c(ω) becomes **non-monotonic** in ω when ν > 1 (nonlinear contagion): it first decreases with mixing (facilitating adoption), reaches an optimal minimum at ω*, then increases again (excessive reshuffling prevents reinforcement buildup)
- The bistable region grows wider as both ω and ν increase

**(ii) Four dynamical phases from heterogeneity + temporality**

When real-world group-size heterogeneity is incorporated, the system organizes into **four distinct dynamical regimes**:

1. **Continuous transition** — classic SIS-like behavior (weak reinforcement)
2. **Discontinuous transition** with absorbing–active bistability (intermediate reinforcement)
3. **Hybrid continuous transition with active bistability** — activity emerges continuously, but two stable active states coexist over a finite range
4. **Discontinuous hybrid transition** — three coexisting stable states: absorbing + two distinct active states

Regimes 3 and 4 (multistable active phases) are **absent in both the quenched and annealed limits** — they are intrinsically temporal phenomena that only emerge at finite plasticity rates. The active states differ in their degree of **mesoscopic localization**: lower active states are sustained by adopters concentrated in the largest groups.

**(iii) Real-world systems require stronger nonlinearity than previously thought**

Empirically measured temporal turnover rates and group-size distributions across seven real social datasets (high school, hospital, workplaces, village, conference, primary school) show that discontinuous transitions require substantially larger nonlinear reinforcement ν than static or mean-field models predict. Static aggregated representations suppress temporal features and can qualitatively mischaracterize both the nature and location of the critical transition.

### Key Conceptual Insight

Temporality is **not a secondary correction** to higher-order structure — it is a **primary control parameter** that determines the macroscopic transition class. The competition between contagion timescale (μ⁻¹) and plasticity timescale (ω⁻¹) governs whether transitions are continuous, discontinuous, or absent.

---

## PASS 2 — Parameters and Their Functionality

### Structural Parameters

| Parameter | Symbol | Role |
|-----------|--------|------|
| Group size | n | Size of an interaction group; governs reinforcement capacity |
| Group-size distribution | {p_n} | Probability distribution of group sizes across the system |
| Number of group memberships per individual | k | How many groups a node belongs to simultaneously |
| Membership distribution | {g_k} | Probability distribution over individual membership counts |

**Special cases used in the paper:**
- p_n = δ_{n,1}: all groups are size 1 (trivially no contagion)
- p_n = δ_{n,2}: pairwise interactions only (dyads)
- p_n = δ_{n,3}: triplets only (simplest higher-order case with nontrivial behavior)
- g_k = δ_{k,1}: each individual belongs to exactly one group (studied in Fig. 2)
- g_k = δ_{k,κ}: each individual belongs to exactly κ groups

### Dynamical Parameters

| Parameter | Symbol | Role |
|-----------|--------|------|
| Intrinsic adoption rate | λ | Per-unit rate at which a susceptible individual can be infected; the primary control parameter in phase diagrams |
| Recovery rate | μ | Rate at which adopters return to the susceptible state; sets the contagion timescale μ⁻¹ |
| Synergy exponent | ν | Controls nonlinearity of reinforcement; ν = 1 gives linear contagion, ν > 1 gives superlinear (synergistic) reinforcement |
| Plasticity rate | ω | Rate at which random individuals are swapped between random groups; ω = 0 is quenched, ω → ∞ is annealed; sets the plasticity timescale ω⁻¹ |

### State Variables

| Variable | Symbol | Definition |
|----------|--------|------------|
| Susceptible density with membership k | s_k | Fraction of individuals with k memberships who are susceptible at time t |
| Group-state distribution | f_{n,i} | Fraction of groups of size n that currently contain exactly i adopters |
| Global prevalence (order parameter) | I(t) | Fraction of the total population that are adopters at time t; I* denotes stationary value |
| Mean-field infection rate | r | Average contagion rate experienced by a randomly chosen susceptible group member |
| External infection pressure | ρ | Mean-field infection pressure from all groups a susceptible belongs to, *excluding* the focal group |

### Derived / Analytical Parameters

| Parameter | Symbol | Definition |
|-----------|--------|------------|
| Infection rate function | β(n,i) = λiᵛ | Rate at which a susceptible in a group of size n with i adopters becomes infected; the iᵛ factor encodes nonlinear reinforcement |
| Invasion threshold | λ_c | Critical spreading rate below which a small seed of adopters cannot spread; left boundary of the active phase |
| Persistence threshold | λ_p | Critical spreading rate below which an established active state cannot be sustained; left boundary of the bistable region |
| Optimal plasticity rate | ω* | Value of ω that minimizes the invasion threshold λ_c(ω) |
| Effective plasticity rate | ⟨ω⟩ | Empirical plasticity rate = inverse of average residence time of individuals in groups |
| Effective participation ratio | P | Measures how broadly activity is distributed across groups; P near p_{n_max} → localized; P near 1 → delocalized |
| Effective structural coupling | Q | Composite descriptor of inter-group connectivity integrating excess membership, group-size heterogeneity, and plasticity rate |
| Inter-event time | τ_e | Time between consecutive group-change events for the same individual |
| Functionals | F, G, H | Functions of {p_n}, β(n,i), ω appearing in the bistability threshold condition; defined in Supplementary Note 2 |

---

## PASS 3 — Equations Summary

### Eq. (1) — Quenched limit of invasion threshold

$$\lim_{\omega \to 0} \lambda_c(\omega) = \infty$$

**Meaning:** In the fully static (quenched) limit, the invasion threshold diverges — contagion can never spread from a small seed because adopters are trapped in finite groups with no new contacts. No active phase exists.

---

### Eq. (2) — Annealed limit of invasion threshold

$$\lim_{\omega \to \infty} \lambda_c(\omega) = \frac{\mu \langle n \rangle}{\langle k \rangle \langle n(n-1) \rangle}$$

**Meaning:** At infinitely fast reshuffling, individuals effectively sample all possible group configurations. This mean-field critical point holds for both linear (ν=1) and nonlinear (ν>1) contagion — nonlinearity has no effect on *where* the transition is, only on its *nature* (continuous vs. discontinuous).

---

### Eq. (3) — Invasion threshold for triplet structures

$$\lambda_{c,n=3}(\omega) = \frac{\mu + \omega}{2^\nu} \left[ \sqrt{1 + \frac{2^\nu \mu}{\omega}} - 1 \right]$$

**Meaning:** Exact closed-form invasion threshold when all groups are triplets (p_n = δ_{n,3}). A finite minimum exists only when ν > 2, explaining the onset of non-monotonic behavior.

---

### Eq. (4) — Optimal plasticity rate and minimal threshold (triplets)

$$\omega^*_{n=3} = \mu \frac{2^{\nu/2}}{2^{\nu/2} - 2}, \qquad \lambda^*_{c,n=3} = \mu \frac{2^{\nu/2} - 1}{2^{\nu-1}}$$

**Meaning:** The plasticity rate that minimizes the invasion threshold (ω*) and the corresponding minimal threshold value (λ*_c) for triplet structures. Only valid when ν > 2.

---

### Eq. (5) — Asymptotic persistence threshold (triplets, large ω)

$$\lim_{\omega \to \infty} \lambda_{p,n=3}(\omega) = \mu \frac{4(2^\nu - 2)}{2^{2\nu}}$$

**Meaning:** The persistence threshold (left boundary of bistability) for triplet structures in the annealed limit. Unlike λ_c, λ_p is strictly monotonic in ω and approaches this finite asymptote.

---

### Eq. (6) — Width of bistable region (triplets, large ω)

$$\lim_{\omega \to \infty} \Delta\lambda_{n=3}(\omega) = \frac{\mu}{2}\left(1 - 2^{2-\nu}\right)^2$$

**Meaning:** The width of the bistable region (Δλ = λ_c - λ_p) at high plasticity for triplet structures. Grows with ν, confirming that stronger nonlinearity produces wider bistability.

---

### Eq. (7) — Effective plasticity rate

$$\langle \omega \rangle = \frac{1}{\langle \tau \rangle}, \qquad \langle \tau \rangle = \frac{1}{E}\sum_{e=1}^{E} \tau_e$$

**Meaning:** Empirical plasticity rate defined as the inverse of the average inter-event time ⟨τ⟩, where τ_e is the time between consecutive group-change events for the same individual across E total events.

---

### Eq. (8) — AME for susceptible density

$$\frac{ds_k}{dt} = \mu(1 - s_k) - kr s_k$$

**Meaning:** Rate of change of the fraction of susceptibles with k memberships. Increases via recovery at rate μ (infected → susceptible); decreases via infection through any of the k groups at rate r per group.

---

### Eq. (9) — AME for group-state distribution

$$\frac{df_{n,i}}{dt} = (i+1)\bigl(\mu + \omega(1-I)\bigr) f_{n,i+1}$$
$$- \bigl[i(\mu + \omega(1-I)) + (n-i)(\beta(n,i) + \rho + \omega I)\bigr] f_{n,i}$$
$$+ (n-i+1)(\beta(n,i-1) + \rho + \omega I) f_{n,i-1}$$

**Meaning:** Three-term master equation for f_{n,i} (fraction of size-n groups with i adopters):
- **Term 1** (+): transitions i+1 → i via recovery (rate μ) or via a susceptible swapped in (rate ω, weighted by (1-I) = susceptible fraction in population)
- **Term 2** (−): all outflows from state (n,i) — infected recover (rate μ), infected swapped out (rate ω), susceptibles get infected internally (rate β), externally (rate ρ), or via swap (rate ωI)
- **Term 3** (+): transitions i-1 → i via internal infection (rate β(n,i-1)), external pressure (rate ρ), or swap of a susceptible with an infected individual (rate ωI)

---

### Eq. (10) — External infection pressure

$$\rho(r) = r \cdot \frac{\sum_k k(k-1)s_k g_k}{\sum_k k s_k g_k}$$

**Meaning:** The mean-field infection pressure from all groups a susceptible node belongs to *other than* the focal group. The fraction is the mean excess membership of a susceptible node (expected number of other groups it belongs to). Multiplied by r (the infection rate per external group contact).

---

### Eq. (11) — Mean-field infection rate

$$r = \frac{\sum_{n,i} \beta(n,i)(n-i) f_{n,i} p_n}{\sum_{n,i}(n-i) f_{n,i} p_n}$$

**Meaning:** The average infection rate β(n,i) experienced by a randomly chosen susceptible group member, averaged over the distribution of group states weighted by the number of susceptible members (n-i) in each group.

---

### Eq. (12) — Global prevalence (order parameter)

$$I(t) = \sum_k (1 - s_k(t)) g_k$$

**Meaning:** The fraction of the total population that are adopters, computed as the complement of the total susceptible fraction weighted by membership distribution. This is the primary observable tracking macroscopic contagion activity.

---

### Eq. (13) — Self-consistency functional

$$\mathcal{M}[\rho(r), I(r)] = \frac{\sum_{n,i} \beta(n,i)(n-i) f_{n,i}(\rho, I) p_n}{\sum_{n,i}(n-i) f_{n,i}(\rho, I) p_n}$$

**Meaning:** At stationarity, the entire system reduces to a functional of three scalar variables r, ρ, and I. The f_{n,i} are now the stationary group-state distributions parameterized by ρ and I.

---

### Eq. (14) — Fixed-point / self-consistency condition

$$r = \mathcal{M}[\rho(r), I(r)]$$

**Meaning:** At equilibrium, the mean-field infection rate r must equal its own self-consistent prediction M. Solutions to this implicit equation are the fixed points (stationary states) of the dynamics. Solved numerically; stability checked via the Jacobian.

---

### Eq. (15) — Invasion threshold condition

$$1 = \left\langle \sum_{i=1}^{n} \frac{n!}{(n-i-1)!\, i!} \left(\frac{\lambda_c}{\mu + \omega}\right)^i \prod_{j=1}^{i} j^\nu \right\rangle \times \left(\frac{\langle k(k-1) \rangle}{\langle k \rangle \langle n \rangle} + \frac{\omega}{\mu} \frac{\langle k \rangle}{\langle n \rangle}\right)$$

**Meaning:** Derived from the condition ∂M/∂r|_{r→0} = 1 (tangency at the trivial state). Determines the critical spreading rate λ_c at which a small seed can invade. The angle brackets denote averaging over group-size distribution {p_n}. The combinatorial factor counts reinforcement pathways; the second factor captures structural and temporal connectivity.

---

### Eq. (16) — Bistability (tricritical) threshold condition

$$0 = F\frac{\langle k(k-1)\rangle^2}{\langle k\rangle^2} + 2G\frac{\langle k(k-1)\rangle}{\langle k\rangle}\frac{\langle k\rangle}{\mu} + H\frac{\langle k\rangle^2}{\mu^2} + 2\frac{1}{\mu}\frac{\langle k^2\rangle^2 - \langle k^3\rangle\langle k\rangle - \frac{\omega}{\mu^2}\langle k^2\rangle\langle k\rangle^2}{\langle k(k-1)\rangle\langle k\rangle + \frac{\omega}{\mu}\langle k\rangle^3}$$

**Meaning:** Derived from ∂²M/∂r²|_{r→0} = 0 — the condition at which the invasion transition changes from second-order (continuous) to first-order (discontinuous). The tricritical line in (ω, λ) or (ν, λ) space is defined by the simultaneous satisfaction of Eqs. (14), (15), and (16). F, G, H are functionals of {p_n}, β(n,i), ω specified in Supplementary Note 2.

---

### Eq. (17) — General critical line for fixed k memberships

$$\lambda_c(\omega; k) = \frac{\mu + \omega}{2^\nu} \left[\sqrt{1 + \frac{2^\nu \mu}{k\omega + \mu(k-1)}} - 1\right]$$

**Meaning:** Generalization of Eq. (3) for individuals with exactly k group memberships (g_κ = δ_{κ,k}). Reduces to Eq. (3) when k=1 (single-group membership). Multiple memberships shift the minimum ω* toward zero.

---

### Eq. (18) — Optimal plasticity rate for k memberships

$$\omega^*(k) = \mu \frac{2^{\nu/2+1} + 2^{\nu+1} - 4 - k(2^\nu - 4)}{k(2^\nu - 4)}, \qquad (\nu > 2)$$

**Meaning:** Optimal plasticity rate (minimizer of invasion threshold) when each individual belongs to exactly k groups. As k grows, ω*(k) → 0 and eventually becomes non-positive, meaning the minimum disappears and λ_c(ω;k) becomes strictly increasing — plasticity only suppresses contagion.

---

### Eq. (19) — Minimal invasion threshold for k memberships

$$\lambda^*_c(k) = \frac{1}{k}\lambda^*_c(1) = \frac{\mu}{k}\frac{2^{\nu/2}-1}{2^{\nu-1}}$$

**Meaning:** The minimal achievable invasion threshold for individuals with k memberships. Scales as 1/k — more memberships lower the threshold by providing more pathways for contagion.

---

### Eqs. (20)–(23) — Persistence threshold (triplets, large ω)

$$\lambda_p(\omega) = \frac{-B_p(\omega) + \sqrt{B_p(\omega)^2 - 4A_p(\omega)C_p(\omega)}}{2A_p(\omega)}$$

where:
- **A_p(ω, μ) = 2^{ν+1}(2^ν ω + 3μ)**
- **B_p(ω, μ) = 4^ν ω² + 2^{ν+1} μω + 8μω + 9μ²**
- **C_p(ω, μ) = -4(2^ν - 2)(ω + μ)²**

**Meaning:** Closed-form expression for the persistence threshold (the left boundary of the bistable region — the value of λ below which even a fully adopted population collapses to the absorbing state) for triplet structures at large ω. Given as a root of a quadratic in λ with coefficients A_p, B_p, C_p that depend on ω, μ, and ν.

---

### Eq. (24) — Effective structural coupling

$$Q = \left[\frac{\langle k(k-1)\rangle}{\langle k\rangle} + \frac{\langle k\rangle\, \omega}{\mu}\right] \frac{\langle n(n-1)\rangle}{\langle n\rangle}$$

**Meaning:** A scalar descriptor that integrates structural heterogeneity and plasticity into a single measure of effective inter-group connectivity experienced during contagion. The bracket contains: (1) the structural excess membership ⟨k(k-1)⟩/⟨k⟩, and (2) a plasticity-driven mixing enhancement ⟨k⟩ω/μ. The multiplicative factor ⟨n(n-1)⟩/⟨n⟩ encodes excess group size (potential reinforcement interactions). Used to rank empirical systems: high-Q systems transition discontinuously at lower ν.

---

## PASS 4 — Dependencies Among Equations

The equations form a layered hierarchy. The diagram below organizes them from inputs through core dynamics to analytical characterizations.

```
STRUCTURAL INPUTS          DYNAMICAL INPUTS
{p_n}, {g_k}               λ, μ, ν, ω
      |                          |
      +----------+---------------+
                 |
          β(n,i) = λiᵛ    ← Infection function (not numbered, defined in text)
                 |
    ┌────────────────────────────────────┐
    │     CORE AME SYSTEM                │
    │   ds_k/dt  [Eq. 8]                │
    │   df_{n,i}/dt  [Eq. 9]            │
    │       ↑            ↑              │
    │   uses ρ [Eq. 10]  uses r [Eq. 11]│
    │       ↑            ↑              │
    │   uses r, {s_k}    uses {f_{n,i}},│
    │                    {p_n}, β        │
    └────────────────────────────────────┘
                 |
         I(t) [Eq. 12] ← uses {s_k, g_k}
                 |
    ┌─────────────────────────────────────┐
    │     STATIONARY ANALYSIS             │
    │   M functional [Eq. 13]            │
    │       ← uses r, ρ, I, {f_{n,i}}   │
    │   Fixed-point: r = M[ρ(r),I(r)]   │
    │   [Eq. 14] ← uses M, ρ, I         │
    └─────────────────────────────────────┘
                 |
     ┌───────────┴──────────────┐
     |                          |
∂M/∂r|_{r→0} = 1         ∂²M/∂r²|_{r→0} = 0
     |                          |
  Invasion threshold        Bistability threshold
  λ_c [Eq. 15]             (tricritical) [Eq. 16]
     |
  Special cases:
  Triplets [Eq. 3] → ω*, λ*_c [Eq. 4]
  General k [Eq. 17] → ω*(k) [Eq. 18] → λ*_c(k) [Eq. 19]
     |
  Persistence threshold λ_p:
  [Eqs. 20–23] ← A_p, B_p, C_p(ω, μ, ν)
     |
  Δλ = λ_c - λ_p [Eq. 6 for triplets]

Empirical systems:
  ⟨ω⟩ [Eq. 7] → locates system on tricritical line
  Q [Eq. 24] ← from Eq. 15 structure + {k}, {n}, ω, μ
```

### Detailed Dependency Chain

**β(n,i) → r and Eq. (9):** β(n,i) = λiᵛ is the fundamental transmission function. It appears directly inside the AME Eq. (9) and is the numerator weight in Eq. (11) for r.

**r ↔ {f_{n,i}} and {s_k} (mutual):** r [Eq. 11] is computed from f_{n,i} and p_n; but r also drives ds_k/dt [Eq. 8]. Similarly, ρ [Eq. 10] is computed from s_k and g_k, and also appears in Eq. (9). This creates a **closed coupled system**: (s_k, f_{n,i}) → (r, ρ) → back into the ODEs.

**I(t) [Eq. 12]:** Derived purely from {s_k} and {g_k}. It then feeds back into Eq. (9) (the ω(1-I) and ωI swap terms) — so I is not just an output but also a dynamical coupling variable.

**Stationary reduction [Eqs. 13–14]:** At steady state, the full ODE system collapses to the scalar self-consistency equation r = M[ρ(r), I(r)]. The functional M has the same form as r in Eq. (11), but with stationary f_{n,i}(ρ, I). This is the critical simplification that makes the model analytically tractable.

**Invasion threshold [Eq. 15] ← Eq. 14 linearized:** Obtained by differentiating Eq. (14) with respect to r and evaluating at r → 0 (where the trivial absorbing state lives). The condition ∂M/∂r|_{r→0} = 1 defines the boundary of stability of the disease-free state. Contains λ_c, μ, ω, ν through the combinatorial reinforcement factor, and {k}, {n} moments through the connectivity factor.

**Bistability threshold [Eq. 16] ← Eq. 14, second-order:** Obtained from the second derivative condition ∂²M/∂r²|_{r→0} = 0. When combined with Eq. (15), it identifies the tricritical point in (λ, ν, ω) space where the nature of the transition switches from continuous to discontinuous.

**Equations (3), (17) are specializations of Eq. (15):** For p_n = δ_{n,3} and g_k = δ_{k,1}, Eq. (15) closes to Eq. (3). For general δ_{k,κ}, it gives Eq. (17). These analytic closures are only possible for monodisperse distributions.

**Equations (4), (18), (19) from differentiating (3) and (17):** The optimal ω* comes from dλ_c/dω = 0 applied to Eqs. (3) and (17). Evaluating those expressions at ω* gives the minimal thresholds Eqs. (4) and (19).

**Persistence threshold [Eqs. 20–23] ← complementary to invasion threshold:** While λ_c comes from linearizing around the disease-free state (r → 0), λ_p comes from analyzing the stability of the fully active state (r large). The quadratic form in Eqs. (20–23) is derived via a separate perturbative analysis (Supplementary Note 5) specific to triplet structures at large ω.

**Eq. (6) ← Eqs. (2) and (5):** The bistable region width Δλ = λ_c - λ_p at large ω is obtained by substituting the annealed invasion threshold Eq. (2) and the asymptotic persistence threshold Eq. (5) for triplets.

**Effective structural coupling Q [Eq. 24] ← inspection of Eq. 15:** By examining the factored structure of the invasion condition, the authors identify that {k} moments, {n} moments, and ω enter through a single grouped combination, motivating the definition of Q as a composite scalar.

**Effective plasticity rate ⟨ω⟩ [Eq. 7]:** Empirically measured from real-data inter-event times. Serves as input into all analytical expressions (Eqs. 15–16, 24) to locate real systems on the theoretical phase diagram.

---

## Summary Table of Equation Purposes

| Equation | Type | Purpose |
|----------|------|---------|
| β(n,i) = λiᵛ | Definition | Nonlinear infection rate function |
| (1) | Limiting case | Invasion threshold diverges at ω=0 |
| (2) | Limiting case | Mean-field invasion threshold at ω→∞ |
| (3) | Special case | Closed-form threshold for triplet groups |
| (4) | Optimization | Optimal ω and minimal threshold for triplets |
| (5) | Asymptotic | Persistence threshold at large ω for triplets |
| (6) | Derived | Bistable region width at large ω for triplets |
| (7) | Measurement | Empirical plasticity rate from data |
| (8) | ODE | AME for susceptible population with k memberships |
| (9) | ODE | AME for group-state distribution f_{n,i} |
| (10) | Closure | External infection pressure (mean-field term) |
| (11) | Closure | Mean-field infection rate r |
| (12) | Observable | Global prevalence I(t) — order parameter |
| (13) | Stationary | Self-consistency functional M at equilibrium |
| (14) | Fixed-point | Self-consistency equation (implicit); solved numerically |
| (15) | Critical | Invasion threshold condition (1st-order perturbation of Eq. 14) |
| (16) | Critical | Bistability/tricritical condition (2nd-order perturbation of Eq. 14) |
| (17) | Generalization | Critical line for arbitrary k memberships |
| (18) | Optimization | Optimal ω*(k) for k memberships |
| (19) | Optimization | Minimal threshold λ*_c(k) |
| (20)–(23) | Analytic | Persistence threshold λ_p(ω) for triplets |
| (24) | Composite | Effective structural coupling Q for ranking systems |

---

## Quick Reference: Key Physics

| Regime | ω | ν | Transition Type |
|--------|---|---|-----------------|
| Quenched | 0 | any | None (λ_c → ∞) |
| Low temporality, linear | low | 1 | Continuous |
| High temporality, linear | high | 1 | Continuous (threshold decreases monotonically) |
| Low temporality, nonlinear | low | >1 | Continuous (threshold high) |
| Intermediate temporality, nonlinear | ω* | >2 | Continuous (threshold minimized at ω*) |
| High temporality, nonlinear | high | >1 | Discontinuous + bistability |
| Any temporality, nonlinear + heterogeneous {p_n} | finite | >1 | Possibly: hybrid continuous with active bistability OR discontinuous with 3 coexisting states |
| Annealed | ∞ | any | Discontinuous (no multistable active phases) |
