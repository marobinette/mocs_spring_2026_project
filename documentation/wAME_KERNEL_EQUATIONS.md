# wAME Kernel Equations

Equations implemented in `wAME.py` (`vector_field_w_kernel`) compared against the paper's AMEs (Eqs. 8–12). The kernel version replaces the scalar plasticity rate ω with a state-dependent switching function w(n, i) and introduces a new derived quantity S_w.

---

## Kernel Function

The diversity-tension kernel implemented in `w_diversity_tension`:

$$w(n, i) = \alpha \cdot \frac{i}{n} \cdot \left(1 - \frac{i}{n}\right)$$

where α = `scale` is the overall switching amplitude. Letting φ = i/n denote the infected fraction within a group:

- φ = 0 or φ = 1 (homogeneous groups): w = 0 — nodes **do not switch** out of groups where everyone shares the same state
- φ = 0.5 (maximally mixed groups): w = α/4 (maximum) — nodes **flee** diverse groups

This encodes **allegiance / homophily**: individuals prefer to stay in groups whose composition matches their own state.

---

## Kernel-Weighted Susceptible Fraction S_w

$$S_w = \frac{\displaystyle\sum_{n,i}(n - i)\, w(n,i)\, f_{n,i}}{\displaystyle\sum_{n,i} n\, w(n,i)\, f_{n,i}}$$

**Physical meaning:** S_w is the probability that a node entering a group (via switching) is susceptible. It is a mean-field average over all group states — the same type of calculation as r — built from first principles:

- **Numerator:** rate at which susceptible nodes leave their groups, summed over all group states. A group in state (n, i) has (n−i) susceptibles each leaving at rate w(n,i), weighted by f_{n,i}.
- **Denominator:** rate at which any node leaves its group, summed over all group states. Same group has n total nodes each leaving at rate w(n,i).

Since w(n,i) is a **group trait** (all nodes in a group switch at the same rate regardless of state), groups are fully distinguished by (n, i) alone. The sum is only over n and i — no p_n needed.

| Case | S_w |
|------|-----|
| w(n,i) = ω (constant, paper's model) | 1 − I (global susceptible fraction) |
| Diversity-tension kernel | Biased toward ~0.5 because switchers come from mixed groups |

**Fallback:** if the total switching flux is zero (degenerate case — e.g., all groups are homogeneous), S_w = 1 − I.

---

## Dynamical Equations

### Susceptible population — unchanged from paper Eq. (8)

$$\frac{ds_k}{dt} = \mu(1 - s_k) - k\, r\, s_k$$

This equation is identical to the paper. The kernel does not affect the individual-level susceptible equation directly; its effect enters only through r.

---

### Group-state distribution — kernel version of paper Eq. (9)

$$\frac{df_{n,i}}{dt} =
(i+1)\bigl(\mu + w(n,\, i+1)\cdot S_w\bigr)\, f_{n,i+1}$$

$$-\Bigl[i\bigl(\mu + w(n,i)\cdot S_w\bigr)
+ (n-i)\bigl(\beta(n,i) + \rho + w(n,i)\cdot(1 - S_w)\bigr)\Bigr]\, f_{n,i}$$

$$+ (n-i+1)\bigl(\beta(n,i-1) + \rho + w(n,i-1)\cdot(1 - S_w)\bigr)\, f_{n,i-1}$$

**Term 1 — transitions i+1 → i** (infected count decreases by 1):
An infected individual either recovers (rate μ) or is swapped out of the group and replaced by a susceptible (rate w(n, i+1) · S_w). The factor S_w reflects that the incoming individual is susceptible with probability S_w.

**Term 2 — outflows from state (n, i)**:
- Each of the i infected individuals recovers (rate μ) or is swapped out and replaced by a susceptible (rate w(n,i) · S_w)
- Each of the n − i susceptible individuals becomes infected internally (rate β(n,i)), externally (rate ρ), or is swapped out and replaced by an infected individual (rate w(n,i) · (1 − S_w))

**Term 3 — transitions i−1 → i** (infected count increases by 1):
One of the n − i + 1 susceptible individuals in a group with i − 1 infected becomes infected via internal transmission (rate β(n, i−1)), external pressure (rate ρ), or by being swapped out and replaced by an infected individual (rate w(n, i−1) · (1 − S_w)).

---

## Auxiliary Quantities — unchanged from paper

These three quantities are computed identically to the paper (Eqs. 10–12):

**Mean-field infection rate r** (paper Eq. 11):

$$r = \frac{\displaystyle\sum_{n,i} \beta(n,i)\,(n-i)\,f_{n,i}\,p_n}{\displaystyle\sum_{n,i}(n-i)\,f_{n,i}\,p_n}$$

**External infection pressure ρ** (paper Eq. 10):

$$\rho = r \cdot \frac{\displaystyle\sum_k k(k-1)\,s_k\,g_k}{\displaystyle\sum_k k\,s_k\,g_k}$$

**Global prevalence I** (paper Eq. 12):

$$I = \sum_k (1 - s_k)\, g_k$$

---

## Comparison: Paper vs. Kernel

| Quantity | Paper (scalar ω) | Kernel w(n,i) |
|----------|-----------------|---------------|
| Switching rate | ω — same for all nodes and groups | w(n,i) = α·φ(1−φ) — depends on group state |
| Susceptible replacement probability | 1 − I — global fraction | S_w — flux-weighted fraction |
| Infected replacement probability | I — global fraction | 1 − S_w |
| i+1→i swap rate | ω(1−I) | w(n, i+1)·S_w |
| i−1→i swap rate | ωI | w(n, i−1)·(1−S_w) |
| Outflow, infected | μ + ω(1−I) | μ + w(n,i)·S_w |
| Outflow, susceptible (swap) | ωI | w(n,i)·(1−S_w) |

**Key structural change:** In the paper, every node switches at the same rate ω and the composition of incoming replacements tracks the global infected fraction I. In the kernel version, nodes in mixed groups switch faster, so the pool of available replacements is enriched with nodes from mixed groups — captured by S_w instead of (1−I).

---

## Full Kernel wAME System (collected)

$$w(n,i) = \alpha \cdot \frac{i}{n}\left(1 - \frac{i}{n}\right)$$

$$S_w = \frac{\sum_{n,i}(n-i)\,w(n,i)\,f_{n,i}}{\sum_{n,i} n\,w(n,i)\,f_{n,i}}$$

$$r = \frac{\sum_{n,i}\beta(n,i)(n-i)\,f_{n,i}\,p_n}{\sum_{n,i}(n-i)\,f_{n,i}\,p_n}, \qquad \beta(n,i) = \lambda i^\nu$$

$$\rho = r\cdot\frac{\sum_k k(k-1)\,s_k\,g_k}{\sum_k k\,s_k\,g_k}$$

$$I = \sum_k (1-s_k)\,g_k$$

$$\frac{ds_k}{dt} = \mu(1-s_k) - k\,r\,s_k$$

$$\frac{df_{n,i}}{dt} =
(i+1)\bigl(\mu + w(n,i+1)\,S_w\bigr)f_{n,i+1}$$
$$-\Bigl[i\bigl(\mu + w(n,i)\,S_w\bigr) + (n-i)\bigl(\beta(n,i)+\rho+w(n,i)(1-S_w)\bigr)\Bigr]f_{n,i}$$
$$+(n-i+1)\bigl(\beta(n,i-1)+\rho+w(n,i-1)(1-S_w)\bigr)f_{n,i-1}$$

---

## Resolved: S_w does not include p_n

**Confirmed with professor.** The scratch paper had ν (not μ) as a trailing term — a carryover from a different project where groups have different synergy exponents. It is not part of this model.

S_w has no p_n because w(n,i) is a group trait: all nodes in a group switch at the same rate, so groups are fully identified by (n, i). The sum over f_{n,i} is the correct average — no structural reweighting by p_n is needed. This differs from r, where p_n is required because f_{n,i} is conditional on n and r is a per-susceptible infection rate that must account for how many groups of each size exist.
