# Network Features and Their Role in the wAME Model

## The two state variables

The ODE tracks two things simultaneously:

- `sm[m]`: the susceptible fraction among nodes that are simultaneously in `m` groups
- `fni[n,i]`: the fraction of size-`n` groups that currently have `i` infected members

`p(n)` governs the `fni` half. `g(k)` governs the `sm` half. They couple through shared quantities computed from both.

---

## How `p(n)` enters the model

**The force of infection `r` (vacc_sweep.py:140–142):**

```
r = Σ_{n,i} [λiᵛ · (n−i) · fni[n,i] · p(n)] / Σ_{n,i} [(n−i) · fni[n,i] · p(n)]
```

This is the per-capita infection rate a susceptible node faces inside one group. `p(n)` is the weight: size-`n` groups contribute to both the infection numerator (how many infected contacts are generating transmission) and the denominator (how many susceptible slots exist). The group size distribution determines which group sizes dominate this average.

**The switching fraction `S_w` (vacc_sweep.py:154–161):**

```
S_w = Σ_{n,i} [(n−i) · w(n,i) · fni[n,i] · p(n)] / Σ_{n,i} [n · w(n,i) · fni[n,i] · p(n)]
```

This is the susceptible fraction among nodes that are actively switching — weighted by how much each group is switching. `p(n)` determines which group sizes contribute most.

**The diversity-tension kernel `w(n,i) = α·(i/n)·(1−i/n)` (vacc_sweep.py:85–87):**

This kernel only takes a meaningful range of values if groups have enough members. For a group of size `n=2`, `i` can only be 0, 1, or 2, giving switching rates of 0, `0.25α`, or 0. For larger `n`, the kernel traces a smooth parabola with many intermediate states. So `p(n)` — specifically `nmax` and the tail of the distribution — shapes how richly the switching kernel can express itself.

**The `fni` dynamics themselves (vacc_sweep.py:166–177):**

Groups transition between states `[n,i]` → `[n,i±1]` through three mechanisms: within-group infection at rate `λiᵛ`, background infection at rate `ρ` (see below), and switching at rate `w(n,i)`. The full shape of `p(n)` determines how these processes are distributed across group sizes.

---

## How `g(k)` enters the model

**The total infected fraction `I` (vacc_sweep.py:122–123):**

```
I = Σ_m (1 − sm[m]) · g(m)
```

`g(k)` is the weight: the overall infected fraction is a weighted average of the infected fraction within each membership class. If `g(k)` is concentrated at `k=1`, `I` essentially tracks just `sm[1]`. If `g(k)` has mass at large `k`, high-membership nodes contribute proportionally.

**The per-node infection rate in `sm` (vacc_sweep.py:163):**

```
dsm[m]/dt = μ·(1 − sm[m]) − sm[m]·m·r
```

A node in `m` simultaneous groups faces infection pressure `m·r` — each group contributes `r` independently. This is the core reason multi-group membership matters: every doubling of `m` doubles the infection rate that class faces. If `g(k)` has mass at `k=10`, those nodes are infected at 10 times the rate of singly-connected nodes.

**The background infection term `ρ` (vacc_sweep.py:145–149):**

```
ρ = r · Σ_m [m·(m−1)·sm[m]·g(m)] / Σ_m [m·sm[m]·g(m)]
```

This is the most structurally important piece. `ρ` represents infection spreading through shared group memberships — if you're in `m` groups, you have `m−1` additional memberships beyond any one group where you could encounter infected nodes. The key ratio is approximately `⟨m(m−1)⟩ / ⟨m⟩` over susceptible nodes.

Crucially: if `m=1` for all nodes, `m·(m−1) = 0` and `ρ = 0` exactly. The shared-membership infection pathway vanishes entirely. Multi-group membership is what turns `ρ` on, and the second moment `⟨m(m−1)⟩` amplifies it faster than linearly.

---

## The problem with CNS

CNS has `kmax=233`. Even if `g(233)` is tiny, the term `m·(m−1) = 233·232 = 54,056` means those nodes contribute ~54,000 times more to the `ρ` numerator than a singly-connected node. The fat tail of CNS's membership distribution makes `ρ` enormous. Meanwhile in the `sm` equation, those nodes face infection rate `233r` — they are practically always infected.

The deeper problem: `g(k)` for CNS was computed cumulatively across the full observation window, not as simultaneous membership at any given moment. A node's `k=233` means it appeared in 233 distinct groups over the entire dataset — sequentially, not at the same time. None of this reflects genuine simultaneous network structure. The large `ρ` and strong high-`m` dynamics are artifacts of aggregation, not physics.

---

## The problem with Thiers13

Thiers13 has `⟨k⟩ = 1.043`. Most nodes have `m=1`. That means:

- `ρ ≈ 0` — the shared-membership pathway is essentially off
- `sm[m]` for `m≥2` carries almost no weight in `I`
- The `m·r` term is just `1·r` for nearly everyone

The `fni` dynamics (switching, within-group infection) are fully active, but the `sm`/`g(k)` side of the model is doing almost nothing. You're running a model with a rich multi-group membership structure using a network where nodes are almost never simultaneously in more than one group.

---

## What genuine multi-group membership means mechanically

Having nodes in more than one group simultaneously means:

1. `ρ > 0` — the model's second infection pathway is active
2. `dsm[m]/dt` with `m≥2` carries real weight — different membership classes evolve differently
3. The `⟨m(m−1)⟩` term creates a nonlinear amplification that can drive bistability distinct from what the within-group kernel alone produces

None of the SocioPatterns datasets give you this cleanly. They all have near-zero `ρ`. CNS has nonzero `ρ` but for spurious reasons. The honest option is either to construct a synthetic `g(k)` with the desired properties, or to find a different class of empirical network (affiliation networks, organizational memberships) where simultaneous multi-group membership is structurally guaranteed.

---

## Network summary

### Empirical networks

| Network    | `⟨k⟩` | `kmax` | `⟨k(k-1)⟩/⟨k⟩` | `⟨n⟩` | `nmax` | Timesteps |
|------------|--------|--------|-----------------|--------|--------|-----------|
| CNS        | 2.330  | 233    | 10.691          | 3.578  | 20     | 8,064     |
| Thiers13   | 1.043  | 5      | 0.090           | 2.312  | 6      | 18,179    |
| LyonSchool | 1.094  | 5      | —               | 2.499  | —      | 5,846     |
| SFHH       | 1.047  | 4      | —               | 2.422  | —      | 5,707     |
| InVS15     | 1.023  | 6      | —               | 2.185  | —      | 49,677    |
| Malawi     | 1.001  | 3      | —               | 2.055  | —      | 57,789    |
| LH10       | 1.032  | 4      | —               | 2.411  | —      | 12,278    |
| InVS13     | 1.022  | 4      | —               | 2.010  | —      | 49,322    |

The `⟨k(k-1)⟩/⟨k⟩` column is the ρ/r prefactor — how strongly the shared-membership infection pathway is activated relative to within-group infection. It is exactly zero when all nodes have k=1.

### Synthetic networks

All synthetic networks use Thiers13's empirical `p(n)` for group sizes. Only `g(k)` varies. Generated by `make_synthetic_networks.py`.

| Network                 | `⟨k⟩` | `kmax` | `⟨k(k-1)⟩/⟨k⟩` | `g(k)` type              |
|-------------------------|--------|--------|-----------------|--------------------------|
| Synthetic_delta_k1      | 1.000  | 1      | 0.000           | δ_{k,1} (baseline)       |
| Synthetic_delta_k2      | 2.000  | 2      | 1.000           | δ_{k,2}                  |
| Synthetic_delta_k3      | 3.000  | 3      | 2.000           | δ_{k,3}                  |
| Synthetic_delta_k5      | 5.000  | 5      | 4.000           | δ_{k,5}                  |
| Synthetic_poisson_k2    | 2.000  | 15     | 1.594           | zero-truncated Poisson(2) |
| Synthetic_poisson_k3    | 3.000  | 19     | 2.821           | zero-truncated Poisson(3) |
| Synthetic_poisson_k5    | 5.000  | 25     | 4.965           | zero-truncated Poisson(5) |

**Delta vs Poisson at the same `⟨k⟩`:** The Poisson always has a higher ρ/r prefactor because its variance amplifies `⟨k(k-1)⟩`. Comparing delta and Poisson at the same mean separates the effect of the *mean* from the *spread* of multi-group membership.
