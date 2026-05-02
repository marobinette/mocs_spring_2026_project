# Summary of Discussion with Laurent Hébert-Dufresne
## Diversity-Tension Kernel — Project Context and Code Implications

---

## 1. What Figure 3 Is and Why It Matters

Laurent walked through Figure 3 of the paper in detail. It is the central result — a collection of phase transition portraits for the high school network showing how the stationary prevalence I* depends on the transmission rate λ, for seven different values of the synergy exponent ν.

The vertical axis is I*, the equilibrium fraction of adopters. The horizontal axis is λ, the baseline transmission rate. The key insight is that the *shape* of this curve changes dramatically depending on ν:

- **ν ≈ 1 (linear contagion):** Classic continuous SIS-like transition. I* grows smoothly from zero as λ crosses the critical point. Delocalized — activity is spread broadly across groups.
- **ν moderate (2–4):** Discontinuous transition with bistability. I* jumps from zero to a large value at λc. This is the classic absorbing-active bistability the paper is most known for.
- **ν larger still (5–7):** Something surprising happens — a *hybrid continuous* transition appears before the discontinuous jump. Activity first emerges continuously in the largest groups (localized), then jumps discontinuously to a delocalized state. Two stable active branches coexist.
- **ν ≈ 10:** The richest regime — three coexisting stable states. The absorbing state, a localized active state sustained by the largest groups, and a delocalized active state. Getting to the localized branch requires fine-tuning the initial condition to a narrow window.

Laurent emphasized that the last two regimes had never been observed before in contagion models. He called them "new phases of contagious matter."

---

## 2. How Santi Gets All Three Lines

This connects directly to the parameter sweep code. Laurent explained that to recover all stable branches in the richest regime (three stable states), you need to sample many initial conditions densely because the middle branch is only reachable from a narrow range of I0.

However, for most regimes two initial conditions are sufficient:

- **I0 close to 0** — the system climbs to the first stable active branch it encounters from below.
- **I0 close to 1** — the system relaxes to the first stable active branch it encounters from above.

These two traces together recover the upper and lower stable branches. The middle unstable branch (dashed lines in the figure) is not directly observable — it is the boundary between basins of attraction. The fact that I*(high) ≠ I*(low) at a given λ is itself the signal that bistability exists.

**This is exactly why the sweep code uses two initial conditions per cell** (`I0_low = 1e-3`, `I0_high = 0.99`) and records both final values. The bistability indicator `delta_I = I*(high) - I*(low)` is large wherever two branches coexist and collapses to zero where there is only one stable state.

---

## 3. The Parameter Sweep Strategy

Laurent's recommendation for the sweep was clear and is directly reflected in the code:

> "Start with two initial conditions all the time — one close to zero, one close to one. Sweep lambda for different extreme values of nu. That should give you an idea of what kind of system your network is producing."

The sweep grid covers:
- **λ:** logspace from 1e-4 to 1e-1 (15 points locally, to be expanded on VACC)
- **ν:** linear from 1.0 to 12.0 (15 points locally)

Laurent's bet on where to look is at the two extremes of ν:
- ν ≈ 1: either classic continuous or S-shaped if localization is present
- ν ≈ 10: the rich multistable regime, if the network's heterogeneity supports it

The intermediate values of ν will fill in the picture between those extremes.

---

## 4. The Tricritical Line — Why We Are Not Computing It

The paper's analytical shortcut for knowing where to look in parameter space is the tricritical line — the boundary in (ω, λ, ν) space separating continuous transitions from discontinuous ones. This is Equation 16 in the paper, derived by setting the second derivative of the self-consistency function M to zero.

Laurent discussed at length why this is not straightforward with the diversity-tension kernel. In the paper, ω is a scalar constant that gets pulled outside of sums during the derivative calculations. With the kernel, ω is w(n, i) — a function of the group state — so it cannot be pulled outside. Every sum that previously simplified now has an extra (n, i) dependency woven through it.

The conclusion was:

> "This is like a big physics project. Maybe we stick with the parameter sweep. That could be like a summer thing — or for the master's student in Spain with a physics background."

**Practical implication:** The parameter sweep is not a fallback or a shortcut. It is the right approach given that the analytical tricritical condition does not transfer cleanly to the kernel model. The heatmap of δI across (λ, ν) space is the computational equivalent of the tricritical line — it maps where the interesting regimes live without requiring the analytical derivation.

---

## 5. The Network Choice — CNS vs Thiers13

Laurent's discussion of Figure 3 was based on the high school network (Thiers13), which is the paper's main empirical example. However, there is an important structural difference relevant to the kernel.

Thiers13 has gk ≈ δk,1 — almost everyone belongs to exactly one group at a time. The paper's clean analytical results are derived in this homogeneous membership regime. The diversity-tension kernel, by contrast, only becomes dynamically interesting when people belong to multiple groups simultaneously — the tension in one group can push someone toward other groups they already belong to, creating competition between memberships.

**The CNS network** (loaded in Cell 8) has genuine multi-membership structure, making it the more natural testbed for the kernel. This also means direct comparison to the paper's analytical curves is not the goal. The baseline comparison is the constant-ω model run on the same CNS network — that is what isolates the kernel's contribution from structural differences.

---

## 6. The Framing — Reversing the Causal Arrow

The deepest contribution Laurent helped articulate is a reversal of the paper's causal logic:

**Paper:** Group dynamics (plasticity rate ω) → shapes → contagion dynamics (phase portrait)

**Your model:** Contagion dynamics (local prevalence i/n) → shapes → group dynamics (switching rate w(n,i)) → feeds back into → contagion dynamics

In the paper, ω is a control parameter set from outside. In your model, the effective plasticity the system experiences is an emergent property — it depends on where the system is in the phase portrait at every moment. The system is coevolutionary rather than driven.

This has concrete dynamical consequences:
- Early spread (i/n small everywhere): w ≈ 0, system is quasi-quenched, invasion must happen with almost no mixing benefit
- Endemic state (groups pushed to endpoints by kernel): w → 0 again, system self-quenches
- The transient (broad i/n distribution, high tension): maximum fluidity, reinforcement hardest to accumulate

---

## 7. The Diversity-Tension Framing

Laurent articulated the social mechanism clearly near the end of the discussion:

> "If the group is 50-50, people have to do too much code switching. Whereas if you're the only one in the minority, you just use the code of the other people. Or you find your people and stick with that."

This maps directly onto the kernel:
- Homogeneous groups (i/n = 0 or 1): w = 0, no tension, people stay
- Maximally mixed groups (i/n = 0.5): w = α/4, maximum tension, people leave

The name settled on is **diversity tension**. The mechanism is allegiance and homophily — people prefer to be in groups where their state is the norm, and flee groups where they are caught between two roughly equal factions.

---

## 8. Immediate Next Steps

In order of priority:

1. **Let the current sweep finish** — kernel run and constant-ω baseline both complete, results saved to `Files/sweep_local.npz`
2. **Read the heatmap** — look for where δI is nonzero in the kernel panel, compare to baseline panel, examine the difference plot
3. **Validate the extremes** — does ν ≈ 1 give continuous transition? Does ν ≈ 10 give any bistability signal?
4. **Calibrate the λ range** — current results suggest the transition boundary may sit below the grid's lower end; may need to extend λ downward on the VACC sweep
5. **Design the VACC sweep** — larger grid, multiple α values, denser I0 sampling in cells that show bistability signal
6. **Draft the abstract** — once at least one concrete finding about bistability (preserved, suppressed, or shifted) is in hand

---

## 9. Open Questions

- Does bistability survive the diversity-tension kernel, or does self-quenching at the endpoints stabilize the system into a single active state?
- Is the kernel's effect on the invasion threshold a suppression (quasi-quenched early dynamics) or enhancement (homogeneous groups freeze once active)?
- Does the kernel produce any regime that looks qualitatively different from what the paper found — or does it primarily shift where in (λ, ν) space the known regimes appear?
- Is the parabolic kernel shape strong enough, or does a threshold or asymmetric kernel produce more distinct behavior?
- What is the right α to use — and how does the phase portrait change as α varies?
