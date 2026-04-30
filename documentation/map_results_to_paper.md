# Mapping VACC Sweep Results onto Figure 3

## Fig 3a and the δI heatmap are the same information, viewed from perpendicular angles

Fig 3a plots the tricritical line in **(ω, λ) space**, with ν color-coded along the line. It answers: "for a given (ω, λ) pair, which ν values produce a discontinuous transition?" The dashed line at ω=5 cuts the tricritical curve at three points — meaning at ω=5, the bistable window in λ opens and closes non-trivially as ν varies.

The **baseline δI heatmap** is a fixed-ω=5 slice of that same space, plotted directly in **(λ, ν) space**. The boundary of the δI>0 region *is* the tricritical line evaluated at ω=5, projected into (λ, ν). Where the dashed ω=5 line in Fig 3a intersects the tricritical curve three times, the heatmap shows a bistable island — a finite λ window at each ν that opens, widens, and closes again. The concordance between the high-δI region at low λ / high ν and the paper's "Discontinuous" upper-left region in Fig 3a is the direct validation that the sweep is recovering the correct phase structure.

---

## Fig 3b and our approach diverge in one important way

Fig 3b shows full I\*(λ) bifurcation curves at each ν — continuous segments (stable branches), dashed segments (unstable), and jump discontinuities. The quantity δI = I\*(high) − I\*(low) collapses each bifurcation diagram to a single number: the gap between the two attractors reachable by the two initial conditions.

This works cleanly for **Fig 3d** (classic absorbing-active bistability) — I0_low lands on the absorbing state and I0_high on the active state. But Fig 3b also shows two richer regimes:

- **Fig 3e (hybrid continuous with active bistability):** two *active* branches coexist over a finite λ range. I0_low may climb to the lower active branch while I0_high falls to the upper. δI is nonzero, but measures a localized-vs-delocalized gap rather than an absorbing-active gap. The two cases are indistinguishable in the heatmap.
- **Fig 3f (three coexisting states):** absorbing + localized active + delocalized active. The middle branch is reachable only from a narrow IC window. With two ICs, it is almost certainly missed entirely. Any δI signal comes from the outer two states only.

**Practical implication:** bistable cell counts are a lower bound on the true richness of the phase structure. At high ν (near the ν=12 grid ceiling), where the three-state regime is most likely, the heatmap may be giving a conservative picture.

---

## The kernel results map onto moving along the ω axis of Fig 3a

For the constant-ω baseline, the system sits at a fixed horizontal slice (ω=5 dashed line) in Fig 3a. Increasing α with the diversity-tension kernel is not simply moving that line upward — the effective ω is state-dependent. But the mechanistic effect is analogous to moving the slice downward: at the endpoints (φ=0 or 1), w→0, which corresponds to ω→0 in Fig 3a's lower-left "Continuous" region. The system cannot remain in the "Discontinuous" upper region because the epidemic dynamics themselves push group states toward the endpoints where the kernel vanishes.

Above the critical α (~5, where peak kernel rate w=α/4 exceeds μ=1), this self-quenching happens faster than the disease can exploit the mixing. The system is effectively locked into the Continuous regime regardless of λ and ν. The α series of heatmaps shows the bistable region shrinking and collapsing as the kernel increasingly confines the effective dynamics to the lower-left corner of Fig 3a.

---

## Structural caveat: the synthetic network likely has a simpler tricritical line

The paper explicitly notes that the non-monotonic, folded tricritical line (which causes the ω=5 dashed line to intersect the curve three times) is a consequence of heterogeneous group-size and membership distributions. For homogeneous structures, the tricritical line is monotonic (paper Fig 2, Eq. 3).

Synthetic_delta_k5 is a regular network. Depending on whether "delta_k5" refers to a delta distribution in group sizes or in membership, the tricritical line in this network may be simpler and more monotonic than Fig 3a's folded curve. This means:

- The three-state regime (Fig 3f) is less likely to appear.
- The bistable region is probably a single connected island rather than the complex multi-lobe structure of the high school dataset.

This is a feature rather than a bug: the synthetic network is a cleaner test of the kernel's suppression effect precisely because the background phase structure is simpler and easier to interpret.
