# Project Outline — *The Tension of Diversity in Group Interactions in Complex Contagion*

This is a planning document, not the paper. It captures (a) where the current draft is weak, (b) what edits to make to the introduction and methods, (c) the recommended results panel, and (d) how to handle the two kernels.

---

## 1. Overall assessment of the current draft

Strengths:
- The motivation is genuinely novel — making ω depend on group composition is a clean, well-posed extension of Lamata-Otín et al.
- The derivation of $S_w$ (susceptible-share of the switching flux) is the right conceptual fix once ω becomes state-dependent. This is the part of the methods that earns the paper its keep — make sure it stays prominent.
- You have already produced enough figures (ridgelines, α-grid heatmaps, steady-state curves on three networks for the diversity-tension kernel, plus the k=2 tension-shifted kernel) that the results section can be filled in immediately.

Weak spots that need fixing before this is publication-grade:
1. **Abstract is empty.** Write last, but it should land in 4 sentences: (i) higher-order temporal contagion produces bistability that depends on a switching rate ω; (ii) prior work treats ω as a constant; (iii) we make ω depend on local infection composition via two kernels and re-derive the wAME closure with a flux-corrected swap probability $S_w$; (iv) state-dependent switching can suppress, preserve, or qualitatively reshape the bistable regime depending on whether the kernel quenches at homogeneous compositions.
2. **Results and Discussion are placeholders.** Cover photo of a forest / "Condition A/B" template table need to come out.
3. **Several model objects are used but not defined in main text.** $\beta(n,i)$, $s_k$, $g_k$, $p_n$, the membership-state ODE that closes $s_k$, and the synergy exponent $\nu$ all appear in figures and equations but are never introduced. A reader who has not read Lamata-Otín cannot follow §2.
4. **Equation (1) has bracket/formatting glitches** — the inline curly braces don't render cleanly and the second line is missing a closing bracket. Re-typeset carefully.
5. **The introduction ends at the kernel formula but never states a research question or what the paper will show.** Add an explicit "we find that…" paragraph.
6. **"TO-DO: Reference network sources"** — fill in (SocioPatterns for Thiers13/LyonSchool; cite the Poisson-network construction and the CNS dataset if used).
7. **Closing parenthesis** on the kernel expression at end of intro is missing: `α·(i/n)(1−i/n)`.
8. **Only one kernel is defined in §2.** If the tension-shifted kernel is going to appear in results, it must be defined here.

---

## 2. Introduction — recommended rewrite plan

Keep the current three-paragraph spine (HO contagion → temporality → state-dependent extension), but tighten and add a contribution paragraph. Specifically:

- **¶1 (HO contagion):** Trim. The classroom/workplace examples are fine, but the sentence "Repeated exposure to ideas, beliefs, and even products can create pressure to engage and adopt that is completely absent in one on one interactions" is doing work that one phrase ("collective reinforcement") can do. Cite Iacopini et al. and one or two follow-ups.
- **¶2 (Lamata-Otín's contribution):** Good content but compress. The key facts the reader needs are: ω governs composition turnover; in the static limit (ω→0) groups freeze; in the annealed limit (ω→∞) groups are continually re-randomized; intermediate ω produces multistability and discontinuous transitions in regimes where neither limit does. Make those three timescale regimes explicit — they motivate why making ω state-dependent matters.
- **¶3 (gap):** Sharpen the homophily framing. Instead of "a node in a fully susceptible or fully infected group may have little reason to leave," say it more mechanistically: cognitive dissonance, peer-pressure, opinion-dynamics literature on ingroup preference (cite Centola, McPherson on homophily). The point is that the assumption "ω is independent of state" is a *modeling choice*, not a physical fact, and there are two opposite plausible alternatives.
- **¶4 (NEW — contribution):** Add a paragraph that does three things:
  1. Define the two kernels you study side-by-side: the **diversity-tension kernel** ω(n,i) = α·φ·(1−φ), which vanishes at the homogeneous endpoints; and the **tension-shifted kernel** ω(n,i) = α·φ·(1−φ) + ω₀, which preserves a baseline turnover rate.
  2. State the central finding in one sentence: *the diversity-tension kernel collapses bistability above a critical α by freezing the absorbing compositions, while the tension-shifted kernel preserves bistability and continuously deforms the phase boundary.*
  3. Briefly note that other kernel shapes (e.g. avoidance, where switching peaks at the homogeneous extremes) are natural future extensions of the same framework.
- **Notation cleanup at end of intro:** define $S_w$ in one sentence and forward-reference the methods ("derived in §2"), don't try to motivate it inline.

---

## 3. Methods — recommended rewrite plan

The methods do the heavy lifting. Restructure as four short subsections:

1. **Background: the wAME closure.** Half a paragraph stating what $f_{n,i}$, $s_k$, $g_k$, $p_n$, $\beta(n,i) = \lambda i^\nu$, and the order parameter $I = \sum_k (1-s_k) g_k$ are. Currently $\beta$ and $\nu$ appear in the figures but are never defined.
2. **Kernel definitions (NEW subsection).** Explicitly write down the two kernels with a small two-panel figure (`figures/diversity_tension_kernel.png` already exists — use that as Fig. 1, but add a second panel for the tension-shifted form). State the limits: as α→0, both kernels reduce to either ω=0 (DT) or ω=ω₀ (TS), recovering known regimes. As α→∞ with TS, the homogeneous endpoints still receive a finite ω₀ contribution.
3. **State-dependent ME and the swap probability $S_w$.** Keep the current presentation of Eq. (1). Two fixes:
   - Repair the bracket typesetting (use `\left[`, `\right]` or `aligned`).
   - After Eq. (3), explain *why* $S_w$ is needed in one sentence — currently the reader has to figure it out from "switching events are now enriched from groups with high local diversity." Better: "Because the switching flux is biased toward mixed groups, the global susceptible fraction $1-I$ is no longer the correct replacement probability; $S_w$ corrects for this bias."
4. **Numerical analysis.** Keep, but:
   - Move "Sweeps are run on five networks…" into its own paragraph and fill in citations.
   - State explicitly that the baseline ω=5 was chosen because it sits in the middle of the discontinuous regime in Lamata-Otín et al. Fig. 3 (so the reader knows what the comparison is *to*).
   - Add a one-line validation note: at α→0 (DT) and at α=0,ω₀=5 (TS), our integrator reproduces their published δI heatmap, which we use as a regression test.

---

## 4. Should we keep both kernels, or just the tension-shifted one?

**Recommendation: keep both, lead with the diversity-tension kernel, and frame them as a controlled pair.**

The temptation is to drop the diversity-tension kernel because it produces a "trivial" outcome (above critical α the system always hits a quenched limit and bistability dies). That undersells the result. The diversity-tension kernel is the *cleaner mechanistic statement of the homophily hypothesis* — switching is purely tension-driven, with no baseline rate — and the fact that this assumption is enough to destroy the bistability that Lamata-Otín et al. discovered is the paper's most interesting finding. The tension-shifted kernel then plays the role of the *robustness check*: when a baseline turnover ω₀ is added, the qualitative bistability picture survives, so the absorbing-endpoint mechanism is what kills bistability, not the kernel shape per se.

This is a much stronger story than either kernel alone:
- DT alone: "we made up a kernel and it killed the effect" — sounds like a flaw.
- TS alone: "we made up a kernel and the picture deforms a bit" — sounds like a perturbation paper.
- DT + TS together: "*it's the freezing of homogeneous compositions, not the dependence on local diversity, that controls bistability*" — that's a mechanistic claim, and that's what's worth publishing.

Frame the avoidance kernel (peaks at the homogeneous extremes; mentioned in the kernel-comparison figure) explicitly in the discussion as a natural next step that the same wAME-with-$S_w$ machinery handles unchanged. Don't show avoidance results in this paper unless you already have them — promising too much in the discussion is fine; including a half-finished panel is not.

---

## 5. Results section — what to show, in order

You already have most of these figures. The order should march the reader through the mechanism:

**Fig. 1 — Kernels.** Two-panel comparison of ω(n,i) vs i/n for DT and TS at several α. Use `figures/diversity_tension_kernel.png` as the basis and add the TS panel. Annotate the homogeneous endpoints.

**Fig. 2 — Validation against Lamata-Otín baseline.** δI heatmap over (λ, ν) at constant ω=5, reproducing Fig. 3 of Lamata-Otín et al. as a horizontal slice. One panel; small. This is the regression test that earns trust for everything that follows.

**Fig. 3 — Kernel α-sweep on a single network (Synthetic Poisson k=2).** A 2×N grid: top row is DT for α ∈ {0.5, 2, 5, 10, 20, 50}; bottom row is TS for the same α (with ω₀ chosen so the kernel maximum matches DT at each α — be explicit about how ω₀ is set). Use the already-generated `vacc_sweep_alpha_grid.png` figures. Caption should call out the qualitative difference: bistability collapse in DT, continuous deformation in TS.

**Fig. 4 — Mechanism: stationary group-state distributions.** Use the existing `ridgeline_stationary` figures. Show $f_{n,i}$ at the high-IC steady state for one (λ,ν) point inside the bistable region, for α ∈ {0, 5, 50}, both kernels. The reader should see mass piling up at i=0 and i=n in DT as α grows, while in TS the distribution spreads. This figure is the visual proof of the absorbing-endpoint mechanism.

**Fig. 5 — Critical α curve / bistability collapse.** From `vacc_analysis.ipynb`: for each α, plot the number of bistable cells (δI > 0.05) as a function of α. DT and TS as two curves. Mark the critical α range (between 5 and 10 from the existing data). This is the headline quantitative result.

**Fig. 6 — Network dependence.** A small multi-panel showing the bistable-cell-count vs α curve for the three networks where you have data: Synthetic Poisson k=2, k=3, and Thiers13 (and add LyonSchool / k=5 if those runs finish in time). This grounds the claim that the mechanism is robust across structurally different networks rather than an artifact of one $p_n$.

If space is tight (this is a 4-page format), Fig. 4 and Fig. 6 are the candidates to drop or move to a supplement.

### Results not yet generated that would noticeably strengthen the paper

- **Time-series traces** of $I(t)$ from both initial conditions for one bistable point and one non-bistable point, both kernels. One small inset per panel; clarifies what δI is actually measuring.
- **Critical-α scan with finer resolution** between α=5 and α=10 for the DT kernel — the existing sweep jumps from δI≈1 at α=5 to δI≈0 at α=10, which leaves the transition shape unresolved. A 5-point scan in that interval would let you say whether the bistability collapse is itself sharp or smooth.
- **Sensitivity to ω₀** in the TS kernel — pick one (λ,ν) bistable point and sweep ω₀ ∈ [0, 10] at fixed α. Shows that ω₀ → 0 recovers DT continuously, which closes the loop on the "freezing endpoints is what matters" claim.
- **Sw diagnostic.** Plot $S_w$ vs time for a representative trajectory, alongside the global $1-I$. The gap between them is the quantitative justification for the closure correction in §2 — currently the reader has to take it on faith.

The first three are achievable with the existing VACC pipeline and shouldn't take more than a day of compute. The $S_w$ diagnostic is a notebook plot — minutes.

---

## 6. Discussion — bullet skeleton

- Restate the mechanism finding (endpoint-freezing, not kernel shape, controls bistability collapse).
- Connect to homophily / opinion-dynamics literature: this is the first time, to your knowledge, that the contagion-driven feedback on group composition has been written down inside the wAME framework.
- Limitations: (a) deterministic mean-field — no finite-size noise, no stochastic escape from the absorbing compositions; (b) kernel shapes are phenomenological — no derivation from a microscopic utility function; (c) only SIS, not SIR or threshold dynamics.
- Future work: avoidance kernel; coupling kernel parameters to node attributes (heterogeneous tension); empirical estimation of ω(n,i) from a temporal contact dataset where infection state is observable (e.g. flu-tracking studies on the same SocioPatterns deployments).

---

## 7. Concrete next-step checklist

- [ ] Re-typeset Eq. (1); add closing paren on intro kernel formula.
- [ ] Define $\beta(n,i)$, $\nu$, $s_k$, $g_k$, $p_n$ in §2.1.
- [ ] Add tension-shifted kernel definition in §2.2 with explicit ω₀ scaling rule.
- [ ] Add one-paragraph validation note (baseline reproduces Lamata-Otín).
- [ ] Write contribution paragraph at end of §1.
- [x] Generate: TS sweep figures for k=3 and Thiers13 to match the DT coverage.
- [ ] Generate: fine α-scan in [5, 10] for DT kernel.
- [ ] Generate: $S_w$ vs $1-I$ time-series plot.
- [ ] Replace placeholder forest photo and Condition A/B table.
- [ ] Fill citations for Thiers13, LyonSchool, Poisson networks, Iacopini, Centola/McPherson.
- [ ] Write abstract last.
