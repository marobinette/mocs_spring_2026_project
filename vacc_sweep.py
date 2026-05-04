#!/usr/bin/env python3
"""
VACC parameter sweep — diversity-tension kernel or baseline.

Sweeps (lambda x nu) for a single alpha value, or runs a single baseline sweep
with constant omega. Designed to be submitted as a SLURM array job, one task
per alpha value.

Kernel: w(n,i) = 4α·φ·(1−φ)  peaks at φ=0.5 (flees mixed groups)
Baseline: w(n,i) = OMEGA_BASELINE (constant)

Usage:
    python vacc_sweep.py --alpha 3.0
    python vacc_sweep.py --alpha-index 6           # index into ALPHA_VALUES
    python vacc_sweep.py --baseline
    python vacc_sweep.py --alpha 3.0 --network Synthetic_poisson_k5
    python vacc_sweep.py --baseline --network Synthetic_poisson_k5
    python vacc_sweep.py --alpha 3.0 --workers 16
"""

import argparse
import json
import os
import time
from multiprocessing import Pool, cpu_count

import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import binom

# ---------------------------------------------------------------------------
# Grid
# ---------------------------------------------------------------------------
N_LAM = 50
N_NU  = 40

LAM_GRID = np.logspace(-3, 0, N_LAM)
NU_GRID  = np.linspace(1.0, 10.0, N_NU)

I0_LOW  = 1e-3
I0_HIGH = 0.99

MU           = 1.0
TRAJ_POINTS  = 50
T_MAX        = 1000.0

DATA_PATH = "Data/group_statistics.txt"
OUT_DIR   = "Files/vacc"

# 13 log-spaced alpha values from 0.1 to 100 (index = $SLURM_ARRAY_TASK_ID)
ALPHA_VALUES = np.logspace(-1, 2, 13).tolist()

OMEGA_BASELINE = 5.0

# ---------------------------------------------------------------------------
# Network loading
# ---------------------------------------------------------------------------
def load_group_statistics(network, path=DATA_PATH):
    with open(path) as f:
        results = json.load(f)
    if network not in results:
        raise KeyError(f"Network '{network}' not in {path}. Available: {list(results.keys())}")
    data = results[network]
    nmax = int(max(data["group_size_n"]))
    mmax = int(max(data["membership_k"]))
    pn = np.zeros(nmax + 1)
    gm = np.zeros(mmax + 1)
    for n, p in zip(data["group_size_n"], data["group_size_p"]):
        pn[int(n)] = float(p)
    for k, g in zip(data["membership_k"], data["membership_g"]):
        gm[int(k)] = float(g)
    state_meta = _get_state_meta(mmax, nmax, gm, pn)
    return gm, pn, mmax, nmax, state_meta


def _get_state_meta(mmax, nmax, gm, pn):
    m = np.arange(mmax + 1)
    imat = np.zeros((nmax + 1, nmax + 1))
    nmat = np.zeros((nmax + 1, nmax + 1))
    for n in range(2, nmax + 1):
        imat[n, :n + 1] = np.arange(n + 1)
        nmat[n, :n + 1] = n
    pnmat = np.outer(pn, np.ones(nmax + 1))
    return (mmax, nmax, m, np.array(gm), np.array(pn), imat, nmat, pnmat)


# ---------------------------------------------------------------------------
# Kernel
# ---------------------------------------------------------------------------
def w_diversity_tension(n, i, alpha):
    phi = i / n
    return alpha * 4 * phi * (1 - phi)


# ---------------------------------------------------------------------------
# Physics
# ---------------------------------------------------------------------------
def _infection_matrix(lam, nu, nmax):
    mat = np.zeros((nmax + 1, nmax + 1))
    for n in range(2, nmax + 1):
        for i in range(n):
            mat[n, i] = lam * float(i) ** nu
    return mat


def _switching_matrix(nmax, alpha):
    mat = np.zeros((nmax + 1, nmax + 1))
    for n in range(2, nmax + 1):
        for i in range(n + 1):
            mat[n, i] = w_diversity_tension(n, i, alpha)
    return mat


def _baseline_switching_matrix(nmax):
    mat = np.zeros((nmax + 1, nmax + 1))
    for n in range(2, nmax + 1):
        mat[n, :n + 1] = OMEGA_BASELINE
    return mat


def _initialize(state_meta, I0):
    mmax, nmax = state_meta[0], state_meta[1]
    sm = np.full(mmax + 1, 1.0 - I0)
    fni = np.zeros((nmax + 1, nmax + 1))
    for n in range(2, nmax + 1):
        fni[n, :n + 1] = binom.pmf(np.arange(n + 1), n, I0)
    return sm, fni


def _infected_fraction(sm, gm):
    return float(np.sum((1.0 - sm) * gm))


def _vector_field(v, _t, inf_mat, w_mat, state_meta, use_baseline=False):
    mmax, nmax = state_meta[0], state_meta[1]
    m, gm = state_meta[2], state_meta[3]
    imat, nmat, pnmat = state_meta[5], state_meta[6], state_meta[7]

    sm  = v[:mmax + 1]
    fni = v[mmax + 1:].reshape((nmax + 1, nmax + 1))
    fni_field = np.zeros_like(fni)

    denom_r = np.sum((nmat[2:, :] - imat[2:, :]) * fni[2:, :] * pnmat[2:, :])
    if denom_r < 1e-14:
        r = 0.0
    else:
        r = np.sum(
            inf_mat[2:, :] * (nmat[2:, :] - imat[2:, :]) * fni[2:, :] * pnmat[2:, :]
        ) / denom_r

    denom_rho = np.sum(m * sm * gm)
    if denom_rho < 1e-14:
        rho = 0.0
    else:
        rho = r * np.sum(m * (m - 1) * sm * gm) / denom_rho

    S_w_denom = np.sum(nmat[2:, :] * w_mat[2:, :] * fni[2:, :])
    if S_w_denom > 1e-14 and not use_baseline:
        S_w = (
            np.sum((nmat[2:, :] - imat[2:, :]) * w_mat[2:, :] * fni[2:, :])
            / S_w_denom
        )
    else:
        I = _infected_fraction(sm, gm)
        S_w = 1.0 - I

    sm_field = MU * (1.0 - sm) - sm * m * r

    fni_field[2:, :nmax] += (
        imat[2:, 1:] * (MU + w_mat[2:, 1:] * S_w) * fni[2:, 1:]
    )
    fni_field[2:, :] += (
        -imat[2:, :] * (MU + w_mat[2:, :] * S_w)
        - (nmat[2:, :] - imat[2:, :]) * (inf_mat[2:, :] + rho + w_mat[2:, :] * (1.0 - S_w))
    ) * fni[2:, :]
    fni_field[2:, 1:nmax + 1] += (
        (nmat[2:, :nmax] - imat[2:, :nmax])
        * (inf_mat[2:, :nmax] + rho + w_mat[2:, :nmax] * (1.0 - S_w))
        * fni[2:, :nmax]
    )

    return np.concatenate((sm_field, fni_field.reshape((nmax + 1) ** 2)))


def _integrate(lam, nu, alpha, state_meta, I0, use_baseline=False):
    mmax, nmax = state_meta[0], state_meta[1]
    gm = state_meta[3]

    inf_mat = _infection_matrix(lam, nu, nmax)
    w_mat   = _baseline_switching_matrix(nmax) if use_baseline else _switching_matrix(nmax, alpha)
    sm, fni = _initialize(state_meta, I0)
    v0 = np.concatenate((sm, fni.reshape((nmax + 1) ** 2)))

    sol = solve_ivp(
        lambda t, v: _vector_field(v, t, inf_mat, w_mat, state_meta, use_baseline=use_baseline),
        t_span=(0.0, T_MAX),
        y0=v0,
        method="LSODA",
        t_eval=np.linspace(0.0, T_MAX, TRAJ_POINTS),
    )

    return _infected_fraction(sol.y[:mmax + 1, -1], gm)


# ---------------------------------------------------------------------------
# Parallel sweep — one worker per nu slice
# ---------------------------------------------------------------------------
def _sweep_nu_slice(args):
    j, nu, alpha, state_meta, use_baseline = args
    I_low_row  = np.zeros(N_LAM)
    I_high_row = np.zeros(N_LAM)
    for i, lam in enumerate(LAM_GRID):
        I_low_row[i]  = _integrate(lam, nu, alpha, state_meta, I0_LOW,  use_baseline=use_baseline)
        I_high_row[i] = _integrate(lam, nu, alpha, state_meta, I0_HIGH, use_baseline=use_baseline)
        print(
            f"  nu={nu:.2f} [{j+1}/{N_NU}]  lam={lam:.2e} [{i+1}/{N_LAM}]"
            f"  I*(low)={I_low_row[i]:.4f}  I*(high)={I_high_row[i]:.4f}",
            flush=True,
        )
    return j, I_low_row, I_high_row


def run_sweep(alpha, state_meta, n_workers, label, use_baseline=False):
    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"  Grid: {N_LAM} lambda x {N_NU} nu")
    print(f"  T_MAX={T_MAX}  TRAJ_POINTS={TRAJ_POINTS}  workers={n_workers}")
    print(f"{'='*60}\n")
    t0 = time.time()

    tasks = [(j, nu, alpha, state_meta, use_baseline) for j, nu in enumerate(NU_GRID)]

    I_low  = np.zeros((N_LAM, N_NU))
    I_high = np.zeros((N_LAM, N_NU))

    with Pool(n_workers) as pool:
        for j, row_low, row_high in pool.imap_unordered(_sweep_nu_slice, tasks):
            I_low[:, j]  = row_low
            I_high[:, j] = row_high

    delta = I_high - I_low
    print(f"\n  Done in {time.time() - t0:.1f}s")
    return I_low, I_high, delta


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument("--alpha", type=float,
                            help="Diversity-tension amplitude α")
    mode_group.add_argument("--alpha-index", type=int,
                            help=f"Index into ALPHA_VALUES: {ALPHA_VALUES}")
    mode_group.add_argument("--baseline", action="store_true",
                            help=f"Run baseline sweep with constant ω={OMEGA_BASELINE}")
    parser.add_argument("--workers", type=int, default=None,
                        help="Parallel workers (default: all available CPUs)")
    parser.add_argument("--network", type=str, default="Thiers13",
                        choices=["Thiers13", "Synthetic_poisson_k5", "Synthetic_poisson_k3"],
                        help="Network to sweep (default: Thiers13)")
    args = parser.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)

    network   = args.network
    n_workers = args.workers or cpu_count()
    date_str  = time.strftime('%Y-%m-%d')

    gm, pn, mmax, nmax, state_meta = load_group_statistics(network)

    if args.baseline:
        outfile = os.path.join(OUT_DIR, f"{date_str}_{network}_baseline.npz")
        I_low, I_high, delta = run_sweep(
            None, state_meta, n_workers,
            f"Baseline (ω={OMEGA_BASELINE})  network={network}",
            use_baseline=True,
        )
        np.savez_compressed(
            outfile,
            lam_grid=LAM_GRID, nu_grid=NU_GRID,
            I_low=I_low, I_high=I_high, delta=delta,
            mu=MU,
            I0_low=I0_LOW, I0_high=I0_HIGH,
            omega_baseline=OMEGA_BASELINE,
            kernel="baseline",
        )
    else:
        alpha = args.alpha if args.alpha is not None else ALPHA_VALUES[args.alpha_index]
        outfile = os.path.join(OUT_DIR, f"{date_str}_{network}_kernel_alpha_{alpha:g}.npz")
        I_low, I_high, delta = run_sweep(
            alpha, state_meta, n_workers,
            f"Diversity-tension kernel  alpha={alpha}  network={network}",
        )
        np.savez_compressed(
            outfile,
            lam_grid=LAM_GRID, nu_grid=NU_GRID,
            I_low=I_low, I_high=I_high, delta=delta,
            alpha=alpha, mu=MU,
            I0_low=I0_LOW, I0_high=I0_HIGH,
            kernel="diversity_tension",
        )

    print(f"\nSaved → {outfile}")


if __name__ == "__main__":
    main()
