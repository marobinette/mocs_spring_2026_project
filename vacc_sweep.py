#!/usr/bin/env python3
"""
VACC parameter sweep — diversity-tension kernel vs constant-omega baseline.

Sweeps (lambda x nu) for a single alpha value. The baseline (constant omega)
can also be run with --baseline. Designed to be submitted as a SLURM array
job, one task per alpha value.

Usage:
    python vacc_sweep.py --alpha 30.0
    python vacc_sweep.py --alpha 30.0 --workers 16
    python vacc_sweep.py --baseline
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
# Grid — scale these up from the local 15x15 test
# ---------------------------------------------------------------------------
N_LAM = 50
N_NU  = 40

LAM_GRID = np.logspace(-3, 0, N_LAM)    # 1e-3 to 1e0, matching fig_3 axes
NU_GRID  = np.linspace(1.0, 10.0, N_NU)

I0_LOW  = 1e-3
I0_HIGH = 0.99

MU           = 1.0
TRAJ_POINTS  = 500
T_MAX        = 200.0   # 200 characteristic times — safe for near-threshold convergence

OMEGA_SCALAR = 5.0     # constant-omega baseline value
NETWORK      = "Thiers13"
# NETWORK      = "Synthetic_poisson_k5"

DATA_PATH    = "Data/group_statistics.txt"
OUT_DIR      = "Files/vacc"

# Alpha values for the SLURM array job (index = $SLURM_ARRAY_TASK_ID)
ALPHA_VALUES = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 50.0]

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
# Kernels
# ---------------------------------------------------------------------------
def w_diversity_tension(n, i, alpha):
    phi = i / n
    return alpha * 4 * phi * (1 - phi)


def w_constant(n, i, omega):
    return omega


# ---------------------------------------------------------------------------
# Physics
# ---------------------------------------------------------------------------
def _infection_matrix(lam, nu, nmax):
    mat = np.zeros((nmax + 1, nmax + 1))
    for n in range(2, nmax + 1):
        for i in range(n):
            mat[n, i] = lam * float(i) ** nu
    return mat


def _switching_matrix(w_func, nmax, w_args):
    mat = np.zeros((nmax + 1, nmax + 1))
    for n in range(2, nmax + 1):
        for i in range(n + 1):
            mat[n, i] = w_func(n, i, *w_args)
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


def _vector_field(v, _t, inf_mat, w_mat, state_meta, mu, is_omega_constant=False):
    mmax, nmax = state_meta[0], state_meta[1]
    m, gm = state_meta[2], state_meta[3]
    imat, nmat, pnmat = state_meta[5], state_meta[6], state_meta[7]

    sm  = v[:mmax + 1]
    fni = v[mmax + 1:].reshape((nmax + 1, nmax + 1))
    fni_field = np.zeros_like(fni)

    # Infection pressure r — guard against all-infected state
    denom_r = np.sum((nmat[2:, :] - imat[2:, :]) * fni[2:, :] * pnmat[2:, :])
    if denom_r < 1e-14:
        r = 0.0
    else:
        r = np.sum(
            inf_mat[2:, :] * (nmat[2:, :] - imat[2:, :]) * fni[2:, :] * pnmat[2:, :]
        ) / denom_r

    # Excess susceptible membership rho — guard against all-infected state
    denom_rho = np.sum(m * sm * gm)
    if denom_rho < 1e-14:
        rho = 0.0
    else:
        rho = r * np.sum(m * (m - 1) * sm * gm) / denom_rho

    I = _infected_fraction(sm, gm)

    # Effective susceptible fraction among switchers
    S_w_denom = np.sum(nmat[2:, :] * w_mat[2:, :] * fni[2:, :])
    if S_w_denom > 1e-14 and not is_omega_constant:
        S_w = (
            np.sum((nmat[2:, :] - imat[2:, :]) * w_mat[2:, :] * fni[2:, :])
            / S_w_denom
        )
    else:
        S_w = 1.0 - I

    sm_field = mu * (1.0 - sm) - sm * m * r

    # Group state transitions
    fni_field[2:, :nmax] += (
        imat[2:, 1:] * (mu + w_mat[2:, 1:] * S_w) * fni[2:, 1:]
    )
    fni_field[2:, :] += (
        -imat[2:, :] * (mu + w_mat[2:, :] * S_w)
        - (nmat[2:, :] - imat[2:, :]) * (inf_mat[2:, :] + rho + w_mat[2:, :] * (1.0 - S_w))
    ) * fni[2:, :]
    fni_field[2:, 1:nmax + 1] += (
        (nmat[2:, :nmax] - imat[2:, :nmax])
        * (inf_mat[2:, :nmax] + rho + w_mat[2:, :nmax] * (1.0 - S_w))
        * fni[2:, :nmax]
    )

    return np.concatenate((sm_field, fni_field.reshape((nmax + 1) ** 2)))


def _integrate(lam, nu, w_func, w_args, state_meta, I0, is_omega_constant=False):
    mmax, nmax = state_meta[0], state_meta[1]
    gm = state_meta[3]

    inf_mat = _infection_matrix(lam, nu, nmax)
    w_mat   = _switching_matrix(w_func, nmax, w_args)
    sm, fni = _initialize(state_meta, I0)
    v0 = np.concatenate((sm, fni.reshape((nmax + 1) ** 2)))

    sol = solve_ivp(
        lambda t, v: _vector_field(v, t, inf_mat, w_mat, state_meta, MU, is_omega_constant),
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
    j, nu, w_func, w_args, state_meta, is_omega_constant = args
    I_low_row  = np.zeros(N_LAM)
    I_high_row = np.zeros(N_LAM)
    for i, lam in enumerate(LAM_GRID):
        I_low_row[i]  = _integrate(lam, nu, w_func, w_args, state_meta, I0_LOW,  is_omega_constant)
        I_high_row[i] = _integrate(lam, nu, w_func, w_args, state_meta, I0_HIGH, is_omega_constant)
        print(
            f"  nu={nu:.2f} [{j+1}/{N_NU}]  lam={lam:.2e} [{i+1}/{N_LAM}]"
            f"  I*(low)={I_low_row[i]:.4f}  I*(high)={I_high_row[i]:.4f}",
            flush=True,
        )
    return j, I_low_row, I_high_row


def run_sweep(w_func, w_args, state_meta, n_workers, label, is_omega_constant=False):
    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"  Grid: {N_LAM} lambda x {N_NU} nu")
    print(f"  T_MAX={T_MAX}  TRAJ_POINTS={TRAJ_POINTS}  workers={n_workers}")
    print(f"{'='*60}\n")
    t0 = time.time()

    tasks = [
        (j, nu, w_func, w_args, state_meta, is_omega_constant)
        for j, nu in enumerate(NU_GRID)
    ]

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
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--alpha", type=float,
                      help="Diversity-tension amplitude α")
    mode.add_argument("--alpha-index", type=int,
                      help=f"Index into ALPHA_VALUES list: {ALPHA_VALUES}")
    mode.add_argument("--baseline", action="store_true",
                      help=f"Run constant-omega baseline (omega={OMEGA_SCALAR})")
    mode.add_argument("--omega", type=float,
                      help="Run constant-omega sweep at a custom omega value")
    parser.add_argument("--workers", type=int, default=None,
                        help="Parallel workers (default: all available CPUs)")
    args = parser.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)

    gm, pn, mmax, nmax, state_meta = load_group_statistics(NETWORK)
    n_workers = args.workers or cpu_count()
    date_str = time.strftime('%Y-%m-%d')

    if args.baseline or args.omega is not None:
        omega = OMEGA_SCALAR if args.baseline else args.omega
        suffix = "baseline" if args.baseline else f"omega_{omega:.4f}"
        outfile = os.path.join(OUT_DIR, f"{date_str}_{NETWORK}_{suffix}.npz")
        I_low, I_high, delta = run_sweep(
            w_constant, (omega,), state_meta, n_workers,
            f"Constant-omega  omega={omega}",
            is_omega_constant=True,
        )
        np.savez_compressed(
            outfile,
            lam_grid=LAM_GRID, nu_grid=NU_GRID,
            I_low=I_low, I_high=I_high, delta=delta,
            omega_scalar=omega, mu=MU,
            I0_low=I0_LOW, I0_high=I0_HIGH,
        )

    else:
        alpha = args.alpha if args.alpha is not None else ALPHA_VALUES[args.alpha_index]
        outfile = os.path.join(OUT_DIR, f"{date_str}_{NETWORK}_kernel_alpha_{alpha:.1f}.npz")
        I_low, I_high, delta = run_sweep(
            w_diversity_tension, (alpha,), state_meta, n_workers,
            f"Kernel  alpha={alpha}",
        )
        np.savez_compressed(
            outfile,
            lam_grid=LAM_GRID, nu_grid=NU_GRID,
            I_low=I_low, I_high=I_high, delta=delta,
            alpha=alpha, mu=MU,
            I0_low=I0_LOW, I0_high=I0_HIGH,
        )

    print(f"\nSaved → {outfile}")


if __name__ == "__main__":
    main()
