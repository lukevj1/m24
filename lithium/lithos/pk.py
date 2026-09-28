"""Linear compartmental pharmacokinetics with piecewise-constant parameters.

Lithium is not metabolised, is cleared almost entirely by the kidney and shows
linear kinetics over the clinical range, so a linear compartment system is an
adequate structural model. This simulator supports

* one- or two-compartment disposition,
* one depot (absorption) compartment per formulation, so immediate- and
  sustained-release products can be mixed in one dosing history,
* parameters that change over time (piecewise constant). This is how the
  engine represents slow drift in clearance, a newly started interacting drug,
  declining renal function or pregnancy, and
* exact superposition through the eigen-decomposition of the rate matrix, so a
  year of twice-daily dosing costs a few vectorised exponentials.

Units: amounts in mmol Li+, volumes in L, clearances in L/h, time in hours,
concentrations in mmol/L.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping, Sequence

import numpy as np

__all__ = ["Dose", "Params", "simulate", "steady_state"]


@dataclass(frozen=True)
class Dose:
    time: float                 # h
    amount: float               # mmol Li+ as administered (bioavailability applied later)
    formulation: str = "IR"


@dataclass
class Params:
    cl: float                               # L/h
    v1: float                               # L, central volume (the only volume if 1-cmt)
    ka: Mapping[str, float]                 # 1/h, per formulation
    f: Mapping[str, float]                  # bioavailability, per formulation
    q: float = 0.0                          # L/h, inter-compartmental clearance (0 -> 1-cmt)
    v2: float = 0.0                         # L, peripheral volume
    tlag: Mapping[str, float] = field(default_factory=dict)  # h, per formulation

    @property
    def two_compartment(self) -> bool:
        return self.q > 0.0 and self.v2 > 0.0

    def terminal_half_life(self) -> float:
        """Terminal elimination half-life in hours."""
        k10 = self.cl / self.v1
        if not self.two_compartment:
            return np.log(2) / k10
        k12, k21 = self.q / self.v1, self.q / self.v2
        s = k10 + k12 + k21
        beta = 0.5 * (s - np.sqrt(s * s - 4 * k10 * k21))
        return np.log(2) / beta


def _rate_matrix(p: Params, forms: Sequence[str]) -> tuple[np.ndarray, int]:
    nd = len(forms)
    n = nd + (2 if p.two_compartment else 1)
    K = np.zeros((n, n))
    c = nd
    for i, name in enumerate(forms):
        ka = p.ka[name]
        K[i, i] = -ka
        K[c, i] = ka
    K[c, c] = -p.cl / p.v1
    if p.two_compartment:
        k12, k21 = p.q / p.v1, p.q / p.v2
        K[c, c] -= k12
        K[c, c + 1] = k21
        K[c + 1, c] = k12
        K[c + 1, c + 1] = -k21
    return K, c


def _eigensystem(p: Params, forms: Sequence[str]):
    """Eigen-decomposition of the rate matrix, nudging absorption rates if the
    system is (nearly) defective, e.g. ka == k. The nudge is < 0.1 %."""
    nudge = 1.0
    for _ in range(6):
        if nudge != 1.0:
            p = Params(p.cl, p.v1, {k: v * nudge ** (i + 1) for i, (k, v) in enumerate(p.ka.items())},
                       p.f, p.q, p.v2, p.tlag)
        K, c = _rate_matrix(p, forms)
        lam, W = np.linalg.eig(K)
        lam, W = lam.real, W.real
        if np.linalg.cond(W) < 1e8:
            return lam, W, np.linalg.inv(W), c
        nudge *= 1.0001
    raise FloatingPointError("could not diagonalise rate matrix")


def _chunks(n: int, size: int = 400):
    for i in range(0, n, size):
        yield slice(i, min(n, i + size))


def simulate(
    doses: Sequence[Dose],
    schedule: Sequence[tuple[float, Params]] | Params,
    times: Sequence[float] | np.ndarray | float,
    *,
    derivative: bool = False,
):
    """Central-compartment concentration (mmol/L) at ``times``.

    ``schedule`` is a list of ``(t_start, Params)`` sorted by ``t_start``;
    parameters are constant from one ``t_start`` to the next and the first entry
    applies from minus infinity. A single ``Params`` means constant parameters.

    With ``derivative=True`` also returns dC/dt (mmol/L/h), which the Bayesian
    layer uses to propagate sample-timing uncertainty into the error model.
    """
    if isinstance(schedule, Params):
        schedule = [(-np.inf, schedule)]
    times = np.atleast_1d(np.asarray(times, dtype=float))
    forms = sorted({d.formulation for d in doses}) or ["IR"]
    fidx = {f: i for i, f in enumerate(forms)}

    lag0 = schedule[0][1].tlag
    arrive = np.array([d.time + lag0.get(d.formulation, 0.0) for d in doses], dtype=float)
    amount = np.array([d.amount for d in doses], dtype=float)
    depot = np.array([fidx[d.formulation] for d in doses], dtype=int)
    order = np.argsort(arrive, kind="stable")
    arrive, amount, depot = arrive[order], amount[order], depot[order]

    conc = np.full(times.shape, np.nan)
    dconc = np.full(times.shape, np.nan)
    starts = [s for s, _ in schedule]
    ends = starts[1:] + [np.inf]
    x = None  # state vector at the start of the current segment (None == all zero)

    for seg, ((a, p), b) in enumerate(zip(schedule, ends)):
        lam, W, Winv, c = _eigensystem(p, forms)
        fvec = np.array([p.f.get(f, 1.0) for f in forms])
        in_seg = (arrive >= a) & (arrive < b) if seg else (arrive < b)
        td = arrive[in_seg]
        H = Winv[:, depot[in_seg]] * (amount[in_seg] * fvec[depot[in_seg]])  # (n, J)
        g = None if x is None else Winv @ x

        obs = (times >= a) & (times < b) if seg else (times < b)
        idx = np.flatnonzero(obs)
        for sl in _chunks(idx.size):
            ii = idx[sl]
            to = times[ii]
            z = np.zeros((to.size, lam.size))
            if g is not None:
                z += np.exp(np.outer(to - a, lam)) * g
            if td.size:
                tau = to[:, None] - td[None, :]
                mask = tau >= 0.0
                E = np.exp(np.where(mask, tau, 0.0)[..., None] * lam) * mask[..., None]
                z += np.einsum("mjn,nj->mn", E, H)
            conc[ii] = z @ W[c] / p.v1
            if derivative:
                dconc[ii] = (z * lam) @ W[c] / p.v1

        if np.isfinite(b):
            zb = np.zeros(lam.size)
            if g is not None:
                zb += np.exp(lam * (b - a)) * g
            if td.size:
                zb += (np.exp(np.outer(b - td, lam)) * H.T).sum(axis=0)
            x = W @ zb

    return (conc, dconc) if derivative else conc


def steady_state(
    regimen: Sequence[tuple[float, float, str]],
    p: Params,
    clock_times: Sequence[float] | np.ndarray | float,
    period: float = 24.0,
) -> np.ndarray:
    """Steady-state concentration for a regimen repeating every ``period`` hours.

    ``regimen`` holds ``(clock_hour, amount_mmol, formulation)`` administrations
    within one period; ``clock_times`` are evaluation times within the period.
    Solved exactly: in eigen-coordinates the one-period map is diagonal, so the
    periodic state is a geometric series.
    """
    clock_times = np.atleast_1d(np.asarray(clock_times, dtype=float)) % period
    forms = sorted({r[2] for r in regimen})
    fidx = {f: i for i, f in enumerate(forms)}
    lam, W, Winv, c = _eigensystem(p, forms)
    s = np.array([(r[0] + p.tlag.get(r[2], 0.0)) % period for r in regimen])
    idx = np.array([fidx[r[2]] for r in regimen])
    amt = np.array([r[1] * p.f.get(r[2], 1.0) for r in regimen])
    H = Winv[:, idx] * amt                                     # (n, J)
    z0 = (np.exp(np.outer(period - s, lam)) * H.T).sum(axis=0) / (1.0 - np.exp(lam * period))
    tau = clock_times[:, None] - s[None, :]
    mask = tau >= 0.0
    E = np.exp(np.where(mask, tau, 0.0)[..., None] * lam) * mask[..., None]
    z = np.exp(np.outer(clock_times, lam)) * z0 + np.einsum("mjn,nj->mn", E, H)
    return z @ W[c] / p.v1
