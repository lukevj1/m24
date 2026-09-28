"""Forecasts built on a posterior: the standardised 12-hour level and profiles.

The standardised 12-hour serum lithium (Li12) is the number every guideline is
written in: the steady-state concentration 12 h after the evening dose. Real
blood tests are drawn when the patient can get to the collection centre, often
days after a dose change and at the wrong time of day. The engine answers the
question clinicians actually have: *given everything we know, what is this
person's Li12 on their current regimen, and how sure are we?*
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .bayes import Posterior
from .case import Administration, li12_clock, regimen_cycle
from .pk import steady_state


@dataclass
class Summary:
    median: float
    lo: float            # 5th percentile
    hi: float            # 95th percentile
    samples: np.ndarray
    unit: str = "mmol/L"
    digits: int = 2

    def prob_between(self, a: float, b: float) -> float:
        return float(np.mean((self.samples >= a) & (self.samples <= b)))

    def prob_above(self, x: float) -> float:
        return float(np.mean(self.samples > x))

    def prob_below(self, x: float) -> float:
        return float(np.mean(self.samples < x))

    def __str__(self) -> str:
        d = self.digits
        return f"{self.median:.{d}f} (90% CrI {self.lo:.{d}f}-{self.hi:.{d}f}) {self.unit}"


def summarise(samples: np.ndarray, unit: str = "mmol/L", digits: int = 2) -> Summary:
    lo, med, hi = np.percentile(samples, [5, 50, 95])
    return Summary(float(med), float(lo), float(hi), np.asarray(samples), unit, digits)


def li12_per_unit(post: Posterior, admins: Sequence[Administration], t: float, n: int = 2000,
                  rng: np.random.Generator | None = None, cl_mult: float | np.ndarray = 1.0,
                  sir: bool = False) -> np.ndarray:
    """Posterior samples of steady-state Li12 per 1 mmol/day of the regimen's
    shape. Kinetics are linear, so Li12 for any daily dose is this times the
    dose; that makes searching over doses essentially free.

    ``cl_mult`` may be an array of per-sample multipliers (e.g. an uncertain
    interaction effect)."""
    cycle = regimen_cycle(admins)
    total = sum(a for _, a, _ in cycle)
    unit_cycle = [(c, a / total, f) for c, a, f in cycle]
    clock = li12_clock(admins)
    draws = post.draw(n, t, rng, sir=sir)
    mult = np.broadcast_to(np.asarray(cl_mult, dtype=float), (n,))
    out = np.empty(n)
    for i, ((eta, drift), m) in enumerate(zip(draws, mult)):
        p = post.params_at(eta, drift, t, cl_mult=m)
        out[i] = steady_state(unit_cycle, p, clock)[0]
    return out


def standardized_li12(post: Posterior, admins: Sequence[Administration], t: float | None = None,
                      n: int = 2000, rng: np.random.Generator | None = None, sir: bool = False) -> Summary:
    """The patient's Li12 on ``admins`` at time ``t`` (default: last level)."""
    if t is None:
        t = max([lv.time for lv in post.case.levels], default=0.0)
    daily = sum(a.mmol for a in admins)
    return summarise(daily * li12_per_unit(post, admins, t, n, rng, sir=sir))


def half_life_h(post: Posterior, t: float, n: int = 500, rng: np.random.Generator | None = None) -> Summary:
    draws = post.draw(n, t, rng)
    return summarise(np.array([post.params_at(e, d, t).terminal_half_life() for e, d in draws]), "h", 0)


def profile(post: Posterior, times: np.ndarray, n: int = 300, rng: np.random.Generator | None = None,
            extra_doses=None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Median and 90% band of the concentration-time curve over ``times``
    (history plus any ``extra_doses`` planned for the future)."""
    rng = rng or np.random.default_rng(1)
    U = post.sample(n, rng)
    doses = list(post.case.doses) + list(extra_doses or [])
    curves = np.array([post.predict(times, u, doses=doses) for u in U])
    lo, med, hi = np.percentile(curves, [5, 50, 95], axis=0)
    return med, lo, hi
