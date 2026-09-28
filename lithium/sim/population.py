"""Virtual patients for in-silico trials.

Each virtual patient has covariates drawn from a population resembling an
adult bipolar-disorder clinic and "true" pharmacokinetics drawn from a TRUTH
model that deliberately differs from the engine's prior (different structure,
covariate effects and variability). The engine never sees the truth; it sees
only what a clinic would: doses as prescribed, levels with assay noise, and
sample times as recorded.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from lithos.case import CovariateTimeline
from lithos.models import Covariates, Formulation, PopModel
from lithos.pk import Dose, Params, simulate
from lithos.renal import egfr_ckd_epi_2021


def _scr_for_egfr(target: float, age: float, sex: str) -> float:
    lo, hi = 20.0, 1500.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if egfr_ckd_epi_2021(mid, age, sex) > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


@dataclass
class VirtualPatient:
    pid: int
    cov: CovariateTimeline
    eta: dict[str, float]
    drift_times: np.ndarray       # h, start of each daily drift segment
    drift_values: np.ndarray      # log-CL drift per segment
    miss_prob: float

    def schedule(self, truth: PopModel, cl_mult_fn=None) -> list[tuple[float, Params]]:
        c0 = self.cov(0.0)
        out = []
        for t, d in zip(self.drift_times, self.drift_values):
            m = cl_mult_fn(t) if cl_mult_fn else 1.0
            out.append((t, truth.individual(c0, self.eta, float(d), cl_mult=m)))
        out[0] = (-np.inf, out[0][1])
        return out

    def take(self, doses: list[Dose], rng: np.random.Generator) -> list[Dose]:
        """Doses actually swallowed (random misses)."""
        return [d for d in doses if rng.random() >= self.miss_prob]


def make_population(n: int, truth: PopModel, seed: int = 1, horizon_days: int = 400,
                    older_fraction: float = 0.2) -> list[VirtualPatient]:
    rng = np.random.default_rng(seed)
    pts = []
    names = truth.eta_names
    Om = truth.omega_matrix()
    for i in range(n):
        older = rng.random() < older_fraction
        age = float(np.clip(rng.normal(74, 5) if older else rng.normal(42, 12), 18, 90))
        sex = "F" if rng.random() < 0.55 else "M"
        height = float(rng.normal(163, 7) if sex == "F" else rng.normal(177, 7))
        bmi = float(np.clip(rng.normal(29, 5.5), 18, 48))            # higher BMI is common in this population
        weight = bmi * (height / 100) ** 2
        egfr = float(np.clip(rng.normal(112 - 0.85 * (age - 20), 14), 35, 135))
        scr = _scr_for_egfr(egfr, age, sex)
        cov = CovariateTimeline(sex=sex, age_at_t0=age, weight=weight, height=height, creatinine=[(0.0, scr)])
        eta_v = rng.multivariate_normal(np.zeros(len(names)), Om)
        eta = dict(zip(names, eta_v))
        # Drift from the truth's two OU processes, simulated exactly on a daily
        # grid. (A weekly grid made clearance jump once a week, which penalised
        # any estimate made just after a jump.)
        days = int(horizon_days) + 2
        times = np.arange(days) * 24.0
        drift = np.zeros(days)
        for sd, tau in (truth.drift_slow, truth.drift_fast):
            rho = np.exp(-1.0 / tau)
            x = rng.normal(0, sd)
            path = []
            for _ in range(days):
                path.append(x)
                x = rho * x + rng.normal(0, sd * np.sqrt(1 - rho * rho))
            drift += np.array(path)
        miss = 0.03 if rng.random() < 0.8 else 0.15                 # most adherent, some partially adherent
        pts.append(VirtualPatient(i, cov, eta, times, drift, miss))
    return pts


def true_level(vp: VirtualPatient, truth: PopModel, doses_taken: list[Dose], t: float, cl_mult_fn=None) -> float:
    return float(simulate(doses_taken, vp.schedule(truth, cl_mult_fn), [t])[0])


def true_li12(vp: VirtualPatient, truth: PopModel, cycle, t: float, cl_mult: float = 1.0) -> float:
    """True steady-state 12-h level on a regimen at time t (full adherence)."""
    from lithos.pk import steady_state
    k = int(np.searchsorted(vp.drift_times, t, side="right")) - 1
    p = truth.individual(vp.cov(0.0), vp.eta, float(vp.drift_values[max(k, 0)]), cl_mult=cl_mult)
    clock = (max(c for c, _, _ in cycle) + 12.0) % 24.0
    return float(steady_state(cycle, p, clock)[0])


def measure(value: float, rng: np.random.Generator, cv: float = 0.04, sd_add: float = 0.02) -> float:
    return round(max(0.01, value * (1 + cv * rng.standard_normal()) + sd_add * rng.standard_normal()), 2)
