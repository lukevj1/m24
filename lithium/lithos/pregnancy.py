"""Pregnancy and the postpartum period.

Renal lithium clearance rises through pregnancy as glomerular filtration
increases, then falls back abruptly after delivery. A fixed dose therefore
drifts LOW during pregnancy (relapse risk) and HIGH after delivery (toxicity
risk, in the weeks when relapse risk is also highest). The engine represents
this as a gestation-dependent clearance multiplier and turns it into a
week-by-week dose and monitoring plan that the perinatal team can review.

The default curve is derived from the dose-corrected level changes reported
in observational cohorts (see ``ANCHORS`` and docs/EVIDENCE.md); individual
patients deviate from it, which is why levels are still checked every 4 weeks
and weekly near term, and why every level updates the forecast.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .bayes import Posterior
from .case import Administration, Effect
from .dosing import recommend
from .units import Product

WEEK = 7 * 24.0

# (gestational week, dose-corrected level change vs pre-conception). Filled
# from Wesseloo et al. 2017 (trimester means) - see docs/EVIDENCE.md.
ANCHORS = [(0.0, 0.00), (6.0, -0.10), (10.0, -0.24), (20.0, -0.36), (32.0, -0.30), (39.0, -0.21)]
# Postpartum: level-to-dose ratio rises back above baseline within days of
# delivery, then settles.
POSTPARTUM = [(0.0, +0.09), (2.0, +0.05), (6.0, 0.0)]  # (weeks after delivery, level change)


def cl_multiplier(gest_week: float | None = None, postpartum_week: float | None = None) -> float:
    """Clearance relative to pre-conception (level change x -> CL x 1/(1+x))."""
    if postpartum_week is not None:
        w, x = zip(*POSTPARTUM)
        change = float(np.interp(postpartum_week, w, x))
    else:
        w, x = zip(*ANCHORS)
        change = float(np.interp(gest_week, w, x))
    return 1.0 / (1.0 + change)


def effects(conception_h: float, delivery_h: float, step_weeks: float = 1.0) -> list[Effect]:
    """Weekly clearance effects from conception to 8 weeks postpartum."""
    out = []
    t = conception_h
    while t < delivery_h:
        gw = (t - conception_h) / WEEK
        out.append(Effect(t, min(t + step_weeks * WEEK, delivery_h), cl_multiplier(gest_week=gw),
                          label=f"pregnancy week {gw:.0f}"))
        t += step_weeks * WEEK
    t = delivery_h
    while t < delivery_h + 8 * WEEK:
        pw = (t - delivery_h) / WEEK
        out.append(Effect(t, t + step_weeks * WEEK, cl_multiplier(postpartum_week=pw),
                          label=f"postpartum week {pw:.0f}"))
        t += step_weeks * WEEK
    return out


# Between-woman variability around the population curve (log-scale SD).
CURVE_SD = 0.15


def _mult_samples(m: float, n: int, rng: np.random.Generator) -> np.ndarray:
    if abs(m - 1.0) < 1e-9:
        return np.ones(n)
    return m * np.exp(rng.normal(0.0, CURVE_SD, size=n))


@dataclass
class PregnancyStep:
    gest_week: float | None
    postpartum_week: float | None
    cl_mult: float
    dose: str
    daily_mg: float
    forecast: str


def plan(pre_pregnancy: Posterior, product: Product, current: list[Administration], *,
         target: tuple[float, float] = (0.6, 0.8), weeks: tuple[int, ...] = (0, 8, 12, 16, 20, 24, 28, 32, 36),
         postpartum_weeks: tuple[int, ...] = (0, 1, 2, 4, 6)) -> list[PregnancyStep]:
    """Forecast the dose that keeps the 12-h level in ``target`` through
    pregnancy and after delivery, starting from the pre-pregnancy posterior.

    This is a planning aid for the perinatal team, not an automatic titration:
    each check level (every 4 weeks, weekly from 36 weeks, within 24 h of
    delivery) replaces the population curve with the patient's own data.
    """
    t = max([lv.time for lv in pre_pregnancy.case.levels], default=0.0)
    rng = np.random.default_rng(17)
    n = 1500
    steps = []
    clocks = [a.clock for a in current]
    for gw in weeks:
        m = cl_multiplier(gest_week=gw)
        rec = recommend(pre_pregnancy, product, target=target, template=clocks, current=current, t=t,
                        cl_mult=_mult_samples(m, n, rng), max_step_ratio=2.0, n=n)
        steps.append(PregnancyStep(gw, None, m, rec.chosen.describe(product), rec.chosen.daily_mg,
                                   str(rec.chosen.li12)))
    late = steps[-1]
    for pw in postpartum_weeks:
        m = cl_multiplier(postpartum_week=pw)
        rec = recommend(pre_pregnancy, product, target=target, template=clocks, current=current, t=t,
                        cl_mult=_mult_samples(m, n, rng), max_step_ratio=2.0, n=n)
        steps.append(PregnancyStep(None, pw, m, rec.chosen.describe(product), rec.chosen.daily_mg,
                                   str(rec.chosen.li12)))
    return steps


def postpartum_hazard(pre_pregnancy: Posterior, late_pregnancy_mg: float, current: list[Administration],
                      n: int = 2000) -> "Summary":
    """Forecast 12-h level in the first postpartum week if the late-pregnancy
    dose is simply continued after delivery."""
    from .forecast import li12_per_unit, summarise
    rng = np.random.default_rng(23)
    t = max([lv.time for lv in pre_pregnancy.case.levels], default=0.0)
    total = sum(a.mg for a in current)
    admins = [Administration(a.clock, a.mg * late_pregnancy_mg / total, a.formulation, a.salt) for a in current]
    daily = sum(a.mmol for a in admins)
    m = _mult_samples(cl_multiplier(postpartum_week=0.5), n, rng)
    return summarise(daily * li12_per_unit(pre_pregnancy, admins, t, n, rng, cl_mult=m))
