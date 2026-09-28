"""Pregnancy and the postpartum period.

Renal lithium clearance rises through pregnancy as glomerular filtration
increases, then falls back abruptly after delivery. A fixed dose therefore
drifts LOW during pregnancy (relapse risk) and HIGH after delivery (toxicity
risk, in the weeks when relapse risk is also highest). The engine represents
this as a gestation-dependent clearance multiplier and turns it into a
week-by-week dose and monitoring plan that the perinatal team can review.

Two published descriptions of the pregnancy curve disagree and both are
offered (docs/EVIDENCE.md):

* ``westin`` (default): dose-adjusted levels fall log-linearly by 1.2% per
  gestational week (-7% at week 6, -22% at week 20, -34% at week 34), so
  clearance rises ~1.2%/week; doses generally need +50% by the third
  trimester (Westin et al. 2017, BMJ Open, full text).
* ``wesseloo``: trimester means of -24%, -36% (nadir) and -21% vs
  preconception, i.e. a partial third-trimester rebound (Wesseloo et al.
  2017, 1101 levels from 113 pregnancies).

After delivery clearance falls back within days (a steep rise in
dose-adjusted levels by day 4), and levels run ~9-11% above preconception
postpartum (Wesseloo; Imaz 2025). Individual women deviate from any curve,
which is why levels are still checked every 4 weeks, weekly from 36 weeks and
within 24 h of birth, and why every level updates the forecast.

Delivery management is contested and is surfaced, not hidden: US labelling
advises reducing or stopping lithium 2-3 days before expected delivery,
whereas Molenaar et al. (2021, 233 perinatal levels) found no intrapartum
rise in dose-corrected levels and advise against pre-delivery reduction.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .bayes import Posterior
from .case import Administration, Effect
from .dosing import recommend
from .units import Product

WEEK = 7 * 24.0

# Wesseloo 2017 trimester means placed at mid-trimester weeks.
WESSELOO = [(0.0, 0.00), (7.0, -0.24), (20.0, -0.36), (33.0, -0.21), (41.0, -0.21)]
WESTIN_SLOPE = 0.012  # log dose-adjusted level per gestational week (Westin 2017)
# Postpartum level change vs preconception (weeks after delivery). The +9%
# plateau is from Wesseloo/Imaz; its duration beyond the first weeks is an
# assumption.
POSTPARTUM = [(0.0, +0.09), (8.0, +0.09), (12.0, 0.0)]


def cl_multiplier(gest_week: float | None = None, postpartum_week: float | None = None,
                  curve: str = "westin") -> float:
    """Clearance relative to pre-conception (a level change x means CL x 1/(1+x))."""
    if postpartum_week is not None:
        w, x = zip(*POSTPARTUM)
        return 1.0 / (1.0 + float(np.interp(postpartum_week, w, x)))
    if curve == "westin":
        return float(np.exp(WESTIN_SLOPE * max(gest_week, 0.0)))
    w, x = zip(*WESSELOO)
    return 1.0 / (1.0 + float(np.interp(gest_week, w, x)))


def effects(conception_h: float, delivery_h: float, step_weeks: float = 1.0, curve: str = "westin") -> list[Effect]:
    """Weekly clearance effects from conception to 8 weeks postpartum."""
    out = []
    t = conception_h
    while t < delivery_h:
        gw = (t - conception_h) / WEEK
        out.append(Effect(t, min(t + step_weeks * WEEK, delivery_h), cl_multiplier(gest_week=gw, curve=curve),
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


def plan(post: Posterior, product: Product, current: list[Administration], *,
         target: tuple[float, float] = (0.6, 0.8), postpartum_target: tuple[float, float] = (0.8, 1.0),
         weeks: tuple[int, ...] = (0, 8, 12, 16, 20, 24, 28, 32, 36),
         postpartum_weeks: tuple[int, ...] = (0, 1, 2, 4, 6), curve: str = "westin",
         conception_h: float | None = None, delivery_h: float | None = None) -> list[PregnancyStep]:
    """Forecast the dose that keeps the 12-h level in ``target`` through
    pregnancy and after delivery (``postpartum_target`` for the first month).

    Two modes:
    * ``conception_h is None``: ``post`` was fitted before pregnancy; the
      population curve is applied as a clearance multiplier with between-woman
      uncertainty.
    * ``conception_h`` given: ``post`` was fitted on a case that already
      carries the pregnancy effects (``effects()``) and any levels measured
      during pregnancy, so the plan updates from the woman's own data; only
      weeks after the latest level carry extra curve uncertainty.

    A planning aid for the perinatal team, not an automatic titration: each
    check level (every 4 weeks, weekly from 36 weeks, within 24 h of birth)
    updates the forecast.
    """
    t_last = max([lv.time for lv in post.case.levels], default=0.0)
    rng = np.random.default_rng(17)
    n = 1500
    steps = []
    clocks = [a.clock for a in current]

    def forecast(t_eval, mult, tgt):
        return recommend(post, product, target=tgt, template=clocks, current=current, t=t_eval,
                         cl_mult=mult, max_step_ratio=2.0, horizon_days=0.0, n=n)

    for gw in weeks:
        if conception_h is None:
            m = cl_multiplier(gest_week=gw, curve=curve)
            rec = forecast(t_last, _mult_samples(m, n, rng), target)
        else:
            t_eval = conception_h + gw * WEEK
            m = cl_multiplier(gest_week=gw, curve=curve)
            ahead = max(0.0, (t_eval - t_last) / WEEK)
            spread = np.exp(rng.normal(0.0, CURVE_SD * min(1.0, ahead / 8.0), size=n))
            rec = forecast(max(t_eval, t_last), spread, target)
        steps.append(PregnancyStep(gw, None, m, rec.chosen.describe(product), rec.chosen.daily_mg,
                                   str(rec.chosen.li12)))
    for pw in postpartum_weeks:
        m = cl_multiplier(postpartum_week=pw)
        tgt = postpartum_target if pw < 4 else target
        if conception_h is None or delivery_h is None:
            rec = forecast(t_last, _mult_samples(m, n, rng), tgt)
        else:
            t_eval = delivery_h + pw * WEEK + 12.0
            rec = forecast(max(t_eval, t_last), _mult_samples(1.0 + 1e-9, n, rng), tgt)
        steps.append(PregnancyStep(None, pw, m, rec.chosen.describe(product), rec.chosen.daily_mg,
                                   str(rec.chosen.li12)))
    return steps


def delivery_notes() -> list[str]:
    """Points the perinatal plan must address explicitly (docs/EVIDENCE.md E35)."""
    return [
        "Before delivery: practice differs. US labelling advises reducing or stopping lithium 2-3 days before the "
        "expected date; an observational study of 233 perinatal levels found no intrapartum rise and advises against "
        "pre-delivery reduction. Decide with the obstetric and perinatal psychiatry team.",
        "Levels: within 24 h of birth, then twice weekly for 2 weeks; restart or continue at the planned postpartum "
        "dose (first month target often 0.8-1.0 for relapse prevention).",
        "Postpartum analgesia: routine NSAIDs raise lithium levels - prefer paracetamol or plan a level.",
        "Breastfeeding: infant exposure is substantial (infant serum levels about a third to a half of maternal); "
        "a shared decision with paediatric input and infant monitoring if breastfeeding.",
    ]


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
