"""Monitoring schedule: guideline floor + uncertainty-driven tightening.

Fixed calendars ("every 3 months, then every 6") treat a 30-year-old with a
rock-steady level of 0.65 and an 80-year-old at 0.95 with falling eGFR the
same. The engine keeps the guideline calendar as a *floor* - it never
recommends monitoring less often than the configured guideline - and then
brings tests forward when the forecast says the patient is likely to drift out
of the safe band before the next routine test.

The guideline rules are data (``RULES``) so services can configure the local
standard (NICE, RANZCP/Australian practice, CANMAT/ISBD, local policy).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np

from .bayes import Posterior
from .case import Administration
from .forecast import li12_per_unit

DAY = 24.0
MONTH = 30.4375 * DAY

# Guideline calendars (months unless stated). Values reflect NICE CG185 as
# commonly applied; see docs/EVIDENCE.md. Other guidelines differ mainly in
# the maintenance interval and in how "high risk" is defined.
RULES = {
    "NICE": {
        "level_after_change_days": 7,
        "level_first_year_months": 3,
        "level_maintenance_months": 6,
        "level_high_risk_months": 3,
        "renal_thyroid_calcium_weight_months": 6,
        "high_risk_level": 0.8,
    },
    "RANZCP": {
        "level_after_change_days": 7,
        "level_first_year_months": 3,
        "level_maintenance_months": 6,
        "level_high_risk_months": 3,
        "renal_thyroid_calcium_weight_months": 6,
        "high_risk_level": 0.8,
    },
}


@dataclass
class Task:
    due_h: float
    test: str
    reason: str
    priority: str = "routine"   # routine | soon | urgent


@dataclass
class Context:
    """Facts the calendar needs; most are derived automatically in a real system."""

    start_h: float                       # when lithium was started
    last_change_h: float                 # last dose change
    stable: bool                         # levels stable on the current dose
    age: float
    interacting_drugs: bool = False
    renal_or_thyroid_risk: bool = False
    raised_calcium: bool = False
    poor_symptom_control: bool = False
    poor_adherence: bool = False
    last_level: float | None = None
    last_level_h: float | None = None
    last_bloods_h: float | None = None   # U&E/eGFR, TFT, calcium, weight
    pregnant_weeks: float | None = None  # gestational age now, if pregnant
    postpartum_days: float | None = None

    def high_risk_reasons(self, threshold: float) -> list[str]:
        reasons = []
        if self.age >= 65:
            reasons.append("older adult")
        if self.interacting_drugs:
            reasons.append("interacting medicines")
        if self.renal_or_thyroid_risk:
            reasons.append("renal/thyroid risk")
        if self.raised_calcium:
            reasons.append("raised calcium")
        if self.poor_symptom_control:
            reasons.append("poor symptom control")
        if self.poor_adherence:
            reasons.append("adherence concerns")
        if self.last_level is not None and self.last_level >= threshold:
            reasons.append(f"last level >= {threshold}")
        return reasons


def guideline_tasks(ctx: Context, now_h: float, guideline: str = "NICE") -> list[Task]:
    r = RULES[guideline]
    tasks: list[Task] = []
    last_level_h = ctx.last_level_h if ctx.last_level_h is not None else ctx.start_h

    if ctx.pregnant_weeks is not None:
        interval = 7 * DAY if ctx.pregnant_weeks >= 36 else 28 * DAY
        tasks.append(Task(max(now_h, last_level_h + interval), "lithium level",
                          "pregnancy: every 4 weeks, weekly from 36 weeks (clearance rises then falls abruptly at delivery)",
                          "soon"))
    elif ctx.postpartum_days is not None and ctx.postpartum_days < 28:
        tasks.append(Task(now_h + DAY if ctx.postpartum_days < 1 else now_h + 7 * DAY, "lithium level",
                          "postpartum: within 24 h of delivery, then weekly for the first month", "soon"))
    elif not ctx.stable or now_h - ctx.last_change_h < 14 * DAY:
        due = max(ctx.last_change_h + r["level_after_change_days"] * DAY, now_h)
        tasks.append(Task(due, "lithium level", "titration: after each dose change until stable", "soon"))
    else:
        first_year = now_h - ctx.start_h < 12 * MONTH
        risks = ctx.high_risk_reasons(r["high_risk_level"])
        if first_year:
            months, why = r["level_first_year_months"], "first year of treatment"
        elif risks:
            months, why = r["level_high_risk_months"], "higher-risk group: " + ", ".join(risks)
        else:
            months, why = r["level_maintenance_months"], "stable maintenance"
        tasks.append(Task(max(now_h, last_level_h + months * MONTH), "lithium level", why))

    bloods_h = ctx.last_bloods_h if ctx.last_bloods_h is not None else ctx.start_h
    m = r["renal_thyroid_calcium_weight_months"]
    if ctx.renal_or_thyroid_risk or ctx.raised_calcium:
        m = min(m, 3)
    tasks.append(Task(max(now_h, bloods_h + m * MONTH), "U&E + eGFR, TSH, calcium, weight/BMI",
                      f"routine safety bloods every {m} months"))
    return tasks


def uncertainty_horizon(post: Posterior, admins: Sequence[Administration], now_h: float, *,
                        band: tuple[float, float] = (0.4, 1.0), p_max: float = 0.10,
                        egfr_slope_per_year: float = 0.0, egfr_now: float | None = None,
                        horizon_months: int = 12, n: int = 1500,
                        rng: np.random.Generator | None = None) -> tuple[float, list[tuple[float, float]]]:
    """How long until the forecast probability of being outside ``band`` exceeds ``p_max``?

    Uncertainty grows with time because clearance drifts (OU process) and, if a
    renal trend is supplied, because clearance follows eGFR down. Returns the
    horizon in hours and the monthly risk curve.
    """
    rng = rng or np.random.default_rng(5)
    daily = sum(a.mmol for a in admins)
    curve = []
    horizon = horizon_months * MONTH
    for m in range(0, horizon_months + 1):
        t = now_h + m * MONTH
        mult = 1.0
        if egfr_now and egfr_slope_per_year:
            future = max(5.0, egfr_now + egfr_slope_per_year * m / 12.0)
            mult = future / egfr_now
        s = daily * li12_per_unit(post, admins, t, n, rng, cl_mult=mult)
        p_out = float(np.mean((s < band[0]) | (s > band[1])))
        curve.append((m, p_out))
        if p_out > p_max and horizon == horizon_months * MONTH:
            horizon = max(1, m) * MONTH if m > 0 else 0.0
    return horizon, curve


@dataclass
class Plan:
    tasks: list[Task]
    horizon_h: float | None = None
    risk_curve: list[tuple[float, float]] = field(default_factory=list)


def plan(ctx: Context, now_h: float, *, post: Posterior | None = None,
         admins: Sequence[Administration] | None = None, guideline: str = "NICE",
         egfr_slope_per_year: float = 0.0, egfr_now: float | None = None) -> Plan:
    tasks = guideline_tasks(ctx, now_h, guideline)
    horizon, curve = None, []
    if post is not None and admins and ctx.stable and ctx.pregnant_weeks is None:
        horizon, curve = uncertainty_horizon(post, admins, now_h, egfr_slope_per_year=egfr_slope_per_year,
                                             egfr_now=egfr_now)
        for tk in tasks:
            if tk.test == "lithium level" and now_h + horizon < tk.due_h:
                tk.due_h = now_h + horizon
                if horizon <= 0:
                    tk.reason += "; DUE NOW: forecast risk of being outside 0.4-1.0 mmol/L already exceeds 10%"
                else:
                    tk.reason += (f"; brought forward: forecast risk of leaving 0.4-1.0 mmol/L exceeds 10% "
                                  f"within ~{horizon / MONTH:.0f} month(s)")
                tk.priority = "soon"
    return Plan(sorted(tasks, key=lambda t: t.due_h), horizon, curve)
