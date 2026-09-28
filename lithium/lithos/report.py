"""The clinician card: one screen, one decision.

Design rule: lead with the answer (is this person in range, what should the
dose be, when is the next test), then show why, then show the uncertainty.
Everything else is one click deeper.
"""

from __future__ import annotations

from typing import Sequence

from .bayes import Posterior
from .case import Administration
from .dosing import Recommendation
from .forecast import Summary
from .interactions import WhatIf
from .renal import RenalTrend
from .safety import Triage
from .schedule import MONTH, Plan

BANNER = "RESEARCH PROTOTYPE - not a medical device, not for clinical use"


def _status(s: Summary, target: tuple[float, float]) -> str:
    lo, hi = target
    p_in = s.prob_between(lo, hi)
    if s.median < lo:
        where = "BELOW target"
    elif s.median > hi:
        where = "ABOVE target"
    else:
        where = "IN target"
    return f"{where} ({lo}-{hi}); P(in range) {p_in:.0%}"


def card(post: Posterior, admins: Sequence[Administration], li12: Summary, *, title: str,
         target: tuple[float, float] = (0.6, 0.8), rec: Recommendation | None = None,
         whatifs: Sequence[WhatIf] = (), renal: RenalTrend | None = None, plan: Plan | None = None,
         triage: Triage | None = None, now_h: float | None = None) -> str:
    case = post.case
    L: list[str] = []
    rule = "=" * 78
    L += [rule, f"{title}", BANNER, f"Model: {post.model.name} [{post.model.status}]", rule]

    regimen = " + ".join(f"{a.mg:g} mg {a.formulation} at {int(a.clock):02d}:{int(round(a.clock % 1 * 60)):02d}"
                         for a in admins)
    L.append(f"Current regimen: {regimen}  ({sum(a.mg for a in admins):g} mg/day, "
             f"{sum(a.mmol for a in admins):.1f} mmol Li+/day)")

    if case.levels:
        L.append("\nLevels used:")
        for d in post.diagnostics.get("levels", []):
            since = d["hours_since_dose"]
            since_s = f"{since:.1f} h after last dose" if since is not None else "no dose recorded"
            L.append(f"  {case.when(d['time'])}: {d['value']:.2f} mmol/L ({since_s}); "
                     f"model fit {d['fitted']:.2f}")
            for f in d["flags"]:
                L.append(f"    ! {f}")
    else:
        L.append("\nNo levels yet: forecast is from the population model and covariates only.")

    if triage is not None and triage.urgency in ("same-day", "emergency"):
        L.append(f"\nIf this regimen were simply continued, the steady-state 12-h level would reach {li12}.")
    else:
        L.append(f"\nSTANDARDISED 12-HOUR LEVEL (steady state, current regimen): {li12}")
        L.append(f"  -> {_status(li12, target)}; P(> 1.0) {li12.prob_above(1.0):.1%}")
    hist = post.clearance_history()
    if len(hist) >= 2:
        rel = hist[-1][1] / hist[-2][1]
        if abs(rel - 1) > 0.2:
            direction = "lower" if rel < 1 else "higher"
            L.append(f"  ! Clearance now ~{abs(rel - 1):.0%} {direction} than at the previous check. "
                     + ("Look for dehydration, interacting drugs, intercurrent illness or renal decline."
                        if rel < 1 else "Look for missed doses, a formulation change, pregnancy or natriuretic drugs."))

    if rec is not None:
        L.append(f"\nRECOMMENDATION: {rec.change.upper()} -> {rec.chosen.describe(rec.product)}")
        for r in rec.rationale:
            L.append(f"  - {r}")
        for w in rec.warnings:
            L.append(f"  ! {w}")
        alts = [o for o in rec.options if o is not rec.chosen][:2]
        if alts:
            L.append("  Alternatives: " + "; ".join(f"{o.daily_mg:g} mg/day -> {o.li12.median:.2f} "
                                                   f"(P in range {o.p_target:.0%})" for o in alts))
        if rec.change != "no change":
            L.append(f"  Next level: in {rec.recheck_days} days (any time of day, at least 6 h after a dose; "
                     f"record the time of the last dose). Conventional practice would wait "
                     f"{rec.recheck_days_conventional} days for steady state.")

    for w in whatifs:
        L.append(f"\nWHAT IF {w.interaction.label.upper()} IS STARTED (no lithium change):")
        L.append(f"  12-h level {w.before.median:.2f} -> {w.after} ; P(> 1.0) {w.p_above_1:.0%}, "
                 f"P(> 1.2) {w.p_above_1_2:.0%}")
        L.append(f"  Advice: {w.interaction.advice}")
        if w.adjusted is not None:
            a = w.adjusted
            L.append(f"  If it must be started: {a.change} lithium to {a.chosen.describe(a.product)} -> forecast "
                     f"{a.chosen.li12} (P in range {a.chosen.p_target:.0%}, P > 1.0 {a.chosen.p_high:.0%}); "
                     f"level {max(5, a.recheck_days)}-7 days after starting.")

    if renal is not None:
        L.append(f"\nKIDNEYS: latest eGFR {renal.latest:.0f} ({renal.latest_stage}); "
                 + (f"trend {renal.slope:+.1f} mL/min/1.73m2/yr (90% range {renal.slope_lo:+.1f} to "
                    f"{renal.slope_hi:+.1f}) over {renal.years:.1f} y" if renal.n >= 2 and renal.slope == renal.slope
                    else "no trend yet"))
        for f in renal.flags:
            L.append(f"  ! {f}")

    if triage is not None:
        L.append(f"\nTOXICITY TRIAGE: {triage.urgency.upper()}")
        for a in triage.actions:
            L.append(f"  - {a}")
        for e in triage.extrip:
            L.append(f"  * {e}")

    if plan is not None:
        L.append("\nMONITORING PLAN:")
        for tk in plan.tasks:
            when = case.when(tk.due_h)
            L.append(f"  [{tk.priority}] {when}: {tk.test} - {tk.reason}")
        if plan.risk_curve and plan.horizon_h is not None and plan.horizon_h < 12 * MONTH:
            L.append(f"  Uncertainty horizon: {plan.horizon_h / MONTH:.0f} month(s) until the forecast "
                     f"risk of leaving 0.4-1.0 mmol/L exceeds 10%.")
        if plan.risk_curve:
            curve = ", ".join(f"{int(m)}m {p:.0%}" for m, p in plan.risk_curve[:7])
            L.append(f"  Forecast risk of leaving 0.4-1.0 mmol/L if nothing changes: {curve}")
    L.append(rule)
    return "\n".join(L)


def patient_message(rec: Recommendation | None, li12: Summary, target: tuple[float, float],
                    next_test_days: int | None) -> str:
    lo, hi = target
    if li12.median < lo:
        where = "a little below the range we are aiming for"
    elif li12.median > hi:
        where = "a little above the range we are aiming for"
    else:
        where = "right where we want it"
    msg = [f"Your lithium level is {where}."]
    if rec is not None and rec.change not in ("no change",):
        msg.append(f"Your prescriber has approved a new dose: {rec.chosen.describe(rec.product)}.")
    if next_test_days:
        msg.append(f"Please have your next blood test in about {next_test_days} days. Any time of day is fine; "
                   "just tap 'I took my dose' in the app when you take your evening tablets so we know the timing.")
    msg.append("If you get vomiting or diarrhoea, can't keep fluids down, or notice a worse tremor, unsteadiness "
               "or confusion, skip your lithium and contact us the same day.")
    return " ".join(msg)
