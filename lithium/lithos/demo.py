"""Run: python -m lithos.demo

Five everyday situations that make clinicians nervous about lithium, each
played through the engine. Patients are simulated from hidden "true"
parameters so the engine's answers can be checked against the truth.
"""

from __future__ import annotations

from datetime import datetime

import numpy as np

from .bayes import individualise
from .case import Administration, Case, CovariateTimeline, Level, expand, li12_clock, regimen_cycle
from .dosing import default_target, initiation_plan, recommend
from .forecast import standardized_li12
from .interactions import what_if
from .models import TRUTH_PERTURBED, get_engine
from .pk import simulate, steady_state
from .renal import egfr_trend
from .report import card, patient_message
from .safety import triage
from .schedule import Context, plan
from .units import PRODUCTS

T0 = datetime(2026, 9, 1)  # a Tuesday at midnight
LITHICARB = PRODUCTS["Lithicarb 250 mg"]
RNG = np.random.default_rng(2026)


def _true_params(cov, cl_factor, v_factor):
    """Hidden 'truth' for a simulated patient, from a different model than the engine uses."""
    return TRUTH_PERTURBED.individual(cov, {"cl": np.log(cl_factor), "v1": np.log(v_factor)})


def _observe(doses, truth, t, cv=0.04):
    return round(float(simulate(doses, truth, [t])[0] * (1 + cv * RNG.standard_normal())), 2)


def _true_li12(truth, admins):
    return float(steady_state(regimen_cycle(admins), truth, li12_clock(admins))[0])


def _naive(obs, target=(0.6, 0.8)):
    """What reading the number at face value, as if it were a 12-h level, would suggest."""
    if obs > target[1]:
        return f"read at face value, {obs:.2f} looks ABOVE target -> 'reduce the dose'"
    if obs < target[0]:
        return f"read at face value, {obs:.2f} looks BELOW target -> 'increase the dose'"
    return f"read at face value, {obs:.2f} looks in range -> 'no change'"


def mistimed(model):
    """Twice-daily patient takes the morning dose, then has blood taken 2 h later."""
    cov = CovariateTimeline(sex="F", age_at_t0=58, weight=72, height=163, creatinine=[(0, 82)])
    admins = [Administration(8.0, 500), Administration(20.0, 500)]
    doses = expand(8.0, 24 * 40 + 8.1, admins)       # includes the 08:00 dose on the test day
    truth = _true_params(cov(0), 0.95, 1.05)
    t_draw = 24 * 40 + 10.0
    obs = _observe(doses, truth, t_draw)
    case = Case(cov, doses, [Level(t_draw, obs, timing_sd=0.25)], t0=T0)
    post = individualise(case, model)
    li12 = standardized_li12(post, admins)
    rec = recommend(post, LITHICARB, template="bd", current=admins)
    title = (f"1. MORNING DOSE TAKEN 2 H BEFORE THE TEST - {_naive(obs)}."
             f"  [simulated truth: 12-h level {_true_li12(truth, admins):.2f}]")
    return card(post, admins, li12, title=title, rec=rec)


def initiation(model):
    """New patient: first dose from the model, one level on day 4 at a convenient time."""
    cov = CovariateTimeline(sex="M", age_at_t0=34, weight=86, height=181, creatinine=[(0, 88)])
    prior = individualise(Case(cov, [], [], t0=T0), model)
    ip = initiation_plan(prior, LITHICARB)
    admins = ip.maintenance.admins
    doses = expand(21.0, 24 * ip.lead_in_days + 21.0, ip.lead_in)
    doses += expand(24 * ip.lead_in_days + 21.0, 24 * 8 + 17.0, admins)
    truth = _true_params(cov(0), 0.75, 0.9)     # clears lithium more slowly than a typical man his age
    t_draw = 24 * 8 + 16.5                               # day 8, 16:30: 19.5 h after the dose, not steady state
    obs = _observe(doses, truth, t_draw)
    case = Case(cov, doses, [Level(t_draw, obs)], t0=T0)
    post = individualise(case, model)
    li12 = standardized_li12(post, admins)
    rec = recommend(post, LITHICARB, template="nocte", current=admins)
    title = (f"2. STARTING LITHIUM - plan from the population model: {ip.describe(LITHICARB)}\n"
             f"   Level on day 8 at 16:30 (4 days on the new dose, not at steady state): {obs:.2f}.  "
             f"[simulated truth: 12-h level on this dose {_true_li12(truth, admins):.2f}]")
    return card(post, admins, li12, title=title, rec=rec), (rec, li12)


def interaction(model):
    """Stable patient; the GP is about to prescribe hydrochlorothiazide."""
    cov = CovariateTimeline(sex="F", age_at_t0=62, weight=68, height=160, creatinine=[(0, 74)])
    admins = [Administration(21.0, 750)]
    doses = expand(21.0, 24 * 200, admins)
    truth = _true_params(cov(0), 1.0, 1.0)
    levels = [Level(24 * d + 9.5, _observe(doses, truth, 24 * d + 9.5)) for d in (60, 150)]
    case = Case(cov, doses, levels, t0=T0)
    post = individualise(case, model)
    li12 = standardized_li12(post, admins)
    wi = what_if(post, admins, "hydrochlorothiazide", product=LITHICARB)
    title = "3. A NEW BLOOD-PRESSURE TABLET - alert fires in the GP's software at the moment of prescribing"
    return card(post, admins, li12, title=title, whatifs=[wi])


def renal(model):
    """Long-term patient whose creatinine is creeping up: trend + adaptive recall."""
    creat = [(0, 92), (24 * 180, 95), (24 * 365, 101), (24 * 545, 104), (24 * 730, 110), (24 * 910, 117),
             (24 * 1095, 124)]
    cov = CovariateTimeline(sex="M", age_at_t0=68, weight=80, height=175, creatinine=creat)
    admins = [Administration(21.0, 750)]
    doses = expand(24 * 900 + 21.0, 24 * 1100, admins)
    truth_t = [(-np.inf, _true_params(cov(24 * 900), 1.0, 1.0)),
               (24 * 1000.0, _true_params(cov(24 * 1095), 1.0, 1.0))]
    levels = [Level(24 * d + 9.0, round(float(simulate(doses, truth_t, [24 * d + 9.0])[0]), 2)) for d in (960, 1090)]
    case = Case(cov, doses, levels, t0=T0)
    post = individualise(case, model)
    li12 = standardized_li12(post, admins)
    yrs = np.array([t for t, _ in creat]) / (24 * 365.25)
    egfrs = np.array([cov(t).egfr for t, _ in creat])
    rt = egfr_trend(yrs, egfrs)
    now = 24 * 1095.0
    ctx = Context(start_h=0.0, last_change_h=24 * 800.0, stable=True, age=cov(now).age,
                  renal_or_thyroid_risk=bool(rt.flags), last_level=levels[-1].value,
                  last_level_h=levels[-1].time, last_bloods_h=now)
    pl = plan(ctx, now, post=post, admins=admins, egfr_slope_per_year=rt.slope, egfr_now=rt.latest)
    target = default_target(cov(now).age)
    rec = recommend(post, LITHICARB, target=target, current=admins)
    title = ("4. SLOWLY FALLING KIDNEY FUNCTION - trend detection, an age-appropriate target (0.4-0.6) "
             "and a prompt for the renal conversation")
    return card(post, admins, li12, title=title, target=target, rec=rec, renal=rt, plan=pl)


def toxicity(model):
    """Gastroenteritis for three days; level 1.9 mmol/L."""
    cov = CovariateTimeline(sex="F", age_at_t0=67, weight=64, height=158, creatinine=[(0, 78), (24 * 63, 128)])
    admins = [Administration(21.0, 750)]
    doses = expand(21.0, 24 * 62 + 22, admins)
    truth = [(-np.inf, _true_params(cov(0), 1.0, 1.0)),
             (24 * 59.0, _true_params(cov(24 * 63), 0.8, 0.9))]   # dehydration: clearance falls
    levels = [Level(24 * 30 + 9, round(float(simulate(doses, truth, [24 * 30 + 9])[0]), 2)),
              Level(24 * 63 + 11, round(float(simulate(doses, truth, [24 * 63 + 11])[0]), 2))]
    case = Case(cov, doses, levels, t0=T0)
    post = individualise(case, model)
    li12 = standardized_li12(post, admins)
    tr = triage(levels[-1].value, {"diarrhoea", "coarse tremor"}, egfr=cov(24 * 63).egfr, post=post,
                t=levels[-1].time)
    title = f"5. GASTROENTERITIS - level {levels[-1].value:.2f} with coarse tremor"
    return card(post, admins, li12, title=title, triage=tr)


def main():
    model = get_engine()
    print(mistimed(model), end="\n\n")
    text, (rec2, li12_2) = initiation(model)
    print(text, end="\n\n")
    for scenario in (interaction, renal, toxicity):
        print(scenario(model), end="\n\n")
    days = rec2.recheck_days if rec2.change != "no change" else 7   # confirm stability (NICE: weekly until stable)
    print("Patient message after scenario 2:\n  " + patient_message(rec2, li12_2, (0.6, 0.8), days))


if __name__ == "__main__":
    main()
