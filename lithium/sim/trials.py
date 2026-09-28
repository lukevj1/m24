"""In-silico trials of model-informed lithium dosing.

Run: python -m sim.run_all  (from the lithium/ directory)

Experiments
1. Titration: usual care (stepwise or proportional adjustments of a level read
   as if it were a correctly timed 12-h level) vs model-informed dosing.
2. Any-time sampling: how well can a level drawn at a convenient time be
   converted to what a correctly timed 12-h level would have shown?
3. Interaction: a thiazide is started in a stable patient; no action vs the
   engine's pre-emptive dose adjustment.
4. Early sampling: how soon after a dose change is a level informative?

Everything the engine sees is what a clinic would see. The truth model is
deliberately different from the engine's prior, patients miss doses, blood is
drawn at realistic times, and assays are noisy.
"""

from __future__ import annotations

import numpy as np

from lithos.bayes import individualise
from lithos.case import Administration, Case, Level, expand, regimen_cycle
from lithos.dosing import initiation_plan, recommend
from lithos.forecast import standardized_li12
from lithos.interactions import CATALOG, what_if
from lithos.models import PopModel  # noqa: F401  (engine may be a PopModel or an ensemble spec)
from lithos.units import PRODUCTS

from .population import VirtualPatient, measure, true_level, true_li12

PRODUCT = PRODUCTS["Lithicarb 250 mg"]
TARGET = (0.6, 0.8)
CLOCK = 21.0            # once-nightly dosing
MAX_DAYS = 56
N_DRAWS = 1000          # posterior draws per decision (enough for 1-2 % precision on probabilities)


def _round250(mg: float) -> float:
    return float(max(250.0, 250.0 * round(mg / 250.0)))


class Adherence:
    """Decides once, and remembers, which prescribed doses the patient swallowed."""

    def __init__(self, vp: VirtualPatient, rng: np.random.Generator):
        self.vp, self.rng, self.until, self.taken = vp, rng, -np.inf, []

    def update(self, prescribed) -> list:
        new = [d for d in prescribed if d.time > self.until]
        self.taken += self.vp.take(new, self.rng)
        if prescribed:
            self.until = max(self.until, max(d.time for d in prescribed))
        return self.taken


def realistic_hours_after_evening_dose(rng: np.random.Generator) -> float:
    """When people actually get to a collection centre after a 21:00 dose."""
    u = rng.random()
    if u < 0.60:
        return float(rng.uniform(11.0, 13.5))   # morning, roughly as instructed
    if u < 0.90:
        return float(rng.uniform(13.5, 19.0))   # late morning / afternoon
    return float(rng.uniform(19.0, 22.0))       # after work


# ---------------------------------------------------------------------------
# Experiment 1: titration
# ---------------------------------------------------------------------------

def _regimen_doses(periods: list[tuple[float, float, float]]):
    doses = []
    for start, end, mg in periods:
        doses += expand(start, end, [Administration(CLOCK, mg)])
    return doses


def _titration_outcomes(vp, truth, periods, tests, confirmed_day):
    """Truth-side outcomes for a sequence of dose periods."""
    # True steady-state Li12 of each dose actually used (the patient's exposure if kept on it).
    li12s = []
    for start, end, mg in periods:
        cycle = regimen_cycle([Administration(CLOCK, mg)])
        li12s.append(true_li12(vp, truth, cycle, start))
    final_li12 = li12s[-1]
    # First day from which the patient stays on doses whose true Li12 is in range.
    therapeutic_day = None
    for (start, _, _), v in zip(periods, li12s):
        if TARGET[0] <= v <= TARGET[1]:
            therapeutic_day = therapeutic_day if therapeutic_day is not None else start / 24.0
        else:
            therapeutic_day = None
    days_over_1 = sum((end - start) / 24.0 for (start, end, _), v in zip(periods, li12s) if v > 1.0)
    return {
        "tests": tests,
        "confirmed_day": confirmed_day,
        "final_true_li12": final_li12,
        "final_in_range": TARGET[0] <= final_li12 <= TARGET[1],
        "final_in_wide_range": 0.5 <= final_li12 <= 0.9,
        "therapeutic_day": therapeutic_day,
        "days_on_dose_over_1": days_over_1,
        "n_dose_changes": len(periods) - 1,
    }


def usual_care(vp: VirtualPatient, truth: PopModel, rng: np.random.Generator, style: str = "stepwise") -> dict:
    """Start low, level at ~1 week (read as a 12-h level), adjust, repeat until two
    consecutive in-range levels on the same dose."""
    mg = 250.0 if vp.cov(0).age >= 65 else 500.0
    periods = [[CLOCK, None, mg]]
    day, tests, streak, confirmed = 7, 0, 0, None
    adherence = Adherence(vp, rng)
    while day <= MAX_DAYS:
        hrs = realistic_hours_after_evening_dose(rng)
        t = (day - 1) * 24.0 + CLOCK + hrs
        periods[-1][1] = t + 1.0
        prescribed = _regimen_doses([tuple(p) for p in periods])
        obs = measure(true_level(vp, truth, adherence.update(prescribed), t), rng)
        tests += 1
        lo, hi = TARGET
        if lo <= obs <= hi:
            streak += 1
            if streak == 2:
                confirmed = day
                break
            day += 7
            continue
        streak = 0
        if style == "proportional":
            new = _round250(mg * 0.7 / obs)
            new = float(np.clip(new, mg - 500, mg + 500))
        else:
            if obs < 0.4:
                new = mg + 500
            elif obs < lo:
                new = mg + 250
            elif obs > 1.0:
                new = max(250.0, mg - 500)
            else:
                new = max(250.0, mg - 250)
        change_t = day * 24.0 + CLOCK            # result back next day; new dose that evening
        periods[-1][1] = change_t
        periods.append([change_t, None, new])
        mg = new
        day += 8
    end = (min(day, MAX_DAYS) + 1) * 24.0
    periods[-1][1] = max(end, periods[-1][0] + 24.0)
    return _titration_outcomes(vp, truth, [tuple(p) for p in periods], tests, confirmed)


def model_informed(vp: VirtualPatient, truth: PopModel, engine, rng: np.random.Generator,
                   high: float = 1.2) -> dict:
    """Model-chosen start (brief lead-in), levels at convenient times with the
    time recorded, Bayesian update and recommendation after every level.

    The stopping rule mirrors usual care: two consecutive levels on the same
    dose that this arm reads as in range (the engine's standardised 12-h
    estimate is in the target and it recommends no change), at least one of
    them drawn 3.3 or more estimated half-lives after the last dose change, so
    stability is confirmed near steady state, as it is in usual care."""
    cov = vp.cov
    prior = individualise(Case(cov, [], []), engine)
    ip = initiation_plan(prior, PRODUCT, target=TARGET, lead_in_days=3, n=N_DRAWS)
    lead_mg = ip.lead_in[0].mg
    mg = ip.maintenance.daily_mg
    last_change = CLOCK + 3 * 24.0
    periods = [[CLOCK, last_change, lead_mg], [last_change, None, mg]]
    day, tests, confirmed = 7, 0, None
    streak: list[bool] = []     # per in-range level on this dose: drawn near steady state?
    levels: list[Level] = []
    adherence = Adherence(vp, rng)
    while day <= MAX_DAYS:
        hrs = realistic_hours_after_evening_dose(rng)
        t = (day - 1) * 24.0 + CLOCK + hrs
        periods[-1][1] = t + 1.0
        prescribed = _regimen_doses([tuple(p) for p in periods])
        obs = measure(true_level(vp, truth, adherence.update(prescribed), t), rng)
        tests += 1
        levels.append(Level(t, obs, timing_sd=0.25))
        post = individualise(Case(cov, prescribed, levels), engine)
        current = [Administration(CLOCK, mg)]
        rec = recommend(post, PRODUCT, target=TARGET, current=current, n=N_DRAWS, high=high,
                        rng=np.random.default_rng(int(rng.integers(1 << 31))))
        est = rec.current.li12.median
        if rec.change == "no change" and TARGET[0] <= est <= TARGET[1]:
            streak.append(t - last_change >= 3.3 * rec.half_life.median)
            if len(streak) >= 2 and any(streak):
                confirmed = day
                break
            day += 7
            continue
        streak = []
        if rec.repeat_first:
            day += 3            # a better-timed level soon; no dose change
            continue
        if rec.change == "no change":
            day += 7            # out of range but no better option (e.g. the safety ceiling): recheck
            continue
        new = rec.chosen.daily_mg
        change_t = day * 24.0 + CLOCK
        periods[-1][1] = change_t
        periods.append([change_t, None, new])
        mg, last_change = new, change_t
        day += 1 + int(np.clip(rec.recheck_days, 4, 7))
    end = (min(day, MAX_DAYS) + 1) * 24.0
    periods[-1][1] = max(end, periods[-1][0] + 24.0)
    return _titration_outcomes(vp, truth, [tuple(p) for p in periods], tests, confirmed)


# ---------------------------------------------------------------------------
# Experiment 2: any-time sampling
# ---------------------------------------------------------------------------

def any_time(vp: VirtualPatient, truth: PopModel, engine: PopModel, rng: np.random.Generator) -> dict:
    """Stable once-nightly patient, blood drawn 2-22 h after the evening dose.
    Truth = what a correctly timed (12 h) level on the same morning would have shown."""
    per_mg = true_li12(vp, truth, regimen_cycle([Administration(CLOCK, 1000.0)]), 0.0) / 1000.0
    mg = _round250(rng.uniform(0.4, 1.1) / per_mg)
    mg = float(min(mg, 2000.0))
    admins = [Administration(CLOCK, mg)]
    day = 90
    hrs = float(rng.uniform(2.0, 22.0))
    t = (day - 1) * 24.0 + CLOCK + hrs
    t12 = (day - 1) * 24.0 + CLOCK + 12.0
    prescribed = expand(CLOCK, t + 1.0, admins)
    taken = vp.take(prescribed, rng)
    obs = measure(true_level(vp, truth, taken, t), rng)
    truth12 = true_level(vp, truth, taken, t12)

    out = {"hours": hrs, "truth12": truth12, "naive": obs, "mg": mg}
    post = individualise(Case(vp.cov, prescribed, [Level(t, obs, timing_sd=0.25)]), engine)
    pred = post.predict_draws([t12], 400, rng)[:, 0]
    out.update(engine_med=float(np.median(pred)), engine_lo=float(np.percentile(pred, 5)),
               engine_hi=float(np.percentile(pred, 95)))
    # Same, with one earlier level (day 45, drawn at a convenient time) in the record.
    t_hist = 44 * 24.0 + CLOCK + realistic_hours_after_evening_dose(rng)
    hist_obs = measure(true_level(vp, truth, taken, t_hist), rng)
    post2 = individualise(Case(vp.cov, prescribed, [Level(t_hist, hist_obs, 0.25), Level(t, obs, 0.25)]), engine)
    pred2 = post2.predict_draws([t12], 400, rng)[:, 0]
    out.update(hist_med=float(np.median(pred2)), hist_lo=float(np.percentile(pred2, 5)),
               hist_hi=float(np.percentile(pred2, 95)))
    return out


# ---------------------------------------------------------------------------
# Experiment 3: thiazide started in a stable patient
# ---------------------------------------------------------------------------

def thiazide(vp: VirtualPatient, truth: PopModel, engine: PopModel, rng: np.random.Generator,
             true_median: float = 0.70, true_range: tuple[float, float] = (0.50, 0.90)) -> dict:
    per_mg = true_li12(vp, truth, regimen_cycle([Administration(CLOCK, 1000.0)]), 0.0) / 1000.0
    mg = _round250(0.7 / per_mg)
    admins = [Administration(CLOCK, mg)]
    prescribed = expand(CLOCK, 24 * 200.0, admins)
    taken = vp.take(prescribed, rng)
    levels = []
    for d in (60, 150):
        t = (d - 1) * 24.0 + CLOCK + realistic_hours_after_evening_dose(rng)
        levels.append(Level(t, measure(true_level(vp, truth, taken, t), rng), 0.25))
    post = individualise(Case(vp.cov, prescribed, levels), engine)
    now = levels[-1].time
    wi = what_if(post, admins, "hydrochlorothiazide", t=now, n=N_DRAWS, product=PRODUCT,
                 rng=np.random.default_rng(int(rng.integers(1 << 31))))
    sd = (np.log(true_range[1]) - np.log(true_range[0])) / (2 * 1.645)
    m_true = float(np.exp(rng.normal(np.log(true_median), sd)))
    base = true_li12(vp, truth, regimen_cycle(admins), now)
    no_action = true_li12(vp, truth, regimen_cycle(admins), now, cl_mult=m_true)
    adj = wi.adjusted.chosen.admins
    with_engine = true_li12(vp, truth, regimen_cycle(adj), now, cl_mult=m_true)
    return {"base": base, "no_action": no_action, "engine_adjusted": with_engine,
            "pred_med": wi.after.median, "pred_lo": wi.after.lo, "pred_hi": wi.after.hi,
            "mg": mg, "adj_mg": wi.adjusted.chosen.daily_mg}


# ---------------------------------------------------------------------------
# Experiment 4: how early after a dose change is a level informative?
# ---------------------------------------------------------------------------

def early_sampling(vp: VirtualPatient, truth: PopModel, engine: PopModel, rng: np.random.Generator,
                   days=(1, 2, 3, 4, 5, 7)) -> list[dict]:
    per_mg = true_li12(vp, truth, regimen_cycle([Administration(CLOCK, 1000.0)]), 0.0) / 1000.0
    mg_a = _round250(0.5 / per_mg)
    mg_b = mg_a + 250.0 * max(1, round(0.4 * mg_a / 250.0))
    change = 30 * 24.0 + CLOCK
    doses = expand(CLOCK, change, [Administration(CLOCK, mg_a)]) + \
        expand(change, change + 12 * 24.0, [Administration(CLOCK, mg_b)])
    taken = vp.take(doses, rng)
    t_prev = 28 * 24.0 + CLOCK + 12.0
    prev = Level(t_prev, measure(true_level(vp, truth, taken, t_prev), rng), 0.25)
    out = []
    for d in days:
        t = change + 24.0 * (d - 1) + 12.0      # 12 h after the d-th dose on the new regimen
        # truth = steady-state 12-h level on the new dose at the time of sampling
        target_truth = true_li12(vp, truth, regimen_cycle([Administration(CLOCK, mg_b)]), t)
        obs = measure(true_level(vp, truth, taken, t), rng)
        row = {"day": d, "truth": target_truth, "naive": obs}
        for label, lv in (("with_history", [prev, Level(t, obs, 0.25)]), ("no_history", [Level(t, obs, 0.25)])):
            post = individualise(Case(vp.cov, doses, lv), engine)
            s = standardized_li12(post, [Administration(CLOCK, mg_b)], t=t, n=600,
                                  rng=np.random.default_rng(int(rng.integers(1 << 31))))
            row[label] = s.median
            row[label + "_width"] = s.hi - s.lo
        out.append(row)
    return out
