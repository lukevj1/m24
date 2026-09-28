import numpy as np
import pytest

from lithos import interactions
from lithos.bayes import fit
from lithos.case import Administration, Case, CovariateTimeline, Level, expand
from lithos.dosing import initiation_plan, recommend
from lithos.forecast import standardized_li12
from lithos.models import get_model
from lithos.pk import simulate
from lithos.renal import crcl_cockcroft_gault, egfr_ckd_epi_2021, egfr_trend
from lithos.safety import hours_to_below, triage
from lithos.schedule import DAY, MONTH, Context, guideline_tasks, plan
from lithos.units import PRODUCTS, mg_to_mmol, mmol_to_mg

MODEL = get_model()
LITHICARB = PRODUCTS["Lithicarb 250 mg"]


def cov_typical(**kw):
    base = dict(sex="F", age_at_t0=45, weight=70, height=165, creatinine=[(0, 75)])
    base.update(kw)
    return CovariateTimeline(**base)


# --- units -----------------------------------------------------------------

def test_salt_conversions():
    assert mg_to_mmol(250) == pytest.approx(6.767, abs=1e-3)
    assert mg_to_mmol(400) == pytest.approx(10.83, abs=1e-2)
    assert mg_to_mmol(520, "citrate") == pytest.approx(5.53, abs=1e-2)
    assert mmol_to_mg(mg_to_mmol(450)) == pytest.approx(450)


# --- renal -----------------------------------------------------------------

def test_ckd_epi_2021_reference_value():
    # 60-year-old man, creatinine 1.0 mg/dL (88.4 umol/L): eGFR 86 by CKD-EPI 2021
    assert egfr_ckd_epi_2021(88.42, 60, "M") == pytest.approx(86, abs=0.6)


def test_cockcroft_gault_reference_value():
    assert crcl_cockcroft_gault(88.42, 60, "M", 70) == pytest.approx(77.8, abs=0.1)


def test_egfr_trend_flags_rapid_decline_only_when_present():
    yrs = np.arange(0, 4.01, 0.5)
    fast = egfr_trend(yrs, 80 - 7 * yrs + np.array([1, -1, 0.5, -0.5, 0, 1, -1, 0.3, 0]))
    assert any("rapid" in f for f in fast.flags)
    stable = egfr_trend(yrs, 88 - 0.8 * yrs)
    assert not any("rapid" in f for f in stable.flags)


# --- Bayesian fit -------------------------------------------------------------

def test_prior_only_fit_is_the_prior():
    post = fit(MODEL, Case(cov_typical(), expand(21, 24 * 10, [Administration(21, 750)]), []))
    assert np.allclose(post.mode, 0.0)
    assert post.log_evidence == 0.0


def test_map_recovers_clearance_with_rich_data():
    cov = cov_typical()
    admins = [Administration(21, 750)]
    doses = expand(21, 24 * 60, admins)
    true_cl = 0.65
    truth = MODEL.individual(cov(0), {"cl": np.log(true_cl), "v1": np.log(1.1)})
    times = [24 * d + h for d in (7, 14, 28, 42) for h in (9.0, 15.0)]
    levels = [Level(t, float(simulate(doses, truth, [t])[0])) for t in times]
    post = fit(MODEL, Case(cov, doses, levels))
    est_cl = np.exp(post.mode[0] + post.mode[2:].mean())
    assert est_cl == pytest.approx(true_cl, rel=0.12)


def test_uncertainty_shrinks_with_more_levels():
    cov = cov_typical()
    admins = [Administration(21, 750)]
    doses = expand(21, 24 * 40, admins)
    truth = MODEL.individual(cov(0), {"cl": np.log(0.8)})
    one = [Level(24 * 10 + 9, float(simulate(doses, truth, [24 * 10 + 9])[0]))]
    two = one + [Level(24 * 30 + 9, float(simulate(doses, truth, [24 * 30 + 9])[0]))]
    w0 = standardized_li12(fit(MODEL, Case(cov, doses, [])), admins, t=24 * 30 + 9)
    w1 = standardized_li12(fit(MODEL, Case(cov, doses, one)), admins, t=24 * 30 + 9)
    w2 = standardized_li12(fit(MODEL, Case(cov, doses, two)), admins, t=24 * 30 + 9)
    widths = [s.hi - s.lo for s in (w0, w1, w2)]
    assert widths[0] > widths[1] > widths[2]


def test_steep_samples_are_down_weighted():
    """A level drawn 1.5 h post-dose with uncertain timing should move the
    estimate less than the same information drawn on the flat part of the curve."""
    cov = cov_typical()
    admins = [Administration(8, 375), Administration(20, 375)]
    doses = expand(8, 24 * 20 + 9.6, admins)
    truth = MODEL.individual(cov(0), {"cl": np.log(0.7)})
    t_steep, t_flat = 24 * 20 + 9.5, 24 * 20 + 8.0 - 0.01
    steep = Level(t_steep, float(simulate(doses, truth, [t_steep])[0]), timing_sd=1.5)
    flat = Level(t_flat, float(simulate(doses, truth, [t_flat])[0]), timing_sd=1.5)
    sd_steep = np.sqrt(fit(MODEL, Case(cov, doses, [steep])).cov[0, 0])
    sd_flat = np.sqrt(fit(MODEL, Case(cov, doses, [flat])).cov[0, 0])
    assert sd_steep > sd_flat


def test_any_time_level_standardises_close_to_truth():
    rng = np.random.default_rng(3)
    cov = cov_typical()
    admins = [Administration(21, 750)]
    doses = expand(21, 24 * 60, admins)
    errors = []
    for cl_f in (0.7, 1.0, 1.3):
        truth = MODEL.individual(cov(0), {"cl": np.log(cl_f)})
        from lithos.case import li12_clock, regimen_cycle
        from lithos.pk import steady_state
        true_li12 = steady_state(regimen_cycle(admins), truth, li12_clock(admins))[0]
        t = 24 * 30 + 21 + rng.uniform(8, 20)
        lv = Level(t, float(simulate(doses, truth, [t])[0]))
        est = standardized_li12(fit(MODEL, Case(cov, doses, [lv])), admins, n=800)
        errors.append(abs(est.median - true_li12) / true_li12)
        assert est.lo <= true_li12 <= est.hi
    assert max(errors) < 0.08


def test_drift_forecast_uncertainty_grows_with_horizon():
    cov = cov_typical()
    admins = [Administration(21, 750)]
    doses = expand(21, 24 * 40, admins)
    truth = MODEL.individual(cov(0), {})
    lv = [Level(24 * 30 + 9, float(simulate(doses, truth, [24 * 30 + 9])[0]))]
    post = fit(MODEL, Case(cov, doses, lv))
    near = standardized_li12(post, admins, t=24 * 31, n=1500)
    far = standardized_li12(post, admins, t=24 * 30 + 12 * MONTH, n=1500)
    assert (far.hi - far.lo) > (near.hi - near.lo)


# --- dosing ----------------------------------------------------------------

def test_recommendation_hits_target_for_typical_patient():
    cov = cov_typical()
    prior = fit(MODEL, Case(cov, [], []))
    rec = recommend(prior, LITHICARB, target=(0.6, 0.8), t=0.0)
    assert rec.chosen.p_high <= 0.05
    assert 0.45 <= rec.chosen.li12.median <= 0.85
    assert rec.chosen.daily_mg % 250 == 0


def test_bd_split_puts_remainder_at_night():
    cov = cov_typical()
    prior = fit(MODEL, Case(cov, [], []))
    rec = recommend(prior, LITHICARB, template="bd", t=0.0)
    for o in rec.options:
        assert o.admins[1].mg >= o.admins[0].mg


def test_initiation_plan_leads_in_below_maintenance():
    prior = fit(MODEL, Case(cov_typical(), [], []))
    ip = initiation_plan(prior, LITHICARB)
    assert sum(a.mg for a in ip.lead_in) < ip.maintenance.daily_mg
    assert ip.first_level_day > ip.lead_in_days


# --- interactions ------------------------------------------------------------

def test_interaction_lookup_and_direction():
    assert interactions.lookup("Hydrochlorothiazide").key == "thiazide"
    assert interactions.lookup("amlodipine") is None
    cov = cov_typical()
    admins = [Administration(21, 750)]
    doses = expand(21, 24 * 40, admins)
    truth = MODEL.individual(cov(0), {})
    lv = [Level(24 * 30 + 9, float(simulate(doses, truth, [24 * 30 + 9])[0]))]
    post = fit(MODEL, Case(cov, doses, lv))
    up = interactions.what_if(post, admins, "ramipril", n=800, product=LITHICARB)
    down = interactions.what_if(post, admins, "empagliflozin", n=800)
    assert up.after.median > up.before.median
    assert down.after.median < down.before.median
    assert up.adjusted.chosen.daily_mg <= 750
    pd = interactions.what_if(post, admins, "carbamazepine", n=400)
    assert pd.after.median == pytest.approx(pd.before.median)


# --- schedule --------------------------------------------------------------

def test_guideline_calendar_phases():
    titr = guideline_tasks(Context(start_h=0, last_change_h=10 * DAY, stable=False, age=40), 11 * DAY)
    assert titr[0].due_h == pytest.approx(17 * DAY)
    first_year = guideline_tasks(Context(start_h=0, last_change_h=30 * DAY, stable=True, age=40,
                                         last_level_h=60 * DAY), 70 * DAY)
    assert first_year[0].due_h == pytest.approx(60 * DAY + 3 * MONTH)
    maint = guideline_tasks(Context(start_h=0, last_change_h=400 * DAY, stable=True, age=40, last_level=0.65,
                                    last_level_h=800 * DAY), 810 * DAY)
    assert maint[0].due_h == pytest.approx(800 * DAY + 6 * MONTH)
    risky = guideline_tasks(Context(start_h=0, last_change_h=400 * DAY, stable=True, age=70, last_level=0.65,
                                    last_level_h=800 * DAY), 810 * DAY)
    assert risky[0].due_h == pytest.approx(800 * DAY + 3 * MONTH)
    preg = guideline_tasks(Context(start_h=0, last_change_h=0, stable=True, age=30, last_level_h=800 * DAY,
                                   pregnant_weeks=37), 801 * DAY)
    assert preg[0].due_h == pytest.approx(807 * DAY)


def test_uncertainty_horizon_brings_forward_a_patient_near_the_ceiling():
    cov = cov_typical(age_at_t0=50)
    admins = [Administration(21, 1000)]
    doses = expand(21, 24 * 400, admins)
    truth = MODEL.individual(cov(0), {"cl": np.log(0.72)})
    lv = [Level(24 * d + 9, float(simulate(doses, truth, [24 * d + 9])[0])) for d in (200, 380)]
    post = fit(MODEL, Case(cov, doses, lv))
    now = 24 * 381.0
    ctx = Context(start_h=0, last_change_h=24 * 100, stable=True, age=50, last_level=0.75,
                  last_level_h=lv[-1].time, last_bloods_h=now)
    p = plan(ctx, now, post=post, admins=admins)
    level_task = [t for t in p.tasks if t.test == "lithium level"][0]
    assert level_task.due_h < lv[-1].time + 6 * MONTH


# --- safety ----------------------------------------------------------------

def test_triage_levels():
    assert triage(0.7).urgency == "routine"
    assert triage(1.3).urgency == "review"
    assert triage(1.7).urgency == "same-day"
    assert triage(1.6, {"confusion"}).urgency == "emergency"
    assert any("RECOMMENDED" in e for e in triage(4.5, egfr=30).extrip)
    assert any("SUGGESTED" in e for e in triage(5.5).extrip)


def test_hours_to_below_matches_first_order_decay():
    cov = cov_typical()
    admins = [Administration(21, 750)]
    doses = expand(21, 24 * 30, admins)
    truth = MODEL.individual(cov(0), {})
    t = 24 * 29 + 21 + 14
    obs = float(simulate(doses, truth, [t])[0])
    post = fit(MODEL, Case(cov, doses, [Level(t, obs)]))
    s = hours_to_below(post, t, observed=1.8, n=200)
    k = post.params_at(post.eta(post.mode), post.drift_fitted(post.mode, t), t)
    expected = np.log(1.8) * k.v1 / k.cl
    assert s.median == pytest.approx(expected, rel=0.15)
