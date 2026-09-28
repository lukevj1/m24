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
from lithos import pregnancy

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


def test_citrate_pack_mmol_is_authoritative():
    assert PRODUCTS["Li-Liquid 509 mg/5 mL"].mmol == 5.4
    assert PRODUCTS["Lithium citrate oral solution 8 mEq/5 mL"].mmol == 8.0


def test_pregnancy_curves():
    # Westin: -34% dose-adjusted level at week 34 -> clearance x ~1.50
    assert 1.0 / pregnancy.cl_multiplier(gest_week=34) == pytest.approx(0.665, abs=0.01)
    # Wesseloo: second-trimester nadir of -36%
    assert 1.0 / pregnancy.cl_multiplier(gest_week=20, curve="wesseloo") == pytest.approx(0.64, abs=0.01)
    assert pregnancy.cl_multiplier(postpartum_week=1) < 1.0


# --- renal -----------------------------------------------------------------

def test_ckd_epi_2021_reference_value():
    # 60-year-old man, creatinine 1.0 mg/dL (88.4 umol/L): eGFR 86 by CKD-EPI 2021
    assert egfr_ckd_epi_2021(88.42, 60, "M") == pytest.approx(86, abs=0.6)


def test_cockcroft_gault_reference_value():
    assert crcl_cockcroft_gault(88.42, 60, "M", 70) == pytest.approx(77.8, abs=0.1)


def test_egfr_trend_needs_three_values_over_90_days():
    assert "trend needs" in " ".join(egfr_trend(np.array([0.0, 0.1]), np.array([80.0, 78.0])).flags)


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
    assert interactions.lookup("paracetamol") is None
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
    assert triage(2.1).urgency == "emergency"                     # chronic accumulation
    assert triage(1.6, {"confusion"}).urgency == "emergency"
    assert any("RECOMMENDED" in e for e in triage(4.5, egfr=30).extrip)
    assert any("SUGGESTED" in e for e in triage(5.5).extrip)


def test_triage_escalates_fluid_loss_and_kidneys_and_gates_extrip():
    assert triage(1.4, {"vomiting", "diarrhoea"}, egfr=35).urgency == "same-day"
    assert triage(1.3, creatinine_baseline=80, creatinine_now=130).urgency == "same-day"
    # EXTRIP is for lithium poisoning: not raised for symptoms at therapeutic levels
    assert triage(0.5, {"confusion"}).extrip == []
    low_seizure = triage(0.3, {"seizure"})
    assert low_seizure.urgency == "emergency" and low_seizure.extrip == []


def test_hours_to_below_matches_first_order_decay():
    cov = cov_typical()
    admins = [Administration(21, 750)]
    doses = expand(21, 24 * 30, admins)
    truth = MODEL.individual(cov(0), {})
    t = 24 * 29 + 21 + 14
    obs = float(simulate(doses, truth, [t])[0])
    post = fit(MODEL, Case(cov, doses, [Level(t, obs)]))
    s = hours_to_below(post, t, observed=1.8, n=200)
    p = post.params_at(post.eta(post.mode), post.drift_fitted(post.mode, t), t)
    expected = np.log(1.8) * p.terminal_half_life() / np.log(2)   # post-distribution decline
    assert s.median == pytest.approx(expected, rel=0.2)


# --- behaviour added with the ensemble -----------------------------------------

def test_kidney_blind_model_is_excluded_outside_its_population():
    from lithos.bayes import Posterior, individualise
    from lithos.models import RENAL_2CMT, get_engine
    cov = cov_typical(age_at_t0=70, creatinine=[(0, 160)])     # older, poor kidney function
    admins = [Administration(21, 500)]
    doses = expand(21, 24 * 40, admins)
    truth = RENAL_2CMT.individual(cov(0), {})
    lv = [Level(24 * d + 9, float(simulate(doses, truth, [24 * d + 9])[0])) for d in (20, 35)]
    post = individualise(Case(cov, doses, lv), get_engine())
    assert isinstance(post, Posterior) and post.model.name == RENAL_2CMT.name


def test_ensemble_weights_are_a_distribution_for_eligible_adults():
    from lithos.bayes import Ensemble, individualise
    from lithos.models import get_engine
    cov = cov_typical(age_at_t0=35, creatinine=[(0, 70)])
    admins = [Administration(21, 750)]
    doses = expand(21, 24 * 40, admins)
    truth = MODEL.individual(cov(0), {})
    lv = [Level(24 * d + 9, float(simulate(doses, truth, [24 * d + 9])[0])) for d in (20, 35)]
    ens = individualise(Case(cov, doses, lv), get_engine())
    assert isinstance(ens, Ensemble) and len(ens.members) == 2
    assert ens.weights.sum() == pytest.approx(1.0)


def test_uninformative_sample_prompts_a_repeat_not_a_dose_change():
    from lithos.bayes import individualise
    from lithos.models import get_engine
    cov = cov_typical(age_at_t0=58, weight=72, height=163, creatinine=[(0, 82)])
    admins = [Administration(8.0, 500), Administration(20.0, 500)]
    doses = expand(8.0, 24 * 40 + 8.1, admins)
    t = 24 * 40 + 10.0                                          # 2 h after the morning dose
    post = individualise(Case(cov, doses, [Level(t, 1.24, timing_sd=0.25)]), get_engine())
    rec = recommend(post, LITHICARB, template="bd", current=admins)
    assert rec.repeat_first
    assert rec.chosen.daily_mg == 1000


def test_starting_dose_is_more_conservative_than_titration_rule():
    prior = fit(MODEL, Case(cov_typical(), [], []))
    start = initiation_plan(prior, LITHICARB)                   # P(> 1.0) <= 5%
    titr = recommend(prior, LITHICARB, t=0.0, high=1.2)          # P(> 1.2) <= 5%
    assert start.maintenance.daily_mg <= titr.chosen.daily_mg
    assert start.maintenance.p_high <= 0.05


def test_pregnancy_effects_and_in_pregnancy_update():
    from lithos.bayes import fit as fit_
    from lithos.case import Case as Case_
    cov = cov_typical(age_at_t0=31, weight=66, height=167, creatinine=[(0, 68)])
    admins = [Administration(21, 1000)]
    conception, delivery = 24 * 100.0, 24 * (100 + 280.0)
    eff = pregnancy.effects(conception, delivery)
    assert eff and all(e.cl_mult > 0 for e in eff)
    doses = expand(21, conception + 24 * 7 * 20, admins)
    truth = MODEL.individual(cov(0), {})
    t20 = conception + 24 * 7 * 20 - 12
    case = Case_(cov, doses, [Level(24 * 60 + 9, float(simulate(doses, truth, [24 * 60 + 9])[0])),
                              Level(t20, 0.55)], effects=eff)
    post = fit_(MODEL, case)
    steps = pregnancy.plan(post, LITHICARB, admins, conception_h=conception, delivery_h=delivery,
                           weeks=(20, 28, 36), postpartum_weeks=(0, 6))
    assert steps[0].daily_mg >= 1000                     # a low mid-pregnancy level supports a higher dose
    assert steps[-1].daily_mg <= steps[2].daily_mg       # after delivery the dose comes back down
    assert len(pregnancy.delivery_notes()) == 4


def test_interaction_lookup_handles_brands_combinations_and_unknowns():
    assert interactions.lookup("Nurofen").key == "nsaid"
    assert interactions.lookup("perindopril arginine").key == "acei"
    assert {i.key for i in interactions.lookup_all("Coversyl Plus")} == {"acei", "thiazide"}
    assert {i.key for i in interactions.lookup_all("irbesartan/hydrochlorothiazide")} == {"arb", "thiazide"}
    found, unknown = interactions.screen(["Nurofen", "paracetamol", "mysterymab"])
    assert [f.key for f in found] == ["nsaid"] and unknown == ["mysterymab"]
    assert interactions.CATALOG["acei"].followup_days == (7, 28)
