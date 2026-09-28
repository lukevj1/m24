"""Dose recommendation: the smallest decision a prescriber has to make.

For each candidate regimen built from real tablets, the engine computes the
posterior distribution of the steady-state 12-hour level and picks the regimen
most likely to land in the target window while keeping the probability of an
excessive level low. Because kinetics are linear, Li12 = G @ dose where G holds
per-sample, per-administration unit contributions, so thousands of candidates
cost one matrix product.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import ceil
from typing import Sequence

import numpy as np

from .bayes import Posterior
from .case import Administration, li12_clock
from .forecast import Summary, half_life_h, summarise
from .pk import steady_state
from .units import Product, mg_to_mmol

def default_target(age: float, *, postpartum_first_month: bool = False) -> tuple[float, float]:
    """Starting-point targets for the standardised 12-hour level (mmol/L).

    Adults: 0.60-0.80 (ISBD/IGSLi 2019; NICE CG185 for people starting
    lithium); individualise to 0.40-0.60 for good response with poor
    tolerability, or 0.80-1.00 for poor response with good tolerability.
    Age 65+: 0.40-0.60 (ISBD/IGSLi majority view; NICE NG222), with maxima of
    0.70-0.80 at 65-79 and 0.70 over 80. First postpartum month: 0.80-1.00
    for relapse prevention (Poels 2018).
    """
    if postpartum_first_month:
        return (0.8, 1.0)
    if age >= 65:
        return (0.4, 0.6)
    return (0.6, 0.8)


TEMPLATES = {
    "nocte": [21.0],          # once daily at night: the usual choice and kinder to kidneys
    "bd": [8.0, 20.0],
}


@dataclass
class DoseOption:
    admins: list[Administration]
    li12: Summary
    p_target: float
    p_high: float
    p_low: float

    @property
    def daily_mg(self) -> float:
        return sum(a.mg for a in self.admins)

    def describe(self, product: Product | None = None) -> str:
        parts = []
        for a in self.admins:
            if product is not None:
                n = a.mg / product.strength_mg
                count = f"{n:g} x {product.strength_mg:g} mg"
                parts.append(f"{a.mg:g} mg ({count}) at {int(a.clock):02d}:{int(round(a.clock % 1 * 60)):02d}")
            else:
                parts.append(f"{a.mg:g} mg at {int(a.clock):02d}:00")
        return " + ".join(parts)


@dataclass
class Recommendation:
    target: tuple[float, float]
    chosen: DoseOption
    options: list[DoseOption]
    current: DoseOption | None
    product: Product
    recheck_days: int
    recheck_days_conventional: int
    half_life: Summary
    rationale: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    repeat_first: bool = False   # evidence too weak to justify a change: repeat a well-timed level first
    check_adherence: bool = False
    high: float = 1.0

    @property
    def change(self) -> str:
        if self.repeat_first:
            return "repeat level before changing"
        if self.current is None:
            return "start"
        d0, d1 = self.current.daily_mg, self.chosen.daily_mg
        return "no change" if abs(d1 - d0) < 1e-6 else ("increase" if d1 > d0 else "decrease")


def unit_contributions(post: Posterior, clocks: Sequence[float], formulation: str, t: float,
                       n: int = 3000, rng: np.random.Generator | None = None,
                       cl_mult: float | np.ndarray = 1.0) -> tuple[np.ndarray, float]:
    """G[i, a]: steady-state Li12 (mmol/L) per 1 mmol/day given at ``clocks[a]``
    for posterior draw i, with Li12 taken 12 h after the last clock time."""
    li_clock = (max(clocks) + 12.0) % 24.0
    params = post.param_draws(n, t, rng, cl_mult=cl_mult)
    G = np.array([[steady_state([(c, 1.0, formulation)], p, li_clock)[0] for c in clocks] for p in params])
    return G, li_clock


def _candidates(product: Product, clocks: Sequence[float], max_daily_mg: float, allow_half: bool):
    step = product.strength_mg / 2 if (allow_half and product.scored) else product.strength_mg
    k_max = int(max_daily_mg // step)
    for k in range(1, k_max + 1):
        if len(clocks) == 1:
            yield [k * step]
        else:
            # Split as evenly as tablets allow; any remainder goes to the evening dose.
            base, rem = divmod(k, len(clocks))
            if base == 0:
                continue
            split = [base * step] * len(clocks)
            for r in range(rem):
                split[-1 - r] += step
            yield split


# Model-independent guardrails (local clinical governance should set these).
# They hold whatever the model says, because every variance term in the model
# is itself an estimate.
GUARDRAILS = {
    "max_start_mg": 1000.0,             # starting maintenance dose, adults < 65 with CrCl >= 60
    "max_start_mg_older_or_ckd": 750.0,  # age >= 65 or CrCl 30-59
    "no_start_below_crcl": 30.0,        # US labelling: not recommended below CrCl 30 mL/min
    "max_step_mg": 500.0,               # largest single change in daily dose
    "max_step_ratio": 1.5,              # largest relative change in daily dose
    "adherence_check_ratio": 1.4,       # clearance this far above kidney-predicted -> confirm adherence
}


def recommend(
    post: Posterior,
    product: Product,
    *,
    target: tuple[float, float] = (0.6, 0.8),
    template: str | Sequence[float] = "nocte",
    current: Sequence[Administration] | None = None,
    t: float | None = None,
    high: float | None = None,
    p_high_max: float = 0.05,
    max_daily_mg: float = 2000.0,
    max_step_ratio: float | None = None,
    allow_half: bool = False,
    cl_mult: float | np.ndarray = 1.0,
    horizon_days: float = 28.0,
    n: int = 3000,
    rng: np.random.Generator | None = None,
) -> Recommendation:
    """Recommend a regimen for the clinician's ``target`` 12-hour level (mmol/L).

    * The safety ceiling follows the target: by default the probability of a
      12-hour level above ``target[1] + 0.2`` must stay under ``p_high_max``
      (so 1.0 for a 0.6-0.8 target, 1.2 for a 0.8-1.0 target). If the ceiling
      stops the forecast reaching the target, the card says so.
    * Doses are chosen for the clearance expected over the next
      ``horizon_days`` (short-term drift relaxes), not a transient dip.
    * Steps are limited (ratio and absolute mg); larger moves are offered as
      alternatives but not auto-selected.
    * If the patient appears to clear lithium much faster than their kidneys
      predict, an increase comes with "confirm adherence first".
    """
    rng = rng or np.random.default_rng(7)
    if t is None:
        t = max([lv.time for lv in post.case.levels], default=max([d.time for d in post.case.doses], default=0.0))
    lo, hi = target
    high = hi + 0.2 if high is None else high
    max_step_ratio = GUARDRAILS["max_step_ratio"] if max_step_ratio is None else max_step_ratio
    clocks = list(TEMPLATES[template]) if isinstance(template, str) else list(template)
    t_draw = t + horizon_days * 24.0
    G, _ = unit_contributions(post, clocks, product.release, t_draw, n, rng, cl_mult)

    def evaluate(mgs: Sequence[float], clks: Sequence[float], Gm: np.ndarray) -> DoseOption:
        mmol = np.array([mg_to_mmol(m, product.salt) for m in mgs])
        s = Gm @ mmol
        return DoseOption([Administration(c, m, product.release, product.salt) for c, m in zip(clks, mgs)],
                          summarise(s), float(np.mean((s >= lo) & (s <= hi))),
                          float(np.mean(s > high)), float(np.mean(s < lo)))

    options = [evaluate(mgs, clocks, G) for mgs in _candidates(product, clocks, max_daily_mg, allow_half)]
    current_opt = None
    if current:
        cur_clocks = [a.clock for a in current]
        forms = {a.formulation for a in current}
        if cur_clocks == clocks and forms == {product.release}:
            Gc = G
        else:
            Gc, _ = unit_contributions(post, cur_clocks, forms.pop(), t_draw, n, rng, cl_mult)
        mm = np.array([a.mmol for a in current])
        s = Gc @ mm
        current_opt = DoseOption(list(current), summarise(s), float(np.mean((s >= lo) & (s <= hi))),
                                 float(np.mean(s > high)), float(np.mean(s < lo)))

    center = 0.5 * (lo + hi)

    def rank_key(o: DoseOption):
        return (-round(o.p_target, 2), abs(o.li12.median - center))

    safe = [o for o in options if o.p_high <= p_high_max]
    warnings: list[str] = []
    if current_opt is not None:
        d0 = current_opt.daily_mg
        # One tablet up or down is always allowed: with 250 mg tablets a ratio
        # limit alone would lock a patient on 250 mg (or 500 mg) for good.
        one = product.strength_mg / 2 if (allow_half and product.scored) else product.strength_mg
        up = max(d0 * max_step_ratio, d0 + one)
        down = min(d0 / max_step_ratio, d0 - one)
        within = [o for o in safe if down - 1e-6 <= o.daily_mg <= up + 1e-6
                  and abs(o.daily_mg - d0) <= max(GUARDRAILS["max_step_mg"], one) + 1e-6]
        if within:
            safe_step = within
        else:
            safe_step = safe
            warnings.append("the best option changes the dose by more than the step limit; consider a staged change")
    else:
        safe_step = safe
    if not safe_step:
        chosen = min(options, key=lambda o: o.p_high)
        warnings.append(f"no regimen keeps P(12-h level > {high:.1f}) below {p_high_max:.0%}; lowest-risk option shown")
    else:
        chosen = sorted(safe_step, key=rank_key)[0]
        # Keep the current regimen if it is in range and nearly as good: with
        # tablet-sized steps two doses often straddle a narrow window, and
        # chasing a few points of probability means flip-flopping between them.
        if current_opt is not None and current_opt.p_high <= p_high_max and \
                current_opt.p_target >= chosen.p_target - 0.10 and lo <= current_opt.li12.median <= hi:
            chosen = current_opt
    if chosen.li12.median < lo and any(lo <= o.li12.median <= hi and o.p_high > p_high_max for o in options):
        warnings.append(f"the safety ceiling (P(> {high:.1f}) <= {p_high_max:.0%}) stops the forecast reaching the "
                        "target with this much uncertainty; another level will narrow it before going higher")

    # Know when not to answer: if the current regimen is plausibly fine, nothing
    # is dangerous, and the evidence is simply too uncertain, a better-timed
    # level beats a dose change driven by uncertainty alone. Not for a
    # hypothetical (an interaction or pregnancy forecast): there the width comes
    # from the effect being forecast, and another level today cannot narrow it.
    hypothetical = bool(np.any(np.asarray(cl_mult, dtype=float) != 1.0))
    repeat_first = False
    if current_opt is not None and chosen is not current_opt and not hypothetical:
        c = current_opt.li12
        vague = (c.hi - c.lo) > 0.35
        plausible = lo - 0.1 <= c.median <= hi + 0.1
        safe_now = float(np.mean(c.samples > 1.2)) < 0.10
        if vague and plausible and safe_now:
            repeat_first = True
            chosen = current_opt
            warnings.append(f"evidence too uncertain to justify a change (90% CrI {c.lo:.2f}-{c.hi:.2f}); "
                            "repeat a level at least 10 h after the evening dose, before any morning dose")

    # Low levels with apparently fast clearance are as often missed doses as fast kidneys.
    check_adherence = False
    if current_opt is not None and chosen.daily_mg > current_opt.daily_mg and post.case.levels:
        ratio = post.clearance_ratio(t)
        if ratio > GUARDRAILS["adherence_check_ratio"]:
            check_adherence = True
            warnings.append(f"levels are lower than kidney function predicts (clearance about {ratio:.1f}x expected): "
                            "confirm adherence and dose timing with the patient before increasing")

    hl = half_life_h(post, t, n=400, rng=rng)
    t_half_days = hl.median / 24.0
    recheck_conv = max(5, ceil(5 * t_half_days))
    recheck = max(3, ceil(2.5 * t_half_days))

    rationale = []
    if current_opt is not None:
        rationale.append(f"Current regimen ({current_opt.describe(product)}): standardised 12-h level "
                         f"{current_opt.li12}; P(in {lo}-{hi}) = {current_opt.p_target:.0%}.")
    rationale.append(f"Recommended ({chosen.describe(product)}): predicted 12-h level {chosen.li12}; "
                     f"P(in {lo}-{hi}) = {chosen.p_target:.0%}; P(below {lo}) = {chosen.p_low:.0%}; "
                     f"P(> {high:.1f}) = {chosen.p_high:.1%}.")
    rationale.append(f"Estimated half-life {hl.median:.0f} h (90% CrI {hl.lo:.0f}-{hl.hi:.0f}): "
                     f"~90% of a new steady state is reached after {3.3 * t_half_days:.1f} days.")
    ranked = sorted(options, key=rank_key)[:5]
    if repeat_first:
        recheck = 1
    rec = Recommendation(target, chosen, ranked, current_opt, product, recheck, recheck_conv, hl, rationale,
                         warnings, repeat_first)
    rec.check_adherence = check_adherence
    rec.high = high
    return rec


@dataclass
class InitiationPlan:
    lead_in: list[Administration]
    lead_in_days: int
    maintenance: DoseOption
    first_level_day: int
    recommendation: Recommendation

    def describe(self, product: Product) -> str:
        lead = " + ".join(f"{a.mg:g} mg at {int(a.clock):02d}:00" for a in self.lead_in)
        return (f"Days 1-{self.lead_in_days}: {lead}; then {self.maintenance.describe(product)} "
                f"(forecast 12-h level {self.maintenance.li12}); first level on day {self.first_level_day} "
                f"at any convenient time >= 6 h after a dose.")


def initiation_plan(prior: Posterior, product: Product, *, target: tuple[float, float] = (0.6, 0.8),
                    template: str | Sequence[float] = "nocte", lead_in_fraction: float = 0.5,
                    lead_in_days: int = 3, **kwargs) -> InitiationPlan:
    """Start low for tolerability, then step straight to the model-predicted
    maintenance dose instead of creeping up week by week.

    The first level is scheduled a few days after reaching the maintenance dose;
    the Bayesian update does not need steady state, so it can come early.
    """
    cov = prior.case.covariates(0.0)
    if cov.crcl < GUARDRAILS["no_start_below_crcl"]:
        raise ValueError(f"CrCl {cov.crcl:.0f} mL/min: lithium is not recommended below "
                         f"{GUARDRAILS['no_start_below_crcl']:.0f} mL/min (US labelling); specialist decision")
    cap = GUARDRAILS["max_start_mg_older_or_ckd"] if (cov.age >= 65 or cov.crcl < 60) else GUARDRAILS["max_start_mg"]
    kwargs.setdefault("max_daily_mg", cap)
    rec = recommend(prior, product, target=target, template=template, t=0.0, horizon_days=0.0, **kwargs)
    step = product.strength_mg
    lead: list[Administration] = []
    for a in rec.chosen.admins:
        mg = max(step, np.floor(a.mg * lead_in_fraction / step) * step)
        lead.append(Administration(a.clock, float(mg), a.formulation, a.salt))
    first_level = lead_in_days + rec.recheck_days
    return InitiationPlan(lead, lead_in_days, rec.chosen, first_level, rec)
