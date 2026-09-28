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


def recommend(
    post: Posterior,
    product: Product,
    *,
    target: tuple[float, float] = (0.6, 0.8),
    template: str | Sequence[float] = "nocte",
    current: Sequence[Administration] | None = None,
    t: float | None = None,
    high: float = 1.0,
    p_high_max: float = 0.05,
    max_daily_mg: float = 2000.0,
    max_step_ratio: float = 1.5,
    allow_half: bool = False,
    cl_mult: float | np.ndarray = 1.0,
    n: int = 3000,
    rng: np.random.Generator | None = None,
) -> Recommendation:
    """Recommend a regimen for ``target`` 12-hour levels (mmol/L).

    ``high`` is the level whose probability is capped at ``p_high_max``; for
    maintenance a Li12 above 1.0 mmol/L is the usual line. ``max_step_ratio``
    limits how far one adjustment can move the daily dose (larger jumps are
    offered but not auto-selected). ``cl_mult`` lets callers ask "what if"
    (e.g. an interacting drug's uncertain effect, one value per draw).
    """
    rng = rng or np.random.default_rng(7)
    if t is None:
        t = max([lv.time for lv in post.case.levels], default=max([d.time for d in post.case.doses], default=0.0))
    clocks = list(TEMPLATES[template]) if isinstance(template, str) else list(template)
    G, _ = unit_contributions(post, clocks, product.release, t, n, rng, cl_mult)
    lo, hi = target

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
            Gc, _ = unit_contributions(post, cur_clocks, forms.pop(), t, n, rng, cl_mult)
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
        within = [o for o in safe if o.daily_mg <= d0 * max_step_ratio and o.daily_mg >= d0 / max_step_ratio]
        if within:
            safe_step = within
        else:
            safe_step = safe
            warnings.append("the best option changes the dose by more than the step limit; consider a staged change")
    else:
        safe_step = safe
    if not safe_step:
        chosen = min(options, key=lambda o: o.p_high)
        warnings.append(f"no regimen keeps P(Li12 > {high}) below {p_high_max:.0%}; lowest-risk option shown")
    else:
        chosen = sorted(safe_step, key=rank_key)[0]
        # Keep the current regimen if it is essentially as good (avoid churn).
        if current_opt is not None and current_opt.p_high <= p_high_max and \
                current_opt.p_target >= chosen.p_target - 0.05 and lo <= current_opt.li12.median <= hi:
            chosen = current_opt

    # Know when not to answer: if the current regimen is plausibly fine, nothing
    # is dangerous, and the evidence is simply too uncertain, a better-timed
    # level beats a dose change driven by uncertainty alone.
    repeat_first = False
    if current_opt is not None and chosen is not current_opt:
        c = current_opt.li12
        vague = (c.hi - c.lo) > 0.35
        plausible = lo - 0.1 <= c.median <= hi + 0.1
        safe_now = float(np.mean(c.samples > 1.2)) < 0.10
        if vague and plausible and safe_now:
            repeat_first = True
            chosen = current_opt
            warnings.append(f"evidence too uncertain to justify a change (90% CrI {c.lo:.2f}-{c.hi:.2f}); "
                            "repeat a level at least 10 h after the evening dose, before any morning dose")

    hl = half_life_h(post, t, n=400, rng=rng)
    t_half_days = hl.median / 24.0
    recheck_conv = max(5, ceil(5 * t_half_days))
    recheck = max(3, ceil(2.5 * t_half_days))

    rationale = []
    if current_opt is not None:
        rationale.append(f"Current regimen ({current_opt.describe(product)}): standardised 12-h level "
                         f"{current_opt.li12}; P(in {lo}-{hi}) = {current_opt.p_target:.0%}.")
    rationale.append(f"Recommended ({chosen.describe(product)}): predicted 12-h level {chosen.li12}; "
                     f"P(in {lo}-{hi}) = {chosen.p_target:.0%}; P(> {high}) = {chosen.p_high:.1%}.")
    rationale.append(f"Estimated half-life {hl.median:.0f} h (90% CrI {hl.lo:.0f}-{hl.hi:.0f}): "
                     f"~90% of a new steady state is reached after {3.3 * t_half_days:.1f} days.")
    ranked = sorted(options, key=rank_key)[:5]
    if repeat_first:
        recheck = 1
    return Recommendation(target, chosen, ranked, current_opt, product, recheck, recheck_conv, hl, rationale,
                          warnings, repeat_first)


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
    kwargs.setdefault("high", 1.0)      # a starting dose is chosen conservatively: P(level > 1.0) <= 5%
    rec = recommend(prior, product, target=target, template=template, t=0.0, **kwargs)
    step = product.strength_mg
    lead: list[Administration] = []
    for a in rec.chosen.admins:
        mg = max(step, np.floor(a.mg * lead_in_fraction / step) * step)
        lead.append(Administration(a.clock, float(mg), a.formulation, a.salt))
    first_level = lead_in_days + rec.recheck_days
    return InitiationPlan(lead, lead_in_days, rec.chosen, first_level, rec)
