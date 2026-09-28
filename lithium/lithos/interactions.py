"""Drug interactions as quantitative, uncertain forecasts.

Most lithium toxicity is not caused by the lithium dose. It is caused by
something *else* changing: a thiazide started for blood pressure, an ACE
inhibitor, a course of NSAIDs, gastroenteritis, a heatwave. Existing tools
fire a generic "interaction: monitor levels" alert. This module instead
predicts, for *this* patient, what their 12-hour level will become, how sure
we are, and what dose would keep them in range. The alert can then fire in the
GP's prescribing software at the moment the interacting drug is chosen.

Pharmacokinetic effects are multiplicative changes in lithium clearance with an
uncertainty range (5th-95th percentile). Values are summarised from the
literature cited in ``evidence`` and should be refined as data accrue;
pharmacodynamic interactions (neurotoxicity, serotonin toxicity) carry advice
but no clearance change.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .bayes import Posterior
from .case import Administration
from .forecast import Summary, li12_per_unit, summarise


@dataclass(frozen=True)
class Interaction:
    key: str
    label: str
    examples: tuple[str, ...]
    cl_mult: float                    # median multiplicative change in lithium clearance
    cl_mult_range: tuple[float, float]  # ~5th-95th percentile of that change
    onset_days: float                 # time to (most of) the effect
    kind: str                         # "pk" or "pd"
    advice: str
    evidence: str

    def sample(self, n: int, rng: np.random.Generator) -> np.ndarray:
        """Log-normal draws with the stated median and 90% range."""
        if self.kind != "pk":
            return np.ones(n)
        lo, hi = self.cl_mult_range
        sd = (np.log(hi) - np.log(lo)) / (2 * 1.645)
        return np.exp(rng.normal(np.log(self.cl_mult), sd, size=n))

    @property
    def level_change(self) -> str:
        """Approximate steady-state level change implied by the clearance change."""
        pct = lambda m: 100.0 * (1.0 / m - 1.0)
        lo, hi = sorted((pct(self.cl_mult_range[1]), pct(self.cl_mult_range[0])))
        return f"{pct(self.cl_mult):+.0f}% (range {lo:+.0f}% to {hi:+.0f}%)"


CATALOG: dict[str, Interaction] = {i.key: i for i in [
    Interaction(
        "thiazide", "Thiazide / thiazide-like diuretic",
        ("hydrochlorothiazide", "chlorthalidone", "chlortalidone", "indapamide", "bendroflumethiazide",
         "chlorothiazide", "metolazone"),
        0.72, (0.55, 0.90), 5, "pk",
        "Avoid if possible. If essential: reduce lithium dose to keep the forecast in range, "
        "check a level 5-7 days after starting and after any diuretic dose change.",
        "Clearance falls ~25-40% (levels rise ~25-50%) via proximal sodium/lithium reabsorption; "
        "see docs/EVIDENCE.md."),
    Interaction(
        "acei", "ACE inhibitor",
        ("ramipril", "perindopril", "enalapril", "lisinopril", "captopril", "trandolapril",
         "quinapril", "fosinopril"),
        0.78, (0.55, 1.00), 14, "pk",
        "Prefer an alternative antihypertensive (e.g. a dihydropyridine calcium-channel blocker). "
        "If started: level within 1-2 weeks, then again at ~4 weeks; watch older patients closely.",
        "Levels rise variably (0 to >50%); onset can be delayed weeks; ~7-fold higher risk of "
        "hospitalisation for toxicity in older adults in the first month (population data); see docs/EVIDENCE.md."),
    Interaction(
        "arb", "Angiotensin-II receptor blocker",
        ("candesartan", "irbesartan", "losartan", "telmisartan", "valsartan", "olmesartan"),
        0.82, (0.60, 1.00), 14, "pk",
        "As for ACE inhibitors: level within 1-2 weeks of starting and after dose changes.",
        "Case reports and small series suggest rises similar to ACE inhibitors; weaker evidence."),
    Interaction(
        "nsaid", "NSAID (incl. COX-2 selective)",
        ("ibuprofen", "naproxen", "diclofenac", "indometacin", "indomethacin", "meloxicam",
         "celecoxib", "etoricoxib", "piroxicam", "ketorolac", "parecoxib"),
        0.80, (0.60, 0.98), 4, "pk",
        "Prefer paracetamol. If an NSAID is needed: shortest course, consider a dose reduction, "
        "level after 4-5 days; include over-the-counter ibuprofen in patient education.",
        "Levels rise ~10-60% depending on agent and person (indometacin, diclofenac at the high end); "
        "low-dose aspirin and sulindac have little effect."),
    Interaction(
        "loop", "Loop diuretic",
        ("furosemide", "frusemide", "bumetanide", "torasemide", "torsemide", "ethacrynic acid"),
        0.90, (0.70, 1.05), 5, "pk",
        "Smaller direct effect than thiazides, but volume depletion drives toxicity: level within a week, "
        "reinforce fluid advice, extra caution in older adults.",
        "Direct effect modest and variable; population data show elevated hospitalisation risk in older adults."),
    Interaction(
        "sglt2", "SGLT2 inhibitor",
        ("empagliflozin", "dapagliflozin", "canagliflozin", "ertugliflozin"),
        1.20, (1.00, 1.50), 7, "pk",
        "May LOWER lithium levels (natriuresis): check a level 1-2 weeks after starting; watch for relapse.",
        "Case reports and small studies of falling levels; direction consistent, magnitude uncertain."),
    Interaction(
        "xanthine", "Theophylline / aminophylline / high caffeine",
        ("theophylline", "aminophylline", "caffeine"),
        1.25, (1.05, 1.50), 3, "pk",
        "Lowers levels; stopping them raises levels. Check a level after starting or stopping.",
        "Increased renal lithium clearance with methylxanthines."),
    Interaction(
        "acetazolamide", "Carbonic anhydrase inhibitor / osmotic diuretic",
        ("acetazolamide", "mannitol", "urea"),
        1.25, (1.00, 1.60), 2, "pk",
        "Lowers levels; used therapeutically in toxicity only under specialist care.",
        "Reduced proximal reabsorption increases lithium clearance."),
    Interaction(
        "metronidazole", "Metronidazole",
        ("metronidazole", "tinidazole"),
        0.85, (0.60, 1.00), 5, "pk",
        "Case reports of toxicity: consider a level during longer courses.",
        "Case reports only."),
    Interaction(
        "carbamazepine", "Carbamazepine (neurotoxicity)",
        ("carbamazepine",),
        1.0, (1.0, 1.0), 0, "pd",
        "Neurotoxicity reported at therapeutic lithium levels; monitor for neurological symptoms.",
        "Pharmacodynamic; case series."),
    Interaction(
        "ccb_nondhp", "Verapamil / diltiazem (neurotoxicity)",
        ("verapamil", "diltiazem"),
        1.0, (1.0, 1.0), 0, "pd",
        "Neurotoxicity and bradycardia reported; dihydropyridines (e.g. amlodipine) are preferred.",
        "Pharmacodynamic; case reports."),
    Interaction(
        "serotonergic", "Serotonergic drugs (serotonin toxicity)",
        ("sertraline", "fluoxetine", "paroxetine", "citalopram", "escitalopram", "venlafaxine",
         "duloxetine", "tramadol", "sumatriptan", "linezolid", "methylene blue"),
        1.0, (1.0, 1.0), 0, "pd",
        "Commonly combined safely; educate about serotonin toxicity symptoms, especially with dose escalation.",
        "Pharmacodynamic; case reports."),
    Interaction(
        "antipsychotic_high", "High-dose antipsychotics (neurotoxicity)",
        ("haloperidol",),
        1.0, (1.0, 1.0), 0, "pd",
        "Rare neurotoxicity/NMS-like reactions reported, mainly with high doses; monitor.",
        "Pharmacodynamic; case reports."),
]}

_INDEX = {name: key for key, it in CATALOG.items() for name in it.examples}


def lookup(drug: str) -> Interaction | None:
    return CATALOG.get(_INDEX.get(drug.strip().lower(), ""))


def screen(drugs: Sequence[str]) -> list[Interaction]:
    found = [lookup(d) for d in drugs]
    return [f for f in found if f is not None]


@dataclass
class WhatIf:
    interaction: Interaction
    before: Summary
    after: Summary
    p_above_1: float
    p_above_1_2: float
    adjusted: "Recommendation | None" = None   # dose that keeps the target with the new drug


def what_if(post: Posterior, admins: Sequence[Administration], drug: str, t: float | None = None,
            n: int = 3000, rng: np.random.Generator | None = None, product=None,
            target: tuple[float, float] = (0.6, 0.8)) -> WhatIf | None:
    """Forecast this patient's 12-hour level if ``drug`` is added with no lithium
    change, and (given a ``product``) the lithium dose that would compensate."""
    it = lookup(drug)
    if it is None:
        return None
    rng = rng or np.random.default_rng(11)
    if t is None:
        t = max([lv.time for lv in post.case.levels], default=0.0)
    daily = sum(a.mmol for a in admins)
    seed = int(rng.integers(1 << 31))
    mult = it.sample(n, rng)
    # Identical posterior draws before and after, so the difference is the interaction alone.
    before = daily * li12_per_unit(post, admins, t, n, np.random.default_rng(seed))
    after = daily * li12_per_unit(post, admins, t, n, np.random.default_rng(seed), cl_mult=mult)
    adjusted = None
    if product is not None and it.kind == "pk":
        from .dosing import recommend  # local import: dosing depends on forecast, not on interactions
        template = [a.clock for a in admins]
        adjusted = recommend(post, product, target=target, template=template, current=admins, t=t,
                             cl_mult=mult, n=n, rng=np.random.default_rng(seed), max_step_ratio=2.0)
    return WhatIf(it, summarise(before), summarise(after),
                  float(np.mean(after > 1.0)), float(np.mean(after > 1.2)), adjusted)
