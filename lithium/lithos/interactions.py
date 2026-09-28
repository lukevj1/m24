"""Drug interactions as quantitative, uncertain forecasts.

Most lithium toxicity is not caused by the lithium dose. It is caused by
something *else* changing: a thiazide started for blood pressure, an ACE
inhibitor, a course of NSAIDs, gastroenteritis, a heatwave. Existing tools
fire a generic "interaction: monitor levels" alert. This module instead
predicts, for *this* patient, what their 12-hour level will become, how sure
we are, and what dose would keep them in range. The alert can then fire in the
GP's prescribing software at the moment the interacting drug is chosen.

Pharmacokinetic effects are multiplicative changes in lithium clearance with an
uncertainty range (about the 5th-95th percentile). Values are summarised from
the sources in ``evidence`` (docs/EVIDENCE.md) and should be refined as data
accrue; pharmacodynamic interactions (neurotoxicity, serotonin toxicity) carry
advice but no clearance change.

Two kinds of evidence point in different directions and both are used. The
pharmacokinetic literature gives thiazides the largest mean level rise, while
the epidemiology in older adults gives ACE inhibitors and loop diuretics the
largest toxicity-admission risk in the first month (thiazides and NSAIDs were
not independently associated). So thiazides call for dose arithmetic, and ACE
inhibitors and loop diuretics call for a time-anchored first-month monitoring
plan whatever the forecast says.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import re
from dataclasses import replace

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
    followup_days: tuple[int, ...] = ()   # check levels after starting (days)

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
        0.76, (0.60, 0.92), 5, "pk",
        "Avoid if possible. If essential: reduce the lithium dose to keep the forecast in range and check a "
        "level 5-7 days after starting and after any diuretic dose change.",
        "Levels commonly rise 25-40% (Finley 2016 review, secondary); renal lithium clearance fell 24% on "
        "thiazide therapy; HCTZ raised serum lithium more than furosemide or placebo in a crossover study "
        "(Crabtree 1991). Not independently associated with toxicity admission in older adults (Juurlink 2004), "
        "plausibly because the interaction is well known and managed."),
    Interaction(
        "acei", "ACE inhibitor",
        ("ramipril", "perindopril", "enalapril", "lisinopril", "captopril", "trandolapril",
         "quinapril", "fosinopril"),
        0.75, (0.55, 1.00), 14, "pk",
        "Prefer another antihypertensive where possible. If started: level within 1-2 weeks and again at "
        "about 4 weeks; the first month is the danger window, especially in older adults.",
        "Steady-state levels rose 36.1% and clearance fell 25.5% after starting an ACE inhibitor (Finley 1996, "
        "n=20; 4 with toxicity symptoms). Toxicity admission RR 7.6 (2.6-22.0) within one month of starting in "
        "people aged 66+ (Juurlink 2004)."),
    Interaction(
        "arb", "Angiotensin-II receptor blocker",
        ("candesartan", "irbesartan", "losartan", "telmisartan", "valsartan", "olmesartan"),
        0.80, (0.60, 1.00), 14, "pk",
        "As for ACE inhibitors: level within 1-2 weeks of starting and after dose changes.",
        "Case reports with losartan, valsartan and candesartan; the candesartan label reports increased lithium "
        "levels. No controlled effect size found."),
    Interaction(
        "nsaid", "NSAID (incl. COX-2 selective)",
        ("ibuprofen", "naproxen", "diclofenac", "indometacin", "indomethacin", "meloxicam",
         "celecoxib", "etoricoxib", "piroxicam", "ketorolac", "parecoxib", "etodolac", "nabumetone"),
        0.85, (0.62, 1.00), 4, "pk",
        "Prefer paracetamol. If an NSAID is needed, prefer regular to as-needed use so levels are predictable, "
        "and check a level after 4-5 days (NICE CG185: monthly until stable, then 3-monthly). Include "
        "over-the-counter ibuprofen in patient education.",
        "US NSAID class labelling: mean minimum lithium concentration +15%, renal clearance about -20%; celecoxib "
        "+17%; meloxicam about +20%; indometacin and piroxicam 'significantly' more; case reports up to ~100%. "
        "Low-dose aspirin and sulindac have little effect."),
    Interaction(
        "loop", "Loop diuretic",
        ("furosemide", "frusemide", "bumetanide", "torasemide", "torsemide", "ethacrynic acid"),
        0.92, (0.75, 1.05), 5, "pk",
        "The direct effect is modest, but volume depletion drives toxicity. Level within a week, reinforce "
        "fluid advice and watch older adults closely through the first month.",
        "Less potent than HCTZ at raising levels (Crabtree 1991); the hazard is secondary volume depletion. "
        "Toxicity admission RR 5.5 (1.9-16.1) within one month of starting in people aged 66+ (Juurlink 2004)."),
    Interaction(
        "sglt2", "SGLT2 inhibitor",
        ("empagliflozin", "dapagliflozin", "canagliflozin", "ertugliflozin"),
        1.15, (1.00, 1.40), 7, "pk",
        "May LOWER lithium levels: check a level 1-2 weeks after starting or changing the dose, and watch "
        "for relapse.",
        "US lithium label (2026): may decrease serum lithium concentrations; NZ Medsafe 2023 and Malaysia NPRA "
        "2024 alerts. Magnitude uncertain."),
    Interaction(
        "xanthine", "Theophylline / aminophylline / caffeine",
        ("theophylline", "aminophylline", "caffeine"),
        1.25, (1.05, 1.50), 3, "pk",
        "Lowers levels; stopping them raises levels (caffeine withdrawal about +24%). Check a level after "
        "starting or stopping.",
        "Theophylline increased lithium clearance about 30% (Perry 1984, secondary); the US label lists "
        "xanthines as lowering levels."),
    Interaction(
        "acetazolamide", "Carbonic anhydrase inhibitor / osmotic diuretic / alkalinising agent",
        ("acetazolamide", "mannitol", "urea", "sodium bicarbonate"),
        1.25, (1.00, 1.60), 2, "pk",
        "Lowers levels by increasing urinary lithium excretion; more frequent monitoring.",
        "US lithium label, Table 4 (no effect size)."),
    Interaction(
        "metronidazole", "Metronidazole / nitroimidazoles",
        ("metronidazole", "tinidazole"),
        0.85, (0.60, 1.00), 5, "pk",
        "Reduced renal clearance reported: consider a level during longer courses.",
        "US lithium label: may increase levels via reduced renal clearance (no effect size)."),
    Interaction(
        "carbamazepine", "Carbamazepine / phenytoin / methyldopa",
        ("carbamazepine", "phenytoin", "methyldopa"),
        1.0, (1.0, 1.0), 0, "pd",
        "Increased risk of adverse reactions (including neurotoxicity) at therapeutic lithium levels; monitor "
        "for neurological symptoms.",
        "US lithium label (pharmacodynamic)."),
    Interaction(
        "ccb", "Calcium-channel blockers (neurotoxicity)",
        ("verapamil", "diltiazem", "amlodipine", "nifedipine", "felodipine", "lercanidipine"),
        1.0, (1.0, 1.0), 0, "pd",
        "Neurological adverse reactions reported (ataxia, tremor, nausea, tinnitus), mainly with verapamil and "
        "diltiazem; not a level-raising interaction.",
        "US lithium label; case reports (verapamil, diltiazem)."),
    Interaction(
        "serotonergic", "Serotonergic drugs (serotonin toxicity)",
        ("sertraline", "fluoxetine", "paroxetine", "citalopram", "escitalopram", "venlafaxine",
         "duloxetine", "tramadol", "sumatriptan", "linezolid", "methylene blue"),
        1.0, (1.0, 1.0), 0, "pd",
        "Commonly combined safely; educate about serotonin toxicity, especially at lithium initiation. "
        "Fluoxetine has been reported to raise or lower levels, so check a level after starting it.",
        "US lithium label."),
    Interaction(
        "antipsychotic", "Antipsychotics (neurotoxicity)",
        ("haloperidol", "olanzapine", "quetiapine", "risperidone", "aripiprazole", "chlorpromazine"),
        1.0, (1.0, 1.0), 0, "pd",
        "Usually combined safely; rare neurotoxicity, encephalopathic syndrome or NMS reported, mainly with "
        "high doses. Monitor neurological status.",
        "US lithium label (pharmacodynamic)."),
    Interaction(
        "mra", "Mineralocorticoid-receptor antagonist",
        ("spironolactone", "eplerenone"),
        0.90, (0.70, 1.05), 7, "pk",
        "Case reports of higher levels, and hyperkalaemia with RAAS agents; check a level after starting.",
        "Case reports only; magnitude uncertain."),
    Interaction(
        "cni", "Calcineurin inhibitor (nephrotoxic)",
        ("ciclosporin", "cyclosporine", "tacrolimus"),
        1.0, (1.0, 1.0), 0, "pd",
        "Nephrotoxic: tighter renal and level monitoring; involve the prescribing specialist.",
        "Additive nephrotoxicity; no controlled lithium PK data."),
    Interaction(
        "iodide", "Iodide preparations",
        ("potassium iodide",),
        1.0, (1.0, 1.0), 0, "pd",
        "Extended combined use may cause hypothyroidism; check TSH.",
        "US lithium label."),
]}

# Check levels after starting (days). ACE inhibitors, ARBs and loop diuretics get a
# day-28 level too: their toxicity risk is concentrated in the first month.
_FOLLOWUP = {"thiazide": (5, 14), "acei": (7, 28), "arb": (7, 28), "nsaid": (5,), "loop": (7, 28),
             "sglt2": (10,), "xanthine": (5,), "acetazolamide": (5,), "metronidazole": (5,), "mra": (7,)}
CATALOG = {k: replace(v, followup_days=_FOLLOWUP.get(k, ())) for k, v in CATALOG.items()}

_INDEX = {name: key for key, it in CATALOG.items() for name in it.examples}

# Common Australian/UK/US brand names -> ingredients (illustrative; production use needs AMT / SNOMED CT-AU
# or RxNorm ingredient mapping rather than a hand list).
BRANDS = {
    "nurofen": "ibuprofen", "brufen": "ibuprofen", "advil": "ibuprofen", "voltaren": "diclofenac",
    "naprosyn": "naproxen", "naprogesic": "naproxen", "celebrex": "celecoxib", "mobic": "meloxicam",
    "arcoxia": "etoricoxib", "indocid": "indomethacin", "toradol": "ketorolac",
    "coversyl": "perindopril", "coversyl plus": "perindopril/indapamide", "tritace": "ramipril",
    "renitec": "enalapril", "zestril": "lisinopril", "accupril": "quinapril",
    "avapro": "irbesartan", "karvea": "irbesartan", "karvezide": "irbesartan/hydrochlorothiazide",
    "avapro hct": "irbesartan/hydrochlorothiazide", "atacand": "candesartan",
    "atacand plus": "candesartan/hydrochlorothiazide", "micardis": "telmisartan",
    "micardis plus": "telmisartan/hydrochlorothiazide", "cozaar": "losartan",
    "hyzaar": "losartan/hydrochlorothiazide", "diovan": "valsartan", "co-diovan": "valsartan/hydrochlorothiazide",
    "olmetec": "olmesartan", "natrilix": "indapamide", "hygroton": "chlortalidone",
    "dithiazide": "hydrochlorothiazide", "lasix": "furosemide", "urex": "furosemide", "burinex": "bumetanide",
    "aldactone": "spironolactone", "inspra": "eplerenone", "jardiance": "empagliflozin",
    "forxiga": "dapagliflozin", "invokana": "canagliflozin", "flagyl": "metronidazole", "diamox": "acetazolamide",
    "tegretol": "carbamazepine", "isoptin": "verapamil", "cardizem": "diltiazem", "norvasc": "amlodipine",
}
_SALTS = {"arginine", "erbumine", "sodium", "potassium", "hydrochloride", "hcl", "maleate", "besylate",
          "mesylate", "cilexetil", "medoxomil", "tromethamine", "sr", "cr", "xr", "mr", "er", "tablets", "tablet",
          "capsules", "capsule", "mg"}
# Checked and not expected to change lithium levels (label or study data), so not reported as unknown.
KNOWN_NON_INTERACTING = {"paracetamol", "acetaminophen", "nefazodone", "gabapentin", "venlafaxine", "paroxetine",
                         "aripiprazole", "ziprasidone", "levothyroxine", "atorvastatin", "metformin"}


def _ingredients(drug: str) -> list[str]:
    text = drug.strip().lower()
    text = BRANDS.get(text, text)
    out = []
    for part in re.split(r"\s*(?:/|\+|,|&|\band\b|\bwith\b)\s*", text):
        part = BRANDS.get(part.strip(), part.strip())
        for sub in re.split(r"\s*/\s*", part):
            words = [w for w in re.findall(r"[a-z-]+", sub) if w not in _SALTS]
            if words:
                out.append(" ".join(words))
    return out


def lookup_all(drug: str) -> list[Interaction]:
    """Every interaction implied by a drug name, including combination products."""
    found: list[Interaction] = []
    for ing in _ingredients(drug):
        key = _INDEX.get(ing) or next((k for name, k in _INDEX.items() if name in ing.split()), None)
        if key and CATALOG[key] not in found:
            found.append(CATALOG[key])
    return found


def lookup(drug: str) -> Interaction | None:
    found = lookup_all(drug)
    return found[0] if found else None


def screen(drugs: Sequence[str]) -> tuple[list[Interaction], list[str]]:
    """Returns (interactions found, names not recognised). Unrecognised names
    must be shown to the clinician rather than silently passed."""
    found: list[Interaction] = []
    unknown: list[str] = []
    for d in drugs:
        hits = lookup_all(d)
        if hits:
            found += [h for h in hits if h not in found]
        elif " ".join(_ingredients(d)) not in KNOWN_NON_INTERACTING:
            unknown.append(d)
    return found, unknown


@dataclass
class WhatIf:
    interaction: Interaction
    before: Summary
    after: Summary
    p_above_1: float
    p_above_1_2: float
    adjusted: "Recommendation | None" = None   # dose that keeps the target with the new drug
    all_interactions: list = field(default_factory=list)

    def followup_tasks(self, start_h: float) -> list[tuple[float, str]]:
        """Check levels after starting, e.g. day 7 and day 28 for ACE inhibitors."""
        hits = self.all_interactions or [self.interaction]
        days = sorted({d for h in hits for d in h.followup_days})
        names = " + ".join(h.label for h in hits)
        return [(start_h + 24.0 * d, f"lithium level + U&E: day {d} after starting {names}") for d in days]


def what_if(post: Posterior, admins: Sequence[Administration], drug: str, t: float | None = None,
            n: int = 3000, rng: np.random.Generator | None = None, product=None,
            target: tuple[float, float] = (0.6, 0.8)) -> WhatIf | None:
    """Forecast this patient's 12-hour level if ``drug`` is added with no lithium
    change, and (given a ``product``) the lithium dose that would compensate."""
    hits = lookup_all(drug)
    if not hits:
        return None
    pk_hits = [h for h in hits if h.kind == "pk"]
    it = pk_hits[0] if pk_hits else hits[0]
    rng = rng or np.random.default_rng(11)
    if t is None:
        t = max([lv.time for lv in post.case.levels], default=0.0)
    daily = sum(a.mmol for a in admins)
    seed = int(rng.integers(1 << 31))
    mult = np.ones(n)
    for h in hits:
        mult = mult * h.sample(n, rng)
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
                  float(np.mean(after > 1.0)), float(np.mean(after > 1.2)), adjusted, hits)
