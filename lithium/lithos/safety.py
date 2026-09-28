"""Toxicity triage and elimination forecasting.

Thresholds here are illustrative defaults that a service must align with its
local toxicology pathway (in Australia, the Poisons Information Centre on
13 11 26). The distinctive part is the *forecast*: from the patient's own
posterior, how long will it take to fall below 1.0 mmol/L if lithium is
withheld? EXTRIP (Decker et al., CJASN 2015) lists an expected time to
[Li+] < 1.0 mEq/L of more than 36 h with optimal management as one criterion
for which extracorporeal treatment is suggested, so an individual forecast is
directly useful to the treating team.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .bayes import Posterior
from .forecast import Summary, summarise

SEVERE = {"seizure", "reduced consciousness", "coma", "arrhythmia", "hypotension"}
MODERATE = {"confusion", "ataxia", "slurred speech", "coarse tremor", "myoclonus", "muscle twitching",
            "drowsiness", "agitation", "hyperreflexia"}
MILD = {"worsening tremor", "nausea", "vomiting", "diarrhoea", "diarrhea", "weakness", "polyuria"}


@dataclass
class Triage:
    urgency: str                     # routine | review | same-day | emergency
    actions: list[str]
    extrip: list[str] = field(default_factory=list)
    hours_to_below_1: Summary | None = None


def hours_to_below(post: Posterior, t: float, observed: float, threshold: float = 1.0,
                   n: int = 400, max_h: float = 240.0, rng: np.random.Generator | None = None) -> Summary:
    """Posterior hours from ``t`` until the level falls below ``threshold`` with
    all further doses withheld. Each posterior curve is scaled to pass through
    the observed level, so the forecast is anchored to the measurement."""
    rng = rng or np.random.default_rng(3)
    grid = t + np.linspace(0.0, max_h, 481)
    doses = [d for d in post.case.doses if d.time <= t]
    U = post.sample(n, rng)
    out = np.empty(n)
    for i, u in enumerate(U):
        c = post.predict(grid, u, doses=doses)
        c = c * (observed / max(c[0], 1e-9))
        below = np.flatnonzero(c < threshold)
        out[i] = grid[below[0]] - t if below.size else max_h
    return summarise(out, "h", 0)


def triage(level: float, symptoms: set[str] | None = None, *, acute_ingestion: bool = False,
           egfr: float | None = None, post: Posterior | None = None, t: float | None = None) -> Triage:
    symptoms = {s.lower() for s in (symptoms or set())}
    severe, moderate, mild = symptoms & SEVERE, symptoms & MODERATE, symptoms & MILD
    actions: list[str] = []
    extrip: list[str] = []

    if severe or level >= 2.5 or (moderate and level >= 1.5):
        urgency = "emergency"
        actions += ["Stop lithium. Emergency department now (ambulance if drowsy, confused or seizing).",
                    "Discuss with toxicology / Poisons Information Centre.",
                    "Repeat level, U&E/creatinine, ECG; IV fluids per toxicology advice."]
    elif moderate or level >= 1.5:
        urgency = "same-day"
        actions += ["Withhold lithium. Same-day medical assessment.",
                    "Repeat level with U&E/creatinine within 24 h; look for the cause "
                    "(interaction, dehydration, renal decline, dose error)."]
    elif level > 1.2 or mild:
        urgency = "review"
        actions += ["Withhold the next dose(s) pending review; repeat the level and U&E within 24-48 h.",
                    "Check for precipitants: new NSAID/ACE inhibitor/ARB/diuretic, vomiting or diarrhoea, "
                    "reduced intake, heat, intercurrent illness."]
    elif level > 1.0:
        urgency = "review"
        actions += ["Above the usual maintenance ceiling: review dose and precipitants; repeat a level."]
    else:
        urgency = "routine"
        actions += ["No toxicity signal from the level alone; interpret with symptoms and trend."]

    if acute_ingestion:
        actions.append("Acute ingestion: early levels mislead (distribution); follow toxicology protocol "
                       "with serial levels. Sustained-release products may need whole-bowel irrigation.")

    renal_impaired = egfr is not None and egfr < 45
    if (renal_impaired and level > 4.0) or severe & {"reduced consciousness", "coma", "seizure", "arrhythmia"}:
        extrip.append("EXTRIP: extracorporeal treatment RECOMMENDED (impaired kidney function with "
                      "[Li+] > 4.0, or decreased consciousness, seizures or life-threatening dysrhythmia).")
    if level > 5.0 or "confusion" in symptoms:
        extrip.append("EXTRIP: extracorporeal treatment SUGGESTED ([Li+] > 5.0 or significant confusion).")

    hrs = None
    if post is not None and t is not None and level > 1.0:
        hrs = hours_to_below(post, t, level)
        actions.append(f"Forecast time to fall below 1.0 mmol/L with lithium withheld: {hrs} "
                       "(assumes current renal function; acute kidney injury prolongs this).")
        if hrs.median > 36 or hrs.hi > 36:
            extrip.append(f"EXTRIP: expected time to [Li+] < 1.0 may exceed 36 h (median {hrs.median:.0f} h, "
                          f"95th percentile {hrs.hi:.0f} h) - ECTR suggested criterion; discuss with toxicology.")
    return Triage(urgency, actions, extrip, hrs)
