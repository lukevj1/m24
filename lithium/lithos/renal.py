"""Kidney function: estimating equations and longitudinal trend analysis.

Lithium clearance tracks glomerular filtration, and chronic kidney disease is
the long-term harm clinicians fear most, so renal function is both a dosing
covariate and a safety outcome here.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

UMOL_PER_MG_DL_CREATININE = 88.42


def egfr_ckd_epi_2021(scr_umol: float, age: float, sex: str) -> float:
    """CKD-EPI 2021 (race-free) creatinine equation, mL/min/1.73 m^2."""
    scr = scr_umol / UMOL_PER_MG_DL_CREATININE
    female = sex.upper().startswith("F")
    kappa, alpha = (0.7, -0.241) if female else (0.9, -0.302)
    r = scr / kappa
    return 142.0 * min(r, 1.0) ** alpha * max(r, 1.0) ** -1.200 * 0.9938 ** age * (1.012 if female else 1.0)


def crcl_cockcroft_gault(scr_umol: float, age: float, sex: str, weight_kg: float) -> float:
    """Cockcroft-Gault creatinine clearance, mL/min (total body weight)."""
    scr = scr_umol / UMOL_PER_MG_DL_CREATININE
    crcl = (140.0 - age) * weight_kg / (72.0 * scr)
    return crcl * (0.85 if sex.upper().startswith("F") else 1.0)


def bsa_du_bois(weight_kg: float, height_cm: float) -> float:
    return 0.007184 * weight_kg ** 0.425 * height_cm ** 0.725


def fat_free_mass(weight_kg: float, height_cm: float, sex: str) -> float:
    """Janmahasatian et al. (2005) fat-free mass, kg."""
    bmi = weight_kg / (height_cm / 100.0) ** 2
    if sex.upper().startswith("F"):
        return 9270.0 * weight_kg / (8780.0 + 244.0 * bmi)
    return 9270.0 * weight_kg / (6680.0 + 216.0 * bmi)


def ckd_stage(egfr: float) -> str:
    for limit, stage in [(90, "G1"), (60, "G2"), (45, "G3a"), (30, "G3b"), (15, "G4")]:
        if egfr >= limit:
            return stage
    return "G5"


@dataclass
class RenalTrend:
    n: int
    years: float
    slope: float              # mL/min/1.73 m^2 per year (Theil-Sen)
    slope_lo: float           # approximate 90% interval
    slope_hi: float
    latest: float
    latest_stage: str
    flags: list[str]


def egfr_trend(times_years: np.ndarray, egfr: np.ndarray, age_now: float | None = None) -> RenalTrend:
    """Robust eGFR slope with clinically framed flags.

    Theil-Sen (median of pairwise slopes) resists the single spurious
    creatinine that plagues routine data. Flags follow widely used definitions:
    KDIGO rapid progression (> 5 mL/min/1.73 m^2 per year sustained) and NICE
    accelerated progression (>= 25 % fall with a category change, or a fall of
    >= 15 mL/min/1.73 m^2 within a year). Age-related decline after 40 is
    roughly 1 mL/min/1.73 m^2 per year, so anything steeper is worth a look.
    """
    t = np.asarray(times_years, dtype=float)
    y = np.asarray(egfr, dtype=float)
    order = np.argsort(t)
    t, y = t[order], y[order]
    flags: list[str] = []
    if t.size < 2 or np.ptp(t) < 0.25:
        latest = float(y[-1])
        return RenalTrend(int(t.size), float(np.ptp(t)) if t.size else 0.0, np.nan, np.nan, np.nan,
                          latest, ckd_stage(latest), ["insufficient data for a trend"])
    i, j = np.triu_indices(t.size, k=1)
    dt = t[j] - t[i]
    keep = dt > 1.0 / 52.0
    slopes = (y[j] - y[i])[keep] / dt[keep]
    slope = float(np.median(slopes))
    lo, hi = (float(v) for v in np.percentile(slopes, [5, 95]))
    latest = float(np.median(y[-2:])) if t.size >= 3 else float(y[-1])
    stage = ckd_stage(latest)

    if slope < -5.0 and hi < -1.0:
        flags.append("rapid decline (KDIGO: > 5 mL/min/1.73m2/yr)")
    elif slope < -2.0:
        flags.append("decline faster than expected for age (~1 mL/min/1.73m2/yr)")
    one_year = t >= t[-1] - 1.0
    if one_year.sum() >= 2:
        first, last = y[one_year][0], latest
        if last <= first - 15.0 or (last <= 0.75 * first and ckd_stage(first) != stage):
            flags.append("accelerated progression within 12 months (NICE definition)")
    if latest < 60:
        flags.append(f"eGFR {latest:.0f}: CKD {stage} - review target level, dosing frequency and monitoring interval")
    if latest < 45:
        flags.append("eGFR < 45: discuss with nephrology; weigh continuing vs alternatives with the patient")
    return RenalTrend(int(t.size), float(np.ptp(t)), slope, lo, hi, latest, stage, flags)
