"""External check of the engine's priors against published real-patient cohorts.

This is not individual-level validation (see validation/evaluate.py for that).
It asks a narrower question with data we can reach: for a patient resembling
each published cohort, does the engine's population prior reproduce what was
actually measured in real people? Every cohort value below is from the
verified research brief (docs/EVIDENCE.md, [E41] and related entries).

Where a paper's covariates are not reported, the matched patient is an
explicit assumption shown in the output, so the reader can judge the match.

    python -m validation.literature
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from lithos.bayes import individualise
from lithos.case import Case, CovariateTimeline
from lithos.models import get_engine
from lithos.pk import Dose, simulate, steady_state

RNG = np.random.default_rng(20260928)
N = 2000


def prior_params(sex, age, wt, ht, scr):
    """Draws of individual parameters from the ensemble prior for a matched patient."""
    cov = CovariateTimeline(sex=sex, age_at_t0=age, weight=wt, height=ht, creatinine=[(0.0, scr)])
    post = individualise(Case(cov, [], []), get_engine())
    return post.param_draws(N, 0.0, RNG), cov(0.0)


def summ(x):
    lo, med, hi = np.percentile(x, [5, 50, 95])
    return med, lo, hi


def vss(p):
    return p.v1 + p.v2


@dataclass
class Row:
    cohort: str
    quantity: str
    published: str
    engine: tuple[float, float, float]
    unit: str
    match: str
    verdict: str = ""
    band: tuple[float, float] = (0.0, 0.0)
    mid: float = 0.0

    @property
    def vs_mean(self) -> float:
        return self.engine[0] / self.mid - 1.0


def verdict(engine_med, pub_lo, pub_hi):
    """Is the engine's median inside the published band (mean +/- 1 SD, a reported
    range, or a stated tolerance when no SD was reported)?"""
    if pub_lo <= engine_med <= pub_hi:
        return "inside"
    rel = engine_med / pub_hi - 1 if engine_med > pub_hi else engine_med / pub_lo - 1
    return f"{abs(rel):.0%} {'above' if rel > 0 else 'below'}"


def run() -> list[Row]:
    rows: list[Row] = []
    _run(rows)
    for r in rows:
        lo, hi = r.band
        r.mid = MIDPOINTS.get((r.cohort, r.quantity), 0.5 * (lo + hi))
    return rows


# Published central values where the band is not symmetric around them.
MIDPOINTS = {
    ("IR vs CR (PMID 3089949)", "Time to peak, IR"): 1.44,
    ("IR vs CR (PMID 3089949)", "Time to peak, CR"): 4.25,
    ("IR vs CR (PMID 3089949)", "Cmax ratio CR/IR"): 0.78,
    ("Once-daily SR (Reddy 2014)", "12-h / 24-h level"): 1.37,
    ("Healthy volunteers (Stip 2001)", "12-h level on 1569 mg/day"): 0.705,
}


def _run(rows: list[Row]) -> None:

    # 1. Older adults on a morning dose (PMID 3110219): CL/F 15.6 +/- 4.0 mL/min, t1/2 26.9 +/- 5.5 h,
    #    V_area 0.64 +/- 0.16 L/kg, peak 0.83 +/- 0.25 mmol/L at 2.2 +/- 1.2 h, distribution complete 10.6 +/- 3.0 h
    ps, c = prior_params("F", 72, 68, 163, 85)
    match = "72 y woman, 68 kg, creatinine 85 umol/L (cohort covariates not reported)"
    cl = np.array([p.cl for p in ps]) * 1000 / 60
    rows.append(Row("Older adults (PMID 3110219)", "Clearance", "15.6 +/- 4.0", summ(cl), "mL/min", match,
                    verdict(np.median(cl), 11.6, 19.6), band=(11.6, 19.6)))
    th = np.array([p.terminal_half_life() for p in ps])
    rows.append(Row("Older adults (PMID 3110219)", "Terminal half-life", "26.9 +/- 5.5", summ(th), "h", match,
                    verdict(np.median(th), 21.4, 32.4), band=(21.4, 32.4)))
    varea = np.array([p.cl / (np.log(2) / p.terminal_half_life()) for p in ps]) / 68.0
    rows.append(Row("Older adults (PMID 3110219)", "V_area", "0.64 +/- 0.16", summ(varea), "L/kg", match,
                    verdict(np.median(varea), 0.48, 0.80), band=(0.48, 0.80)))
    dose_mmol = 0.21 * 68.0
    reg = [(8.0, dose_mmol, "IR")]
    grid = np.linspace(8.0, 32.0, 481)
    peaks, tmaxs, dist_done = [], [], []
    for p in ps[:600]:
        cc = steady_state(reg, p, grid % 24.0)
        i = int(np.argmax(cc))
        peaks.append(cc[i])
        tmaxs.append(grid[i] - 8.0)
        if p.two_compartment:
            k10, k12, k21 = p.cl / p.v1, p.q / p.v1, p.q / p.v2
            s = k10 + k12 + k21
            alpha = 0.5 * (s + np.sqrt(s * s - 4 * k10 * k21))
            dist_done.append(np.log(20) / alpha)       # time for the distribution phase to fall to 5%
    rows.append(Row("Older adults (PMID 3110219)", "Steady-state peak on 0.21 mmol/kg/day", "0.83 +/- 0.25",
                    summ(np.array(peaks)), "mmol/L", match, verdict(np.median(peaks), 0.58, 1.08), band=(0.58, 1.08)))
    rows.append(Row("Older adults (PMID 3110219)", "Time to peak", "2.2 +/- 1.2", summ(np.array(tmaxs)), "h", match,
                    verdict(np.median(tmaxs), 1.0, 3.4), band=(1.0, 3.4)))
    if dist_done:
        rows.append(Row("Older adults (PMID 3110219)", "Distribution phase ~complete (kidney model only)",
                        "10.6 +/- 3.0", summ(np.array(dist_done)), "h", match, verdict(np.median(dist_done), 7.6, 13.6), band=(7.6, 13.6)))

    # 2. Obese vs normal-weight (PMID 8162665): Vss 0.42 +/- 0.09 vs 0.66 +/- 0.16 L/kg; CL 33.9 +/- 7.0 vs 23.0 +/- 6.2 mL/min
    for label, sex, wt, ht, pub_v, pub_cl, v_rng, cl_rng in [
        ("Obese adults (PMID 8162665)", "M", 120.0, 175.0, "0.42 +/- 0.09", "33.9 +/- 7.0", (0.33, 0.51), (26.9, 40.9)),
        ("Normal-weight controls (PMID 8162665)", "M", 72.0, 176.0, "0.66 +/- 0.16", "23.0 +/- 6.2", (0.50, 0.82), (16.8, 29.2)),
    ]:
        ps, c = prior_params(sex, 35, wt, ht, 85)
        match = f"35 y man, {wt:.0f} kg, {ht:.0f} cm, creatinine 85 (cohort covariates not reported)"
        v = np.array([vss(p) for p in ps]) / wt
        rows.append(Row(label, "Vss", pub_v, summ(v), "L/kg", match, verdict(np.median(v), *v_rng), band=v_rng))
        cl = np.array([p.cl for p in ps]) * 1000 / 60
        rows.append(Row(label, "Clearance", pub_cl, summ(cl), "mL/min", match, verdict(np.median(cl), *cl_rng), band=cl_rng))

    # 3. Healthy adults (PMID 3118402): plasma clearance 22-43 mL/min uncorrected
    ps, c = prior_params("M", 30, 72, 178, 85)
    cl = np.array([p.cl for p in ps]) * 1000 / 60
    rows.append(Row("Healthy adults (PMID 3118402)", "Clearance", "22-43 (range)", summ(cl), "mL/min",
                    "30 y man, 72 kg, creatinine 85", verdict(np.median(cl), 22, 43), band=(22, 43)))

    # 4. Inpatients with schizophrenia (PMID 657687): Vd 0.79 +/- 0.34 L/kg; renal CL 24.4 +/- 8.0 mL/min; t1/2 28.9 +/- 7.9 h
    ps, c = prior_params("M", 35, 72, 176, 85)
    match = "35 y man, 72 kg, creatinine 85 (cohort covariates not reported)"
    v = np.array([vss(p) for p in ps]) / 72
    rows.append(Row("Inpatients (PMID 657687)", "Volume", "0.79 +/- 0.34", summ(v), "L/kg", match, verdict(np.median(v), 0.45, 1.13), band=(0.45, 1.13)))
    cl = np.array([p.cl for p in ps]) * 1000 / 60
    rows.append(Row("Inpatients (PMID 657687)", "Renal clearance", "24.4 +/- 8.0", summ(cl), "mL/min", match,
                    verdict(np.median(cl), 16.4, 32.4), band=(16.4, 32.4)))
    th = np.array([p.terminal_half_life() for p in ps])
    rows.append(Row("Inpatients (PMID 657687)", "Half-life", "28.9 +/- 7.9", summ(th), "h", match, verdict(np.median(th), 21.0, 36.8), band=(21.0, 36.8)))

    # 5. Manic-depressive patients, single dose (PMID 2079355): CL 33.2 +/- 15.5 mL/min; V 0.62 +/- 0.26 L/kg; t1/2 15.3 +/- 6.1 h
    ps, c = prior_params("F", 35, 70, 165, 70)
    match = "35 y woman, 70 kg, creatinine 70 (cohort covariates not reported)"
    cl = np.array([p.cl for p in ps]) * 1000 / 60
    rows.append(Row("Single-dose study (PMID 2079355)", "Clearance", "33.2 +/- 15.5", summ(cl), "mL/min", match,
                    verdict(np.median(cl), 17.7, 48.7), band=(17.7, 48.7)))
    th = np.array([p.terminal_half_life() for p in ps])
    rows.append(Row("Single-dose study (PMID 2079355)", "Half-life", "15.3 +/- 6.1", summ(th), "h", match,
                    verdict(np.median(th), 9.2, 21.4), band=(9.2, 21.4)))

    # 6. IR vs controlled release, single dose (PMID 3089949): tpeak 1.44 -> 4.25 h; CR Cmax 22% lower
    ps, c = prior_params("M", 35, 72, 176, 85)
    t = np.linspace(0, 24, 961)
    tir, tcr, ratio = [], [], []
    for p in ps[:600]:
        a = simulate([Dose(0.0, 10.0, "IR")], p, t)
        b = simulate([Dose(0.0, 10.0, "SR")], p, t)
        tir.append(t[np.argmax(a)]); tcr.append(t[np.argmax(b)]); ratio.append(b.max() / a.max())
    match = "35 y man, 72 kg; single dose"
    rows.append(Row("IR vs CR (PMID 3089949)", "Time to peak, IR", "1.44", summ(np.array(tir)), "h", match,
                    verdict(np.median(tir), 1.0, 2.0), band=(1.0, 2.0)))
    rows.append(Row("IR vs CR (PMID 3089949)", "Time to peak, CR", "4.25", summ(np.array(tcr)), "h", match,
                    verdict(np.median(tcr), 3.5, 5.0), band=(3.5, 5.0)))
    rows.append(Row("IR vs CR (PMID 3089949)", "Cmax ratio CR/IR", "0.78", summ(np.array(ratio)), "ratio", match,
                    verdict(np.median(ratio), 0.70, 0.86), band=(0.70, 0.86)))

    # 7. Once-daily SR, 12-h vs 24-h level (Reddy 2014, n=48): 0.82 vs 0.60 -> ratio 1.37
    ps, c = prior_params("F", 40, 70, 165, 75)
    r = []
    for p in ps[:600]:
        cc = steady_state([(21.0, 20.0, "SR")], p, [9.0, 21.0 - 1e-6])
        r.append(cc[0] / cc[1])
    rows.append(Row("Once-daily SR (Reddy 2014)", "12-h / 24-h level", "1.37 (0.82 vs 0.60)", summ(np.array(r)), "ratio",
                    "40 y woman, 70 kg; SR at 21:00", verdict(np.median(r), 1.25, 1.50), band=(1.25, 1.50)))

    # 8. Pepin-dosed healthy volunteers (Stip 2001): mean dose 1569 +/- 291 mg/day -> levels 0.67-0.74
    ps, c = prior_params("M", 25, 72, 178, 85)
    reg = [(8.0, 1569 / 2 * 2 / 73.891, "IR"), (20.0, 1569 / 2 * 2 / 73.891, "IR")]
    lv = np.array([steady_state(reg, p, [8.0])[0] for p in ps[:800]])
    rows.append(Row("Healthy volunteers (Stip 2001)", "12-h level on 1569 mg/day", "0.67-0.74 (weekly means)", summ(lv),
                    "mmol/L", "25 y man, 72 kg, creatinine 85; twice daily (regimen assumed)", verdict(np.median(lv), 0.60, 0.80), band=(0.60, 0.80)))
    return rows


def markdown(rows: list[Row]) -> str:
    L = ["| Cohort | Quantity | Published (real patients) | Band used | Engine prior: median (90% range) | "
         "Engine median vs published mean | Engine median in band? |",
         "|---|---|---|---|---|---|---|"]
    for r in rows:
        med, lo, hi = r.engine
        d = 2 if r.unit in ("mmol/L", "L/kg", "ratio") else 1
        L.append(f"| {r.cohort} | {r.quantity} | {r.published} {r.unit} | {r.band[0]:g}-{r.band[1]:g} | "
                 f"{med:.{d}f} ({lo:.{d}f}-{hi:.{d}f}) | {r.vs_mean:+.0%} | {r.verdict} |")
    inside = sum(r.verdict == "inside" for r in rows)
    within20 = sum(abs(r.vs_mean) <= 0.20 for r in rows)
    L += ["", f"Engine median inside the published band: {inside}/{len(rows)}. Within 20% of the published mean: "
              f"{within20}/{len(rows)}. Band = mean +/- 1 SD where an SD was reported, the reported range, or a stated "
              "tolerance (time to peak, ratios, Stip) where it was not.",
          "", "Matched patients (assumptions where the paper did not report covariates):", ""]
    seen = set()
    for r in rows:
        if r.cohort not in seen:
            L.append(f"- {r.cohort}: {r.match}")
            seen.add(r.cohort)
    return "\n".join(L)


if __name__ == "__main__":
    print(markdown(run()))
