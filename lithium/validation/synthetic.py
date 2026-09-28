"""A synthetic dataset in the canonical format, for checking the pipeline.

Virtual outpatients (from the deliberately misspecified truth model used in
sim/) are started on lithium, titrated the usual way (level read at face
value, +/-250 mg), then followed with three-monthly levels for a year.

What the dataset records is what a clinic would record: the prescribed
regimen (not the missed doses), levels with assay noise, sample times that
are sometimes only approximate, and creatinine with each level. Some patients
start an unrecorded thiazide part-way through. Results on this dataset test
the plumbing and the engine's behaviour; they are not validation.

    python -m validation.synthetic OUT_DIR [--n 60]
    python -m validation.evaluate OUT_DIR
"""

from __future__ import annotations

import argparse
from datetime import datetime, timedelta

import numpy as np

from lithos.case import Administration, expand
from lithos.models import get_model
from sim.population import make_population, measure, true_level
from sim.trials import realistic_hours_after_evening_dose

from .schema import write

CLOCK = 21.0


def make_dataset(out_dir, n: int = 60, seed: int = 11, days: int = 365, thiazide_fraction: float = 0.2) -> None:
    truth = get_model("truth-perturbed")
    pop = make_population(n, truth, seed=seed, horizon_days=days + 30)
    rng = np.random.default_rng(seed + 1)
    base = datetime(2150, 1, 1)
    stamp = lambda h: (base + timedelta(hours=float(h))).strftime("%Y-%m-%d %H:%M")  # noqa: E731
    tables: dict[str, list[dict]] = {k: [] for k in ("patients", "doses", "levels", "creatinine")}
    level_days = [7, 14, 21, 28] + list(range(91, days, 91))
    for vp in pop:
        pid = f"S{vp.pid:04d}"
        c0 = vp.cov(0.0)
        tables["patients"].append({"patient_id": pid, "sex": c0.sex, "age_years": round(c0.age, 1),
                                   "height_cm": round(c0.height), "weight_kg": round(c0.weight, 1)})
        t_thz = None
        if rng.random() < thiazide_fraction:
            t_thz = float(rng.uniform(120, days - 60)) * 24.0
        m_thz = float(np.exp(rng.normal(np.log(0.75), 0.12)))
        cl_fn = (lambda s, t_thz=t_thz, m=m_thz: m if (t_thz is not None and s >= t_thz) else 1.0)
        mg = 250.0 if c0.age >= 65 else 500.0
        seg_start, taken = CLOCK, []

        def give(end_h: float) -> None:
            nonlocal seg_start
            seg = expand(seg_start, end_h, [Administration(CLOCK, mg)])
            taken.extend(vp.take(seg, rng))
            for d in seg:
                tables["doses"].append({"patient_id": pid, "time": stamp(d.time), "mg": mg, "salt": "carbonate",
                                        "formulation": "IR", "source": "prescribed"})
            seg_start = end_h

        for day in level_days:
            t = (day - 1) * 24.0 + CLOCK + realistic_hours_after_evening_dose(rng)
            give(t)
            val = measure(true_level(vp, truth, taken, t, cl_mult_fn=cl_fn), rng)
            exact = rng.random() < 0.7
            t_rec = t if exact else t + float(np.round(rng.normal(0, 0.75) * 4) / 4)
            tables["levels"].append({"patient_id": pid, "time": stamp(t_rec), "value_mmol_l": val,
                                     "time_quality": "exact" if exact else "approximate"})
            tables["creatinine"].append({"patient_id": pid, "time": stamp(t),
                                         "umol_l": round(c0.scr * float(np.exp(rng.normal(0, 0.05))))})
            new = mg
            lo, hi = (0.6, 0.8) if day <= 28 else (0.5, 1.0)
            if val < lo:
                new = mg + 250.0
            elif val > hi:
                new = max(250.0, mg - 250.0)
            if new != mg:
                change = (np.floor(t / 24.0) + 1) * 24.0 + CLOCK   # result next day, new dose that evening
                give(change)
                mg = new
        give(days * 24.0)
    write(out_dir, tables)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--n", type=int, default=60)
    ap.add_argument("--seed", type=int, default=11)
    a = ap.parse_args()
    make_dataset(a.out, n=a.n, seed=a.seed)
    print(f"wrote a synthetic dataset of {a.n} patients to {a.out}")


if __name__ == "__main__":
    main()
