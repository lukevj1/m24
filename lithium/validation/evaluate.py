"""Individual-level validation by iterative Bayesian forecasting.

Each lithium level is forecast from what was known before it was drawn: the
doses recorded up to the draw, all earlier levels, and (by default) only the
creatinine results available by the time of the previous level. The forecast
is then compared with what the laboratory measured. This mirrors how the tool
would be used and is the usual design for external validation of model-informed
precision dosing (e.g. Broeker et al., Clin Microbiol Infect 2019, vancomycin).

Comparators get exactly the same information:

    population    the engine's population prior alone (ignores earlier levels)
    carry         the previous level, unchanged
    proportional  the previous level scaled by the change in average daily dose
                  over the 72 h before each level (linear kinetics)

Reported per method, overall and by stratum: relative bias and relative RMSE,
mean absolute error, share within +/-20 % and within +/-0.1 mmol/L, and, for
model forecasts, coverage of the 90 % prediction interval. 95 % intervals come
from a patient-level bootstrap.

Only aggregate results are printed or written by default. Per-level output
(--rows) must stay where the data are: read docs/VALIDATION.md first.

    python -m validation.evaluate DATA_DIR [--out report.md] [--json metrics.json] [--procs 4]
"""

from __future__ import annotations

import argparse
import csv
import json
import multiprocessing as mp
from pathlib import Path

import numpy as np

from lithos.bayes import individualise
from lithos.case import Case
from lithos.models import get_engine
from lithos.pk import simulate

from .schema import Record, load

METHODS = [
    ("engine", "Engine: Bayesian forecast from earlier levels"),
    ("population", "Engine's population prior only"),
    ("carry", "Previous level, unchanged"),
    ("proportional", "Previous level x change in daily dose"),
]
ASSUMED_WINDOW_H = 120.0   # levels with assumed (undocumented) doses in the 5 days before are reported separately


def forecast(post, doses, t: float, timing_sd: float, n: int, rng: np.random.Generator) -> dict:
    """Forecast of a level drawn at ``t``: median individual prediction, and a 90 %
    prediction interval for the measured value (includes the residual error model)."""
    model = post.model
    f, df = np.empty(n), np.empty(n)
    for i, p in enumerate(post.param_draws(n, t, rng)):
        c, dc = simulate(doses, p, [t], derivative=True)
        f[i], df[i] = c[0], dc[0]
    before = [d.time for d in doses if d.time <= t]
    since = t - max(before) if before else np.inf
    dist = model.dist_phase_sd * np.exp(-since / model.dist_phase_tau_h)
    var = model.sigma_add ** 2 + (model.sigma_prop ** 2 + dist ** 2) * f * f + (df * timing_sd) ** 2
    y = f + np.sqrt(var) * rng.standard_normal(n)
    return {"med": float(np.median(f)), "lo": float(np.percentile(y, 5)), "hi": float(np.percentile(y, 95)),
            "p_gt1": float(np.mean(y > 1.0))}


def _mean_daily(doses, t: float, window_h: float = 72.0) -> float:
    return sum(d.amount for d in doses if t - window_h <= d.time < t) * 24.0 / window_h


def evaluate_record(rec: Record, engine, *, n_draws: int = 300, covariates: str = "prospective",
                    lookback_days: float = 365.0, seed: int = 0) -> list[dict]:
    rng = np.random.default_rng(seed)
    levels = sorted(rec.levels, key=lambda lv: lv.time)
    rows = []
    for k, lv in enumerate(levels):
        t = lv.time
        prev = [p for p in levels[:k] if t - lookback_days * 24.0 <= p.time < t - 1e-9]
        start = (prev[0].time if prev else t) - 30 * 24.0
        pick = [(d, s) for d, s in zip(rec.doses, rec.dose_source) if start <= d.time < t]
        doses = [d for d, _ in pick]
        row = {"patient": rec.patient_id, "k": k, "obs": lv.value, "n_prior": len(prev),
               "time_quality": lv.note or "exact"}
        if not any(d.time >= t - 48.0 for d in doses):
            row["skip"] = "no recorded dose in the 48 h before the level"
            rows.append(row)
            continue
        t_cov = (prev[-1].time + 1.0 if prev else t - 24.0) if covariates == "prospective" else t + 1.0
        tl, cr_ok = rec.covariates(t_cov)
        c_now = tl(t)
        last = max(doses, key=lambda d: d.time)
        row.update(hours_since_dose=t - last.time, formulation=last.formulation, age=c_now.age,
                   egfr=c_now.egfr, cr_available=cr_ok, body_imputed=rec.body_size_imputed,
                   assumed_5d=any(s == "assumed" and d.time >= t - ASSUMED_WINDOW_H for d, s in pick))
        try:
            pop = individualise(Case(tl, doses, []), engine)
            row["population"] = forecast(pop, doses, t, lv.timing_sd or 0.25, n_draws, rng)
            if prev:
                post = individualise(Case(tl, doses, prev), engine)
                row["engine"] = forecast(post, doses, t, lv.timing_sd or 0.25, n_draws, rng)
                row["converged"] = bool(post.converged)
                p = prev[-1]
                row["carry"] = p.value
                d0, d1 = _mean_daily(doses, p.time), _mean_daily(doses, t)
                row["proportional"] = p.value * d1 / d0 if d0 > 0 else None
                row["dose_changed"] = d0 > 0 and abs(d1 / d0 - 1.0) > 0.10
                row["days_since_prev"] = (t - p.time) / 24.0
            else:
                row["engine"] = row["population"]
        except Exception as exc:  # report, never hide
            row["skip"] = f"engine error: {type(exc).__name__}"
        rows.append(row)
    return rows


def _worker(args):
    rec, engine_name, n_draws, covariates, seed = args
    return evaluate_record(rec, get_engine(engine_name), n_draws=n_draws, covariates=covariates, seed=seed)


def run(records: list[Record], *, engine: str = "ensemble", n_draws: int = 300, covariates: str = "prospective",
        procs: int = 1) -> list[dict]:
    jobs = [(r, engine, n_draws, covariates, 1000 + i) for i, r in enumerate(records)]
    if procs <= 1:
        parts = [_worker(j) for j in jobs]
    else:
        with mp.Pool(procs) as pool:
            parts = pool.map(_worker, jobs, chunksize=1)
    return [row for part in parts for row in part]


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def _pred(row: dict, method: str):
    v = row.get(method)
    return v["med"] if isinstance(v, dict) else v


def metrics(obs, pred, lo=None, hi=None) -> dict:
    obs, pred = np.asarray(obs, float), np.asarray(pred, float)
    if obs.size == 0:
        return {"n": 0}
    err = pred - obs
    rel = err / obs
    m = {"n": int(obs.size), "rbias": float(100 * rel.mean()), "rrmse": float(100 * np.sqrt(np.mean(rel ** 2))),
         "mae": float(np.abs(err).mean()), "within20": float(100 * np.mean(np.abs(rel) <= 0.2)),
         "within01": float(100 * np.mean(np.abs(err) <= 0.1 + 1e-9))}
    if lo is not None:
        m["coverage90"] = float(100 * np.mean((obs >= np.asarray(lo)) & (obs <= np.asarray(hi))))
    return m


def _method_metrics(rows: list[dict], method: str) -> dict:
    sub = [r for r in rows if _pred(r, method) is not None]
    obs = [r["obs"] for r in sub]
    pred = [_pred(r, method) for r in sub]
    if method in ("engine", "population"):
        return metrics(obs, pred, [r[method]["lo"] for r in sub], [r[method]["hi"] for r in sub])
    return metrics(obs, pred)


def _bootstrap(rows: list[dict], fn, B: int, rng: np.random.Generator):
    """95 % percentile intervals of the vector fn(rows), resampling patients."""
    by: dict[str, list[dict]] = {}
    for r in rows:
        by.setdefault(r["patient"], []).append(r)
    ids = list(by)
    if len(ids) < 5:
        return None
    vals = np.array([fn([r for i in rng.choice(len(ids), size=len(ids)) for r in by[ids[i]]]) for _ in range(B)],
                    dtype=float)
    lo, hi = np.nanpercentile(vals, [2.5, 97.5], axis=0)
    return [[float(a), float(b)] for a, b in zip(np.atleast_1d(lo), np.atleast_1d(hi))]


def summarise(rows: list[dict], boot: int = 500, seed: int = 1) -> dict:
    rng = np.random.default_rng(seed)
    ok = [r for r in rows if "skip" not in r]
    skipped: dict[str, int] = {}
    for r in rows:
        if "skip" in r:
            skipped[r["skip"]] = skipped.get(r["skip"], 0) + 1
    primary = [r for r in ok if not r["assumed_5d"]]
    follow = [r for r in primary if r["n_prior"] >= 1 and r.get("proportional") is not None]
    first = [r for r in primary if r["n_prior"] == 0]

    def table(sub):
        out = {}
        for m, _ in METHODS:
            mm = _method_metrics(sub, m)
            if mm.get("n", 0) and len({r["patient"] for r in sub}) >= 5:
                ci = _bootstrap(sub, lambda rs, m=m: [(_method_metrics(rs, m).get(k, np.nan))
                                                      for k in ("mae", "within20")], boot, rng)
                if ci is not None:
                    mm["mae_ci"], mm["within20_ci"] = ci[0], ci[1]
            out[m] = mm
        return out

    res = {
        "patients": len({r["patient"] for r in rows}),
        "levels": len(rows),
        "evaluable": len(ok),
        "skipped": skipped,
        "with_assumed_doses": len(ok) - len(primary),
        "follow_up": table(follow),
        "first_level": {"population": _method_metrics(first, "population")},
    }
    if len({r["patient"] for r in follow}) >= 5:
        def diff(rs):
            e, c = _method_metrics(rs, "engine"), _method_metrics(rs, "carry")
            return [e["mae"] - c["mae"], e["within20"] - c["within20"]]
        ci = _bootstrap(follow, diff, boot, rng)
        d = diff(follow)
        res["engine_minus_carry"] = {"mae": d[0], "mae_ci": ci[0], "within20": d[1], "within20_ci": ci[1]}

    strata = {
        "Hours since last dose": lambda r: _band(r["hours_since_dose"], [6, 10, 14, 24], "h"),
        "Earlier levels used": lambda r: _band(r["n_prior"], [1, 2, 4], ""),
        "Age": lambda r: "< 65" if r["age"] < 65 else ">= 65",
        "eGFR (mL/min/1.73m2)": lambda r: _band(r["egfr"], [45, 60, 90], ""),
        "Formulation of last dose": lambda r: r["formulation"],
        "Dose changed since previous level": lambda r: "n/a" if r["n_prior"] == 0 else
        ("yes" if r.get("dose_changed") else "no"),
        "Days since previous level": lambda r: "n/a" if r["n_prior"] == 0 else _band(r["days_since_prev"], [8, 91], "d"),
        "Sample time recorded as": lambda r: r["time_quality"],
    }
    res["strata"] = {}
    for name, key in strata.items():
        groups: dict[str, list[dict]] = {}
        for r in primary:
            groups.setdefault(key(r), []).append(r)
        res["strata"][name] = {g: {m: _method_metrics(sub, m) for m in ("engine", "population", "carry")}
                               for g, sub in sorted(groups.items(), key=lambda kv: _order(kv[0]))}

    hi = [r for r in follow if r["obs"] > 1.0]
    lo_ = [r for r in follow if r["obs"] <= 1.0]
    res["high_levels"] = {
        "n_above_1": len(hi),
        "flagged_pct": float(100 * np.mean([r["engine"]["p_gt1"] >= 0.2 for r in hi])) if hi else None,
        "false_flag_pct": float(100 * np.mean([r["engine"]["p_gt1"] >= 0.2 for r in lo_])) if lo_ else None,
    }
    assumed = [r for r in ok if r["assumed_5d"] and r["n_prior"] >= 1]
    res["assumed_dose_levels"] = {m: _method_metrics(assumed, m) for m in ("engine", "population", "carry")}
    res["not_converged"] = sum(1 for r in follow if r.get("converged") is False)
    return res


def _band(x: float, cuts: list[float], unit: str) -> str:
    edges = [-np.inf] + list(cuts) + [np.inf]
    for a, b in zip(edges[:-1], edges[1:]):
        if a <= x < b:
            if a == -np.inf:
                return f"< {b:g}{unit}" if unit else f"< {b:g}"
            if b == np.inf:
                return f">= {a:g}{unit}" if unit else f">= {a:g}"
            return f"{a:g}-{b:g}{unit}" if unit else f"{a:g}-{b:g}"
    return "?"


def _order(label: str):
    import re
    m = re.search(r"-?\d+(\.\d+)?", label)
    lead = 0 if label.startswith("<") else (2 if label.startswith(">") else 1)
    return (label == "n/a", float(m.group()) if m else 0.0, lead, label)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def _fmt(m: dict, key: str, d: int = 1, pct: bool = False) -> str:
    if not m.get("n"):
        return "-"
    v = m.get(key)
    if v is None:
        return "-"
    s = f"{v:.{d}f}{'%' if pct else ''}"
    ci = m.get(key + "_ci")
    if ci is not None and all(np.isfinite(ci)):
        s += f" ({ci[0]:.{d}f}-{ci[1]:.{d}f})"
    return s


def report_md(res: dict, title: str = "Individual-level validation") -> str:
    L = [f"# {title}", "",
         f"Patients (or admissions): {res['patients']}. Levels: {res['levels']}; evaluable: {res['evaluable']}; "
         f"with undocumented (assumed) doses in the 5 days before, reported separately: {res['with_assumed_doses']}.", ""]
    if res["skipped"]:
        L += ["Not evaluable: " + "; ".join(f"{k}: {v}" for k, v in res["skipped"].items()) + ".", ""]
    fu = res["follow_up"]
    n = fu.get("engine", {}).get("n", 0)
    L += [f"## Forecasting the next level (levels with at least one earlier level, n = {n})", "",
          "| Method | Relative bias | Relative RMSE | Mean abs. error, mmol/L (95% CI) | Within +/-20% (95% CI) | "
          "Within +/-0.1 mmol/L | 90% PI coverage |", "|---|---|---|---|---|---|---|"]
    for m, label in METHODS:
        mm = fu.get(m, {})
        L.append(f"| {label} | {_fmt(mm, 'rbias', 1, True)} | {_fmt(mm, 'rrmse', 1, True)} | {_fmt(mm, 'mae', 3)} | "
                 f"{_fmt(mm, 'within20', 0, True)} | {_fmt(mm, 'within01', 0, True)} | {_fmt(mm, 'coverage90', 0, True)} |")
    if "engine_minus_carry" in res:
        d = res["engine_minus_carry"]
        L += ["", f"Engine minus 'previous level, unchanged': mean absolute error {d['mae']:+.3f} mmol/L "
                  f"(95% CI {d['mae_ci'][0]:+.3f} to {d['mae_ci'][1]:+.3f}); within +/-20% {d['within20']:+.0f} "
                  f"points (95% CI {d['within20_ci'][0]:+.0f} to {d['within20_ci'][1]:+.0f})."]
    fl = res["first_level"]["population"]
    L += ["", f"## First level of each patient (population prior only, n = {fl.get('n', 0)})", "",
          f"Relative bias {_fmt(fl, 'rbias', 1, True)}, relative RMSE {_fmt(fl, 'rrmse', 1, True)}, within +/-20% "
          f"{_fmt(fl, 'within20', 0, True)}, 90% PI coverage {_fmt(fl, 'coverage90', 0, True)}.", ""]
    hl = res["high_levels"]
    if hl["n_above_1"]:
        L += [f"Levels above 1.0 mmol/L among forecast levels: {hl['n_above_1']}; the engine gave P(> 1.0) >= 20% "
              f"beforehand for {hl['flagged_pct']:.0f}% of them (and for {hl['false_flag_pct']:.0f}% of levels "
              "that stayed <= 1.0).", ""]
    L += ["## By stratum (all evaluable levels without assumed doses)", "",
          "Mean absolute error in mmol/L, with share within +/-20% and 90% PI coverage for the engine.", ""]
    for name, groups in res["strata"].items():
        L += [f"**{name}**", "", "| Group | n | Engine MAE | Engine within 20% | Engine PI coverage | "
              "Population MAE | Previous level MAE |", "|---|---|---|---|---|---|---|"]
        for g, ms in groups.items():
            e, p, c = ms["engine"], ms["population"], ms["carry"]
            L.append(f"| {g} | {e.get('n', 0)} | {_fmt(e, 'mae', 3)} | {_fmt(e, 'within20', 0, True)} | "
                     f"{_fmt(e, 'coverage90', 0, True)} | {_fmt(p, 'mae', 3)} | {_fmt(c, 'mae', 3)} |")
        L.append("")
    a = res["assumed_dose_levels"]
    if a["engine"].get("n"):
        L += [f"Levels with assumed doses in the previous 5 days (n = {a['engine']['n']}): engine MAE "
              f"{_fmt(a['engine'], 'mae', 3)}, population {_fmt(a['population'], 'mae', 3)}, previous level "
              f"{_fmt(a['carry'], 'mae', 3)}.", ""]
    if res["not_converged"]:
        L += [f"Fits flagged as not converged: {res['not_converged']}.", ""]
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("data", help="directory in the canonical format (validation/schema.py)")
    ap.add_argument("--engine", default="ensemble")
    ap.add_argument("--draws", type=int, default=300)
    ap.add_argument("--covariates", choices=["prospective", "concurrent"], default="prospective",
                    help="prospective: only creatinine available by the previous level (default)")
    ap.add_argument("--procs", type=int, default=max(1, mp.cpu_count() - 1))
    ap.add_argument("--limit", type=int, default=0, help="evaluate only the first N patients (quick check)")
    ap.add_argument("--out", help="write the Markdown report here")
    ap.add_argument("--json", help="write aggregate metrics (JSON) here")
    ap.add_argument("--rows", help="write per-level predictions (CSV). Keep with the data; do not share.")
    a = ap.parse_args()
    records = load(a.data)
    if a.limit:
        records = records[: a.limit]
    rows = run(records, engine=a.engine, n_draws=a.draws, covariates=a.covariates, procs=a.procs)
    res = summarise(rows)
    md = report_md(res)
    print(md)
    if a.out:
        Path(a.out).write_text(md + "\n")
    if a.json:
        Path(a.json).write_text(json.dumps(res, indent=1, default=float))
    if a.rows:
        flat = []
        for r in rows:
            f = {k: v for k, v in r.items() if not isinstance(v, dict)}
            for m in ("engine", "population"):
                if isinstance(r.get(m), dict):
                    f.update({f"{m}_{k}": v for k, v in r[m].items()})
            flat.append(f)
        keys = sorted({k for f in flat for k in f})
        with open(a.rows, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(flat)


if __name__ == "__main__":
    main()
