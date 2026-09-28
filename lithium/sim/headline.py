"""Headline numbers and chart data from sim/results.json.

    python -m sim.headline            # prints the headline table (markdown)
    python -m sim.headline --json F   # writes chart data for the concept page
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ARMS = [("uc_step", "Usual care, stepwise"), ("uc_prop", "Usual care, proportional"), ("mipd", "Model-informed")]


def load():
    d = json.loads((ROOT / "sim" / "results.json").read_text())
    return d["meta"], [r for r in d["rows"] if "error" not in r]


def titration(rows):
    out = {}
    for key, _ in ARMS:
        xs = [r[key] for r in rows]
        tests = np.array([x["tests"] for x in xs], float)
        conf = [x["confirmed_day"] for x in xs if x["confirmed_day"] is not None]
        ther = [x["therapeutic_day"] for x in xs if x["therapeutic_day"] is not None]
        out[key] = {
            "tests_mean": float(tests.mean()),
            "tests_med": float(np.median(tests)),
            "tests_iqr": [float(v) for v in np.percentile(tests, [25, 75])],
            "confirmed_pct": 100 * float(np.mean([x["confirmed_day"] is not None for x in xs])),
            "confirmed_day_med": float(np.median(conf)) if conf else float("nan"),
            "in_range_pct": 100 * float(np.mean([x["final_in_range"] for x in xs])),
            "in_wide_pct": 100 * float(np.mean([x["final_in_wide_range"] for x in xs])),
            "over1_pct": 100 * float(np.mean([x["days_on_dose_over_1"] > 0 for x in xs])),
            "therapeutic_day_med": float(np.median(ther)) if ther else float("nan"),
            "therapeutic_pct": 100 * float(np.mean([x["therapeutic_day"] is not None for x in xs])),
            "changes_med": float(np.median([x["n_dose_changes"] for x in xs])),
        }
    return out


def anytime(rows):
    at = [r["anytime"] for r in rows if "anytime" in r]
    cat = lambda v: 0 if v < 0.6 else (2 if v > 0.8 else 1)
    bins = [(2, 6), (6, 10), (10, 14), (14, 18), (18, 22.01)]
    per_bin = []
    for lo, hi in bins:
        sub = [a for a in at if lo <= a["hours"] < hi]
        per_bin.append({
            "mid": (lo + min(hi, 22)) / 2, "n": len(sub),
            "naive": float(np.median([abs(a["naive"] - a["truth12"]) for a in sub])),
            "engine": float(np.median([abs(a["engine_med"] - a["truth12"]) for a in sub])),
            "history": float(np.median([abs(a["hist_med"] - a["truth12"]) for a in sub])),
            "naive_wrong": 100 * float(np.mean([cat(a["naive"]) != cat(a["truth12"]) for a in sub])),
            "engine_wrong": 100 * float(np.mean([cat(a["engine_med"]) != cat(a["truth12"]) for a in sub])),
            "coverage": 100 * float(np.mean([a["engine_lo"] <= a["truth12"] <= a["engine_hi"] for a in sub])),
        })
    allm = {
        "naive": float(np.median([abs(a["naive"] - a["truth12"]) for a in at])),
        "engine": float(np.median([abs(a["engine_med"] - a["truth12"]) for a in at])),
        "naive_wrong": 100 * float(np.mean([cat(a["naive"]) != cat(a["truth12"]) for a in at])),
        "engine_wrong": 100 * float(np.mean([cat(a["engine_med"]) != cat(a["truth12"]) for a in at])),
        "coverage": 100 * float(np.mean([a["engine_lo"] <= a["truth12"] <= a["engine_hi"] for a in at])),
        "n": len(at),
    }
    return per_bin, allm


def thiazide(rows):
    th = [r["thiazide"] for r in rows if "thiazide" in r]
    na = np.array([t["no_action"] for t in th])
    ea = np.array([t["engine_adjusted"] for t in th])
    f = lambda a, lo, hi: 100 * float(np.mean((a >= lo) & (a <= hi)))
    return {
        "n": len(th),
        "over1": [100 * float(np.mean(na > 1.0)), 100 * float(np.mean(ea > 1.0))],
        "over12": [100 * float(np.mean(na > 1.2)), 100 * float(np.mean(ea > 1.2))],
        "wide": [f(na, 0.5, 0.9), f(ea, 0.5, 0.9)],
        "target": [f(na, 0.6, 0.8), f(ea, 0.6, 0.8)],
        "under05": [100 * float(np.mean(na < 0.5)), 100 * float(np.mean(ea < 0.5))],
        "coverage": 100 * float(np.mean([t["pred_lo"] <= t["no_action"] <= t["pred_hi"] for t in th])),
        "base_med": float(np.median([t["base"] for t in th])),
    }


def early(rows):
    ea = [e for r in rows if "early" in r for e in r["early"]]
    days = sorted({e["day"] for e in ea})
    res = []
    for d in days:
        sub = [e for e in ea if e["day"] == d]
        tr = np.array([e["truth"] for e in sub])
        res.append({"day": d,
                    "naive": float(np.median(np.abs(np.array([e["naive"] for e in sub]) - tr))),
                    "no_history": float(np.median(np.abs(np.array([e["no_history"] for e in sub]) - tr))),
                    "with_history": float(np.median(np.abs(np.array([e["with_history"] for e in sub]) - tr)))})
    return res


def web(meta, rows) -> dict:
    T = titration(rows)
    bins, allm = anytime(rows)
    th = thiazide(rows)
    ea = early(rows)
    lab = dict(ARMS)
    tests_rows = [{"title": None, "bars": [
        {"label": lab[k], "v": T[k]["tests_mean"], "emph": k == "mipd",
         "note": f"median {T[k]['tests_med']:.0f} (IQR {T[k]['tests_iqr'][0]:.0f}-{T[k]['tests_iqr'][1]:.0f})"}
        for k, _ in ARMS]}]
    pct_rows = [
        {"title": "Final dose truly in 0.6–0.8", "bars": [{"label": lab[k], "v": T[k]["in_range_pct"], "emph": k == "mipd"} for k, _ in ARMS]},
        {"title": "Final dose truly in 0.5–0.9", "bars": [{"label": lab[k], "v": T[k]["in_wide_pct"], "emph": k == "mipd"} for k, _ in ARMS]},
        {"title": "Ever on a dose with a true level above 1.0", "bars": [{"label": lab[k], "v": T[k]["over1_pct"], "emph": k == "mipd"} for k, _ in ARMS]},
    ]
    return {
        "intro": (f"{meta['n']} virtual adults (a fifth aged over 65) were generated from a deliberately different "
                  "model than the engine uses. They miss doses, have blood drawn at realistic times and have noisy "
                  "assays. Usual care reads each level at face value as a 12-hour level."),
        "cards": [
            {"kind": "bars", "title": "Starting lithium: blood tests until a stable dose is confirmed",
             "sub": "Mean tests over an 8-week titration to 0.6–0.8 mmol/L; hover for median and IQR.",
             "rows": tests_rows, "labelW": 170,
             "table": {"headers": ["Arm", "Mean tests", "Median (IQR)", "Confirmed stable", "Days to confirmed (median)"],
                       "rows": [[lab[k], f"{T[k]['tests_mean']:.1f}", f"{T[k]['tests_med']:.0f} ({T[k]['tests_iqr'][0]:.0f}-{T[k]['tests_iqr'][1]:.0f})",
                                 f"{T[k]['confirmed_pct']:.0f}%", f"{T[k]['confirmed_day_med']:.0f}"] for k, _ in ARMS]}},
            {"kind": "bars", "pctFmt": True, "max": 100, "title": "Starting lithium: where the final dose really landed",
             "sub": "Judged against each virtual patient's true 12-hour level, which no arm can see.",
             "rows": pct_rows, "labelW": 170,
             "table": {"headers": ["Arm", "Truly 0.6–0.8", "Truly 0.5–0.9", "Ever > 1.0"],
                       "rows": [[lab[k], f"{T[k]['in_range_pct']:.0f}%", f"{T[k]['in_wide_pct']:.0f}%", f"{T[k]['over1_pct']:.0f}%"] for k, _ in ARMS]}},
            {"kind": "line", "title": "A level drawn at a convenient time",
             "sub": "Median error against what a correctly timed 12-hour level would have shown, by hours after the evening dose.",
             "series": [{"short": "As measured", "emph": False, "pts": [[b["mid"], b["naive"]] for b in bins]},
                        {"short": "Engine", "emph": True, "pts": [[b["mid"], b["engine"]] for b in bins]}],
             "xticks": [b["mid"] for b in bins], "xunit": "h", "xlabel": "hours after the evening dose", "band": [11, 13],
             "table": {"headers": ["Hours", "n", "As measured", "Engine", "Wrong category (as measured)", "Wrong category (engine)", "90% interval coverage"],
                       "rows": [[f"{b['mid']:.0f}", str(b["n"]), f"{b['naive']:.3f}", f"{b['engine']:.3f}", f"{b['naive_wrong']:.0f}%", f"{b['engine_wrong']:.0f}%", f"{b['coverage']:.0f}%"] for b in bins]}},
            {"kind": "bars", "pctFmt": True, "max": 100, "title": "A thiazide is started in a stable patient",
             "sub": "True 12-hour level afterwards: no lithium change vs the engine's pre-emptive adjustment.",
             "rows": [
                 {"title": "Above 1.0 mmol/L", "bars": [{"label": "No change", "v": th["over1"][0]}, {"label": "Engine-adjusted", "v": th["over1"][1], "emph": True}]},
                 {"title": "Above 1.2 mmol/L", "bars": [{"label": "No change", "v": th["over12"][0]}, {"label": "Engine-adjusted", "v": th["over12"][1], "emph": True}]},
                 {"title": "Within 0.5–0.9 mmol/L", "bars": [{"label": "No change", "v": th["wide"][0]}, {"label": "Engine-adjusted", "v": th["wide"][1], "emph": True}]},
             ], "labelW": 130,
             "table": {"headers": ["Outcome", "No change", "Engine-adjusted"],
                       "rows": [["> 1.0", f"{th['over1'][0]:.0f}%", f"{th['over1'][1]:.0f}%"], ["> 1.2", f"{th['over12'][0]:.0f}%", f"{th['over12'][1]:.0f}%"],
                                ["0.5–0.9", f"{th['wide'][0]:.0f}%", f"{th['wide'][1]:.0f}%"], ["0.6–0.8", f"{th['target'][0]:.0f}%", f"{th['target'][1]:.0f}%"],
                                ["< 0.5", f"{th['under05'][0]:.0f}%", f"{th['under05'][1]:.0f}%"]]}},
            {"kind": "line", "title": "How soon after a dose change is a level useful?",
             "sub": "Median error in the new steady-state 12-hour level, by day of sampling, with one earlier level on record.",
             "series": [{"short": "As measured", "emph": False, "pts": [[e["day"], e["naive"]] for e in ea]},
                        {"short": "Engine", "emph": True, "pts": [[e["day"], e["with_history"]] for e in ea]}],
             "xticks": [e["day"] for e in ea], "xunit": "day", "xlabel": "day of sampling after the change",
             "table": {"headers": ["Day", "As measured", "Engine, no history", "Engine + earlier level"],
                       "rows": [[str(e["day"]), f"{e['naive']:.3f}", f"{e['no_history']:.3f}", f"{e['with_history']:.3f}"] for e in ea]}},
        ],
    }


def headline_md(meta, rows) -> str:
    T = titration(rows)
    bins, allm = anytime(rows)
    th = thiazide(rows)
    ea = early(rows)
    u, p, m = T["uc_step"], T["uc_prop"], T["mipd"]
    d1 = next(e for e in ea if e["day"] == 1)
    d7 = next(e for e in ea if e["day"] == 7)
    lines = [
        f"In silico, {meta['n']} virtual patients, engine `{meta['engine']}` vs truth `{meta['truth']}` "
        "(see [`SIMULATION.md`](SIMULATION.md) for methods and caveats):",
        "",
        "| Question | Status quo | Model-informed |",
        "|---|---|---|",
        f"| Blood tests during titration (mean) | {u['tests_mean']:.1f} stepwise · {p['tests_mean']:.1f} proportional | **{m['tests_mean']:.1f}** |",
        f"| Stable dose confirmed within 8 weeks | {u['confirmed_pct']:.0f}% · {p['confirmed_pct']:.0f}% | **{m['confirmed_pct']:.0f}%** |",
        f"| Final dose truly in 0.6–0.8 / 0.5–0.9 | {u['in_range_pct']:.0f}/{u['in_wide_pct']:.0f}% · {p['in_range_pct']:.0f}/{p['in_wide_pct']:.0f}% | **{m['in_range_pct']:.0f}/{m['in_wide_pct']:.0f}%** |",
        f"| Ever on a dose with true level > 1.0 | {u['over1_pct']:.0f}% · {p['over1_pct']:.0f}% | **{m['over1_pct']:.0f}%** |",
        f"| Mistimed level: median error vs the correctly timed level | {allm['naive']:.3f} mmol/L | **{allm['engine']:.3f}** |",
        f"| Mistimed level: wrong below/in/above-range call | {allm['naive_wrong']:.0f}% | **{allm['engine_wrong']:.0f}%** (90% interval coverage {allm['coverage']:.0f}%) |",
        f"| Thiazide started: true level > 1.0 / > 1.2 | {th['over1'][0]:.0f}% / {th['over12'][0]:.0f}% (no change) | **{th['over1'][1]:.0f}% / {th['over12'][1]:.0f}%** (engine-adjusted) |",
        f"| Level 1 day after a dose change: median error | {d1['naive']:.3f} | **{d1['with_history']:.3f}** (with one earlier level) |",
        f"| Level 7 days after a dose change: median error | {d7['naive']:.3f} | **{d7['with_history']:.3f}** |",
    ]
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    meta, rows = load()
    if a.json:
        Path(a.json).write_text(json.dumps(web(meta, rows)))
    print(headline_md(meta, rows))


if __name__ == "__main__":
    main()
