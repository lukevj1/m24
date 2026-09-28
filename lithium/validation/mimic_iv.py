"""Convert MIMIC-IV (hosp module) to the canonical validation format.

MIMIC-IV is credentialed data (PhysioNet). Run this where your data use
agreement allows the data to be. Nothing here downloads or uploads anything,
and it prints only aggregate counts. Read docs/VALIDATION.md first.

What it extracts
----------------
* Lithium administrations from the eMAR (emar + emar_detail): documented
  administration times, doses (mg of lithium carbonate, mEq for citrate
  solution) and the product, from which the release type is inferred.
* Serum lithium and creatinine from labevents, selected by label in
  d_labitems (no hard-coded item ids). Creatinine mg/dL -> umol/L (x 88.42).
* Sex and age (anchor_age, shifted by calendar year), weight and height from omr.

One record per hospital admission with at least one lithium level and at least
one documented administration. Doses taken at home before admission are not
recorded. If the first documented administration is within 36 h of admission,
the first documented day's regimen is assumed for the 14 days before it
(source = "assumed"); the evaluator reports levels influenced by assumed doses
separately. Admissions with any level above --exclude-above (default 2.5 mmol/L,
suggesting poisoning rather than maintenance treatment) are excluded.

    python -m validation.mimic_iv /path/to/mimic-iv/3.1 OUT_DIR
    python -m validation.evaluate OUT_DIR --out report.md --json metrics.json
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import Counter, defaultdict
from datetime import datetime, timedelta
from pathlib import Path

from .schema import parse_time, write

MG_DL_TO_UMOL_L = 88.42
LB_TO_KG = 0.45359237
IN_TO_CM = 2.54
CITRATE_MEQ_PER_ML = 8.0 / 5.0          # US lithium citrate oral solution: 8 mEq per 5 mL
_SR = re.compile(r"\b(er|sr|cr|xr|extended|sustained|controlled|slow|lithobid|eskalith[- ]?cr)\b")
_LIQ = re.compile(r"citrate|syrup|solution|soln|liquid|meq\s*/\s*5\s*ml")
_MG = re.compile(r"(\d+(?:\.\d+)?)\s*mg")


def _open(root: Path, name: str):
    for sub in ("hosp", ""):
        for ext in (".csv.gz", ".csv"):
            p = root / sub / f"{name}{ext}" if sub else root / f"{name}{ext}"
            if p.exists():
                return gzip.open(p, "rt", newline="") if ext == ".csv.gz" else p.open(newline="")
    raise FileNotFoundError(f"{name}.csv(.gz) not found under {root} or {root / 'hosp'}")


def _rows(root: Path, name: str, need: list[str]):
    with _open(root, name) as fh:
        rd = csv.DictReader(fh)
        missing = [c for c in need if c not in (rd.fieldnames or [])]
        if missing:
            raise ValueError(f"{name}: missing column(s) {missing}; found {rd.fieldnames}")
        yield from rd


def _float(s) -> float | None:
    try:
        return float(str(s).strip())
    except (TypeError, ValueError):
        return None


def formulation(text: str) -> str:
    s = (text or "").lower()
    if _LIQ.search(s):
        return "LIQ"
    if _SR.search(s):
        return "SR"
    return "IR"


def dose_mmol(dose_given, unit, product: str) -> tuple[float | None, str]:
    """mmol Li+ for one documented administration row, and a note on how it was derived."""
    v = _float(dose_given)
    if v is None or v <= 0:
        return None, "no numeric dose"
    u = (unit or "").strip().lower()
    prod = (product or "").lower()
    if u in ("meq", "mmol"):
        return v, "mEq"
    if u == "ml":
        return v * CITRATE_MEQ_PER_ML, "mL of 8 mEq/5 mL solution"
    if u == "mg":
        if formulation(prod) == "LIQ":
            # US labelling states solution doses as mg carbonate-equivalent (5 mL = 300 mg = 8 mEq)
            return v * 8.0 / 300.0, "mg carbonate-equivalent (solution)"
        return v * 2.0 / 73.891, "mg lithium carbonate"
    if u in ("tab", "cap", "tablet", "capsule"):
        m = _MG.search(prod)
        if m:
            return v * float(m.group(1)) * 2.0 / 73.891, "tablets x strength"
    return None, f"unit {unit!r} not understood"


def convert(root: str | Path, out_dir: str | Path, *, exclude_above: float = 2.5,
            level_time_quality: str = "approximate", assume_home_days: float = 14.0) -> Counter:
    root = Path(root)
    n = Counter()

    # 1. Lab item ids by label (no hard-coded ids).
    li_items, cr_items = set(), set()
    for r in _rows(root, "d_labitems", ["itemid", "label", "fluid"]):
        label, fluid = r["label"].strip().lower(), (r.get("fluid") or "").strip().lower()
        if label == "lithium" and fluid == "blood":
            li_items.add(r["itemid"])
        if label == "creatinine" and fluid == "blood":
            cr_items.add(r["itemid"])
    if not li_items:
        raise ValueError("no blood lithium item found in d_labitems")

    # 2. eMAR lithium administrations.
    admin: dict[str, dict] = {}
    events = Counter()
    for r in _rows(root, "emar", ["subject_id", "hadm_id", "emar_id", "charttime", "medication", "event_txt"]):
        if "lithium" not in (r["medication"] or "").lower():
            continue
        ev = (r["event_txt"] or "").strip()
        events[ev] += 1
        e = ev.lower()
        if "administered" not in e or e.startswith("not"):
            continue
        if not r["hadm_id"]:
            n["administrations without an admission id"] += 1
            continue
        admin[r["emar_id"]] = {"subject": r["subject_id"], "hadm": r["hadm_id"], "time": parse_time(r["charttime"]),
                               "medication": r["medication"]}
    n["eMAR lithium events"] = sum(events.values())
    n["eMAR lithium administrations"] = len(admin)

    detail: dict[str, list[dict]] = defaultdict(list)
    for r in _rows(root, "emar_detail", ["emar_id", "dose_given", "dose_given_unit", "product_description"]):
        if r["emar_id"] in admin:
            detail[r["emar_id"]].append(r)

    doses_by_adm: dict[str, list[dict]] = defaultdict(list)
    for eid, a in admin.items():
        rows = detail.get(eid, [])
        top = [r for r in rows if not (r.get("parent_field_ordinal") or "").strip()]
        use = top or rows
        vals = [dose_mmol(r["dose_given"], r["dose_given_unit"], r.get("product_description") or a["medication"])
                for r in use]
        good = [v for v, _ in vals if v is not None]
        if not good:
            n["administrations without a usable dose"] += 1
            continue
        mmol = good[0] if top else sum(good)
        if top and len({round(v, 3) for v in good}) > 1:
            n["administrations with conflicting doses (skipped)"] += 1
            continue
        text = " ".join([a["medication"]] + [r.get("product_description") or "" for r in rows])
        doses_by_adm[a["hadm"]].append({"subject": a["subject"], "time": a["time"], "mmol": mmol,
                                        "formulation": formulation(text)})

    subjects = {d["subject"] for ds in doses_by_adm.values() for d in ds}

    # 3. Admissions, patients, labs, body size for those subjects.
    adm_info = {}
    for r in _rows(root, "admissions", ["subject_id", "hadm_id", "admittime", "dischtime"]):
        if r["hadm_id"] in doses_by_adm:
            adm_info[r["hadm_id"]] = {"subject": r["subject_id"], "admit": parse_time(r["admittime"]),
                                      "disch": parse_time(r["dischtime"])}
    pts = {}
    for r in _rows(root, "patients", ["subject_id", "gender", "anchor_age", "anchor_year"]):
        if r["subject_id"] in subjects:
            pts[r["subject_id"]] = r
    li_by_subject, cr_by_subject = defaultdict(list), defaultdict(list)
    for r in _rows(root, "labevents", ["subject_id", "hadm_id", "itemid", "charttime", "valuenum", "valueuom"]):
        if r["subject_id"] not in subjects:
            continue
        v = _float(r["valuenum"])
        if v is None:
            continue
        if r["itemid"] in li_items:
            li_by_subject[r["subject_id"]].append((parse_time(r["charttime"]), v, r["hadm_id"]))
        elif r["itemid"] in cr_items:
            uom = (r["valueuom"] or "").lower()
            cr_by_subject[r["subject_id"]].append((parse_time(r["charttime"]),
                                                   v if "umol" in uom else v * MG_DL_TO_UMOL_L))
    body = defaultdict(lambda: {"kg": [], "cm": []})
    try:
        for r in _rows(root, "omr", ["subject_id", "chartdate", "result_name", "result_value"]):
            if r["subject_id"] not in subjects:
                continue
            name, v = r["result_name"].lower(), _float(r["result_value"])
            if v is None:
                continue
            d = parse_time(r["chartdate"])
            if name.startswith("weight"):
                body[r["subject_id"]]["kg"].append((d, v * LB_TO_KG if "lb" in name else v))
            elif name.startswith("height"):
                body[r["subject_id"]]["cm"].append((d, v * IN_TO_CM if "inch" in name else v))
    except FileNotFoundError:
        n["omr table not found: body size imputed"] += 1

    # 4. Assemble one record per admission.
    tables: dict[str, list[dict]] = {k: [] for k in ("patients", "doses", "levels", "creatinine")}
    fmt = lambda t: t.strftime("%Y-%m-%d %H:%M:%S")  # noqa: E731
    for hadm, doses in sorted(doses_by_adm.items()):
        info = adm_info.get(hadm)
        if info is None:
            n["admissions not found in admissions table"] += 1
            continue
        subj = info["subject"]
        lv = [(t, v) for t, v, h in li_by_subject.get(subj, [])
              if h == hadm or (not h and info["admit"] <= t <= info["disch"])]
        if not lv:
            n["admissions with administrations but no level"] += 1
            continue
        if max(v for _, v in lv) > exclude_above:
            n[f"admissions excluded: a level above {exclude_above:g}"] += 1
            continue
        p = pts.get(subj)
        if p is None:
            n["subject missing from patients table"] += 1
            continue
        pid = f"{subj}-{hadm}"
        age = float(p["anchor_age"]) + (info["admit"].year - int(p["anchor_year"]))
        near = lambda xs: [v for d, v in xs if abs((d - info["admit"]).days) <= 365]  # noqa: E731
        kg, cm = near(body[subj]["kg"]), near(body[subj]["cm"]) or [v for _, v in body[subj]["cm"]]
        med = lambda xs: sorted(xs)[len(xs) // 2] if xs else None  # noqa: E731
        tables["patients"].append({"patient_id": pid, "sex": "F" if p["gender"].upper().startswith("F") else "M",
                                   "age_years": age, "height_cm": med(cm), "weight_kg": med(kg)})
        doses.sort(key=lambda d: d["time"])
        first = doses[0]["time"]
        if (first - info["admit"]).total_seconds() / 3600.0 <= 36.0 and assume_home_days > 0:
            day1 = [d for d in doses if (d["time"] - first).total_seconds() < 24 * 3600]
            for k in range(int(assume_home_days), 0, -1):
                for d in day1:
                    t = d["time"] - timedelta(days=k)
                    tables["doses"].append({"patient_id": pid, "time": fmt(t), "mmol": round(d["mmol"], 4),
                                            "formulation": d["formulation"], "source": "assumed"})
            n["admissions with an assumed pre-admission regimen"] += 1
        for d in doses:
            tables["doses"].append({"patient_id": pid, "time": fmt(d["time"]), "mmol": round(d["mmol"], 4),
                                    "formulation": d["formulation"], "source": "administered"})
        for t, v in sorted(lv):
            tables["levels"].append({"patient_id": pid, "time": fmt(t), "value_mmol_l": v,
                                     "time_quality": level_time_quality})
        lo, hi = info["admit"] - timedelta(days=7), info["disch"]
        for t, v in sorted(cr_by_subject.get(subj, [])):
            if lo <= t <= hi:
                tables["creatinine"].append({"patient_id": pid, "time": fmt(t), "umol_l": round(v, 1)})
        n["admissions written"] += 1
        n["levels written"] += len(lv)
    write(out_dir, tables)
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("mimic_root", help="MIMIC-IV version directory (containing hosp/)")
    ap.add_argument("out", help="output directory for the canonical CSVs (keep it inside your DUA environment)")
    ap.add_argument("--exclude-above", type=float, default=2.5)
    ap.add_argument("--level-time-quality", choices=["exact", "approximate", "unknown"], default="approximate")
    ap.add_argument("--assume-home-days", type=float, default=14.0,
                    help="days of the first documented regimen assumed before admission (0 to disable)")
    a = ap.parse_args()
    counts = convert(a.mimic_root, a.out, exclude_above=a.exclude_above, level_time_quality=a.level_time_quality,
                     assume_home_days=a.assume_home_days)
    print("Conversion summary (aggregate counts only):")
    for k, v in counts.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
