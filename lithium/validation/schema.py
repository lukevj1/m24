"""Canonical individual-level dataset format for validation.

A dataset is a directory of CSV files (UTF-8, header row, one row per event).
Times are ISO 8601 date-times ("2150-03-04 08:30" or "2150-03-04T08:30:00").
De-identified, date-shifted data are fine: only intervals and clock times are
used, and the clock time must be real local time (not shifted).

    patients.csv    patient_id, sex (F/M), age_years (at the first event),
                    height_cm (optional), weight_kg (optional)
    doses.csv       patient_id, time, mg, salt (carbonate|citrate),
                    formulation (IR|SR|LIQ), source, mmol (optional)
                    source: administered | prescribed | reported | assumed
                    mmol: mmol of Li+ (= mEq), overrides mg/salt when given
    levels.csv      patient_id, time, value_mmol_l, time_quality
                    time_quality: exact | approximate | unknown
    creatinine.csv  patient_id, time, umol_l
    weights.csv     optional: patient_id, time, kg

``load`` turns each patient into a :class:`Record`. The hours axis starts at
midnight of the patient's first event, so ``t % 24`` is still the clock time.
Converters (e.g. validation/mimic_iv.py) write this format with ``write``.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

from lithos.case import CovariateTimeline, Level
from lithos.pk import Dose
from lithos.units import mg_to_mmol

COLUMNS = {
    "patients": ["patient_id", "sex", "age_years", "height_cm", "weight_kg"],
    "doses": ["patient_id", "time", "mg", "salt", "formulation", "source", "mmol"],
    "levels": ["patient_id", "time", "value_mmol_l", "time_quality"],
    "creatinine": ["patient_id", "time", "umol_l"],
    "weights": ["patient_id", "time", "kg"],
}
REQUIRED = {
    "patients": ["patient_id", "sex", "age_years"],
    "doses": ["patient_id", "time", "formulation"],
    "levels": ["patient_id", "time", "value_mmol_l"],
    "creatinine": ["patient_id", "time", "umol_l"],
}
TIMING_SD_H = {"exact": 0.25, "approximate": 1.0, "unknown": 2.0}
DEFAULT_HEIGHT_CM = {"F": 163.0, "M": 177.0}
DEFAULT_WEIGHT_KG = {"F": 75.0, "M": 88.0}   # used only when no weight is recorded; flagged


def parse_time(s: str) -> datetime:
    s = s.strip().replace("T", " ")
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d %H:%M:%S.%f", "%Y-%m-%d"):
        try:
            return datetime.strptime(s, fmt)
        except ValueError:
            continue
    raise ValueError(f"unrecognised time {s!r}; use ISO 8601, e.g. 2150-03-04 08:30")


def _num(s: str | None) -> float | None:
    if s is None or str(s).strip() == "":
        return None
    return float(s)


def _read(path: Path, table: str, required: bool = True) -> list[dict]:
    if not path.exists():
        if required:
            raise FileNotFoundError(f"{path} is missing (expected columns: {', '.join(COLUMNS[table])})")
        return []
    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    have = set(rows[0].keys()) if rows else set()
    missing = [c for c in REQUIRED.get(table, []) if rows and c not in have]
    if missing:
        raise ValueError(f"{path.name}: missing column(s) {missing}; found {sorted(have)}")
    return rows


@dataclass
class Record:
    """Everything recorded about one patient (or one hospital admission)."""

    patient_id: str
    sex: str
    age_years: float
    height_cm: float
    weight_kg: float
    body_size_imputed: bool
    t0: datetime
    doses: list[Dose]
    dose_source: list[str]
    levels: list[Level]
    creatinine: list[tuple[float, float]]
    weights: list[tuple[float, float]] = field(default_factory=list)

    def covariates(self, creatinine_until: float | None = None) -> tuple[CovariateTimeline, bool]:
        """Covariate timeline using creatinine results up to ``creatinine_until``
        (all of them if None). Returns (timeline, True if any creatinine was
        available by then; if not, the earliest result is used and flagged)."""
        scr = sorted(self.creatinine)
        if creatinine_until is not None:
            avail = [c for c in scr if c[0] <= creatinine_until]
            ok = bool(avail)
            scr = avail or scr[:1]
        else:
            ok = bool(scr)
        if not scr:
            scr = [(0.0, 80.0)]
        tl = CovariateTimeline(sex=self.sex, age_at_t0=self.age_years, weight=self.weight_kg,
                               height=self.height_cm, creatinine=scr, weights=sorted(self.weights))
        return tl, ok


def load(directory: str | Path) -> list[Record]:
    d = Path(directory)
    patients = _read(d / "patients.csv", "patients")
    doses = _read(d / "doses.csv", "doses")
    levels = _read(d / "levels.csv", "levels")
    creat = _read(d / "creatinine.csv", "creatinine")
    weights = _read(d / "weights.csv", "weights", required=False)

    by: dict[str, dict[str, list]] = {}
    for table, rows in (("doses", doses), ("levels", levels), ("creatinine", creat), ("weights", weights)):
        for r in rows:
            by.setdefault(r["patient_id"], {}).setdefault(table, []).append(r)

    out: list[Record] = []
    for p in patients:
        pid = p["patient_id"]
        ev = by.get(pid, {})
        stamps = [parse_time(r["time"]) for t in ("doses", "levels", "creatinine") for r in ev.get(t, [])]
        if not stamps or not ev.get("levels"):
            continue
        first = min(stamps)
        t0 = datetime(first.year, first.month, first.day)

        def hours(s: str) -> float:
            return (parse_time(s) - t0).total_seconds() / 3600.0

        sex = p["sex"].strip().upper()[:1]
        if sex not in ("F", "M"):
            raise ValueError(f"patient {pid}: sex must be F or M (got {p['sex']!r})")
        h = _num(p.get("height_cm"))
        w = _num(p.get("weight_kg"))
        wts = [(hours(r["time"]), float(r["kg"])) for r in ev.get("weights", []) if _num(r.get("kg"))]
        imputed = h is None or (w is None and not wts)
        if w is None and wts:
            w = sorted(v for _, v in wts)[len(wts) // 2]

        dose_list, sources = [], []
        for r in sorted(ev.get("doses", []), key=lambda r: parse_time(r["time"])):
            mmol = _num(r.get("mmol"))
            if mmol is None:
                mg = _num(r.get("mg"))
                if mg is None:
                    raise ValueError(f"patient {pid}: a dose needs mg (with salt) or mmol")
                mmol = mg_to_mmol(mg, (r.get("salt") or "carbonate").strip().lower())
            form = r["formulation"].strip().upper()
            if form not in ("IR", "SR", "LIQ"):
                raise ValueError(f"patient {pid}: formulation must be IR, SR or LIQ (got {form!r})")
            dose_list.append(Dose(hours(r["time"]), mmol, form))
            sources.append((r.get("source") or "administered").strip().lower())

        lvls = []
        for r in sorted(ev["levels"], key=lambda r: parse_time(r["time"])):
            v = _num(r["value_mmol_l"])
            if v is None or v <= 0:
                continue
            q = (r.get("time_quality") or "exact").strip().lower()
            lvls.append(Level(hours(r["time"]), v, TIMING_SD_H.get(q, TIMING_SD_H["unknown"]), note=q))

        scr = [(hours(r["time"]), float(r["umol_l"])) for r in ev.get("creatinine", []) if _num(r.get("umol_l"))]
        out.append(Record(pid, sex, float(p["age_years"]), h if h is not None else DEFAULT_HEIGHT_CM[sex],
                          w if w is not None else DEFAULT_WEIGHT_KG[sex], imputed, t0, dose_list, sources,
                          lvls, scr, wts))
    return out


def write(directory: str | Path, tables: dict[str, list[dict]]) -> None:
    """Write converter output in the canonical format (only known columns)."""
    d = Path(directory)
    d.mkdir(parents=True, exist_ok=True)
    for name, rows in tables.items():
        cols = COLUMNS[name]
        with (d / f"{name}.csv").open("w", newline="", encoding="utf-8") as fh:
            wr = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            wr.writeheader()
            for r in rows:
                wr.writerow({c: ("" if r.get(c) is None else r.get(c)) for c in cols})
