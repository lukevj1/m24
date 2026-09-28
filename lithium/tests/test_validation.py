"""The validation harness end to end, on synthetic data only (no real patient data)."""

import csv
import gzip
from pathlib import Path

import numpy as np
import pytest

from validation import evaluate, schema, synthetic
from validation.mimic_iv import convert, dose_mmol, formulation


def test_synthetic_dataset_round_trip_and_forecasting(tmp_path):
    synthetic.make_dataset(tmp_path, n=6, seed=3, days=200)
    records = schema.load(tmp_path)
    assert len(records) == 6
    r = records[0]
    assert r.levels and r.doses and r.creatinine
    assert all(0 <= lv.time for lv in r.levels)
    rows = evaluate.run(records[:3], n_draws=60, procs=1)
    ok = [x for x in rows if "skip" not in x]
    assert ok, "no evaluable levels"
    follow = [x for x in ok if x["n_prior"] >= 1]
    assert follow and all(x["engine"]["lo"] <= x["engine"]["hi"] for x in follow)
    assert all(np.isfinite(x["engine"]["med"]) for x in follow)
    res = evaluate.summarise(rows, boot=20)
    assert res["follow_up"]["engine"]["n"] == len([x for x in follow if x.get("proportional") is not None])
    md = evaluate.report_md(res)
    assert "Forecasting the next level" in md


def test_prospective_covariates_do_not_peek_at_future_creatinine(tmp_path):
    synthetic.make_dataset(tmp_path, n=1, seed=5, days=120)
    rec = schema.load(tmp_path)[0]
    t_prev = rec.levels[1].time
    tl, ok = rec.covariates(t_prev + 1.0)
    assert ok
    assert all(t <= t_prev + 1.0 for t, _ in tl.creatinine)


def _gz(path: Path, header: list[str], rows: list[list]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def _fake_mimic(root: Path):
    h = root / "hosp"
    _gz(h / "d_labitems.csv.gz", ["itemid", "label", "fluid", "category"],
        [[1, "Lithium", "Blood", "Chemistry"], [2, "Creatinine", "Blood", "Chemistry"],
         [3, "Creatinine", "Urine", "Chemistry"]])
    _gz(h / "patients.csv.gz", ["subject_id", "gender", "anchor_age", "anchor_year", "anchor_year_group", "dod"],
        [[10, "F", 50, 2150, "2017 - 2019", ""]])
    _gz(h / "admissions.csv.gz", ["subject_id", "hadm_id", "admittime", "dischtime"],
        [[10, 100, "2152-03-01 10:00:00", "2152-03-12 12:00:00"]])
    emar = [[10, 100, f"e{d}", f"2152-03-{d:02d} 21:00:00", "Lithium Carbonate", "Administered"] for d in range(1, 11)]
    emar.append([10, 100, "e99", "2152-03-11 21:00:00", "Lithium Carbonate", "Not Given"])
    _gz(h / "emar.csv.gz", ["subject_id", "hadm_id", "emar_id", "charttime", "medication", "event_txt"], emar)
    det = [[10, f"e{d}", "", "600", "mg", "Lithium Carbonate SR (Lithobid) 300mg Tab"] for d in range(1, 11)]
    _gz(h / "emar_detail.csv.gz", ["subject_id", "emar_id", "parent_field_ordinal", "dose_given", "dose_given_unit",
                                   "product_description"], det)
    _gz(h / "labevents.csv.gz", ["subject_id", "hadm_id", "itemid", "charttime", "valuenum", "valueuom"],
        [[10, 100, 1, "2152-03-05 09:00:00", "0.62", "mmol/L"], [10, 100, 1, "2152-03-10 09:00:00", "0.71", "mmol/L"],
         [10, 100, 2, "2152-03-01 12:00:00", "1.0", "mg/dL"], [10, "", 2, "2152-03-09 08:00:00", "0.9", "mg/dL"],
         [10, 100, 3, "2152-03-05 09:00:00", "80", "mg/dL"]])
    _gz(h / "omr.csv.gz", ["subject_id", "chartdate", "seq_num", "result_name", "result_value"],
        [[10, "2152-02-01", 1, "Weight (Lbs)", "176"], [10, "2152-02-01", 1, "Height (Inches)", "65"]])


def test_mimic_converter_on_a_fake_extract(tmp_path):
    _fake_mimic(tmp_path / "mimic")
    counts = convert(tmp_path / "mimic", tmp_path / "out")
    assert counts["admissions written"] == 1 and counts["levels written"] == 2
    rec = schema.load(tmp_path / "out")[0]
    assert rec.sex == "F" and rec.age_years == pytest.approx(52.0)
    assert rec.weight_kg == pytest.approx(176 * 0.45359237, rel=1e-3)
    assert rec.height_cm == pytest.approx(165.1, rel=1e-3)
    given = [d for d, s in zip(rec.doses, rec.dose_source) if s == "administered"]
    assert len(given) == 10                                  # "Not Given" excluded
    assert given[0].amount == pytest.approx(600 * 2 / 73.891, rel=1e-3)
    assert given[0].formulation == "SR"
    assert any(s == "assumed" for s in rec.dose_source)      # first dose within 36 h of admission
    assert [v for _, v in rec.creatinine] == [pytest.approx(88.4, abs=0.1), pytest.approx(79.6, abs=0.1)]
    assert [lv.value for lv in rec.levels] == [0.62, 0.71]


def test_formulation_and_dose_parsing():
    assert formulation("Lithium Carbonate 300mg Cap") == "IR"
    assert formulation("Lithium Carbonate SR (Lithobid) 300mg Tab") == "SR"
    assert formulation("Lithium Carbonate ER 450 mg") == "SR"
    assert formulation("Lithium Citrate 8mEq/5mL Syrup") == "LIQ"
    assert formulation("crushed lithium carbonate") == "IR"
    assert dose_mmol("8", "mEq", "Lithium Citrate")[0] == pytest.approx(8.0)
    assert dose_mmol("5", "mL", "Lithium Citrate 8mEq/5mL")[0] == pytest.approx(8.0)
    assert dose_mmol("2", "TAB", "Lithium Carbonate 300mg Tab")[0] == pytest.approx(600 * 2 / 73.891)
    assert dose_mmol("", "mg", "x")[0] is None
