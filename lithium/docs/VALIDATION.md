# Testing the engine against real patients

**Short answer.**
- **Individual-level data.** No real, individual-level lithium dataset could be reached from the environment this prototype was built in. The useful datasets are credentialed or held by research groups.
- **Two things were done instead:**
  1. The engine's population prior was checked against **eight published cohorts of real patients** (section 1).
  2. A **ready-to-run validation harness** was built for individual-level data (section 2). It has a converter for MIMIC-IV (a credentialed ICU and hospital database) and has been tested end to end on synthetic data. Anyone with approved access to a dataset can run it where the data live, and share only the aggregate report.

## 1. Published real-patient cohorts (population level)

For each cohort, a patient resembling it was given to the engine with **no lithium levels**. The engine's prior predictions were then compared with what the paper measured in real patients. Reproduce with `python -m validation.literature`.

| Cohort | Quantity | Published (real patients) | Band used | Engine prior: median (90% range) | Engine median vs published mean | Engine median in band? |
|---|---|---|---|---|---|---|
| Older adults (PMID 3110219) | Clearance | 15.6 +/- 4.0 mL/min | 11.6-19.6 | 13.2 (8.3-20.3) | -16% | inside |
| Older adults (PMID 3110219) | Terminal half-life | 26.9 +/- 5.5 h | 21.4-32.4 | 33.5 (20.2-56.4) | +25% | 3% above |
| Older adults (PMID 3110219) | V_area | 0.64 +/- 0.16 L/kg | 0.48-0.8 | 0.56 (0.40-0.77) | -12% | inside |
| Older adults (PMID 3110219) | Steady-state peak on 0.21 mmol/kg/day | 0.83 +/- 0.25 mmol/L | 0.58-1.08 | 1.09 (0.78-1.52) | +31% | 1% above |
| Older adults (PMID 3110219) | Time to peak | 2.2 +/- 1.2 h | 1-3.4 | 1.9 (1.8-2.1) | -14% | inside |
| Older adults (PMID 3110219) | Distribution phase ~complete (kidney model only) | 10.6 +/- 3.0 h | 7.6-13.6 | 10.0 (6.9-13.6) | -6% | inside |
| Obese adults (PMID 8162665) | Vss | 0.42 +/- 0.09 L/kg | 0.33-0.51 | 0.50 (0.36-0.70) | +19% | inside |
| Obese adults (PMID 8162665) | Clearance | 33.9 +/- 7.0 mL/min | 26.9-40.9 | 31.5 (20.2-48.5) | -7% | inside |
| Normal-weight controls (PMID 8162665) | Vss | 0.66 +/- 0.16 L/kg | 0.5-0.82 | 0.69 (0.50-0.96) | +5% | inside |
| Normal-weight controls (PMID 8162665) | Clearance | 23.0 +/- 6.2 mL/min | 16.8-29.2 | 25.3 (16.5-39.9) | +10% | inside |
| Healthy adults (PMID 3118402) | Clearance | 22-43 (range) mL/min | 22-43 | 26.4 (17.1-41.3) | -19% | inside |
| Inpatients (PMID 657687) | Volume | 0.79 +/- 0.34 L/kg | 0.45-1.13 | 0.69 (0.49-0.97) | -12% | inside |
| Inpatients (PMID 657687) | Renal clearance | 24.4 +/- 8.0 mL/min | 16.4-32.4 | 25.0 (16.4-38.6) | +3% | inside |
| Inpatients (PMID 657687) | Half-life | 28.9 +/- 7.9 h | 21-36.8 | 23.9 (14.1-40.6) | -17% | inside |
| Single-dose study (PMID 2079355) | Clearance | 33.2 +/- 15.5 mL/min | 17.7-48.7 | 25.2 (16.5-39.0) | -24% | inside |
| Single-dose study (PMID 2079355) | Half-life | 15.3 +/- 6.1 h | 9.2-21.4 | 19.8 (11.0-34.5) | +29% | inside |
| IR vs CR (PMID 3089949) | Time to peak, IR | 1.44 h | 1-2 | 1.8 (1.6-1.9) | +22% | inside |
| IR vs CR (PMID 3089949) | Time to peak, CR | 4.25 h | 3.5-5 | 4.2 (3.7-4.7) | -1% | inside |
| IR vs CR (PMID 3089949) | Cmax ratio CR/IR | 0.78 ratio | 0.7-0.86 | 0.71 (0.68-0.73) | -10% | inside |
| Once-daily SR (Reddy 2014) | 12-h / 24-h level | 1.37 (0.82 vs 0.60) ratio | 1.25-1.5 | 1.53 (1.31-2.04) | +12% | 2% above |
| Healthy volunteers (Stip 2001) | 12-h level on 1569 mg/day | 0.67-0.74 (weekly means) mmol/L | 0.6-0.8 | 0.85 (0.49-1.45) | +20% | 6% above |

Engine median inside the published band: 17/21. Within 20% of the published mean: 15/21. Band = mean +/- 1 SD where an SD was reported, the reported range, or a stated tolerance (time to peak, ratios, Stip) where it was not.

Matched patients (assumptions where the paper did not report covariates):

- Older adults (PMID 3110219): 72 y woman, 68 kg, creatinine 85 umol/L (cohort covariates not reported)
- Obese adults (PMID 8162665): 35 y man, 120 kg, 175 cm, creatinine 85 (cohort covariates not reported)
- Normal-weight controls (PMID 8162665): 35 y man, 72 kg, 176 cm, creatinine 85 (cohort covariates not reported)
- Healthy adults (PMID 3118402): 30 y man, 72 kg, creatinine 85
- Inpatients (PMID 657687): 35 y man, 72 kg, creatinine 85 (cohort covariates not reported)
- Single-dose study (PMID 2079355): 35 y woman, 70 kg, creatinine 70 (cohort covariates not reported)
- IR vs CR (PMID 3089949): 35 y man, 72 kg; single dose
- Once-daily SR (Reddy 2014): 40 y woman, 70 kg; SR at 21:00
- Healthy volunteers (Stip 2001): 25 y man, 72 kg, creatinine 85; twice daily (regimen assumed)

**What this shows.**
- The prior is in the right place for most quantities:
  - clearance in normal-weight, obese and healthy adults;
  - volume in normal-weight and obese adults;
  - absorption of immediate- and controlled-release tablets;
  - timing of the distribution phase in older adults.
- **Older adults.** The prior is "slower" than the one older cohort: half-life +25%, peak +31% on the same dose per kg.
  - The model would therefore start older patients on lower doses than needed. That is the safer direction of error, and the first level corrects it.
  - This matches what the simulations showed: the model titrates older patients more slowly.
- **Single-dose cohort.** The prior is slower than this cohort too (clearance −24%, half-life +29%). One possible reason: single-dose designs can under-sample the slow terminal phase.
- **Once-daily modified-release.** The 12-hour to 24-hour ratio is 12% higher than reported (1.53 vs 1.37). The modified-release absorption model may be too fast. For once-daily modified-release this matters, because it sets how much a mistimed level is corrected.

**What it does not show.**
- **Not independent.** Some of these cohorts informed the calibration: volume in obesity and old age, the shape of the peak, and modified-release absorption. This is a consistency check, not validation.
- **Lenient criterion.** The band (mean ± 1 SD) describes the spread between patients, not uncertainty in the mean. The stricter column, within 20% of the published mean, is met for 15 of 21 quantities.
- **Assumed covariates.** Most papers did not report their patients' covariates, so the matched patients are assumptions (listed under the table).
- **Nothing about individuals.** Population agreement says nothing about how well the engine forecasts an individual's next level. Standardising a mistimed level, and the drift model, are also untested by this check. That needs section 2.

Context for why individual-level validation matters: an external evaluation of 10 published lithium population models found median absolute prediction errors of 24–138% [E39].

## 2. Individual-level validation harness

### Design: iterative forecasting

Each patient's levels are taken in time order. Level *k* is forecast using only what was known before it was drawn:
- the doses recorded up to the draw;
- all earlier levels;
- by default, only creatinine results available by the time of the previous level (`--covariates prospective`).

The forecast is then compared with the measured value. This mirrors how the tool would be used. It follows the design used to validate model-informed precision dosing for vancomycin (e.g. Broeker et al., *Clin Microbiol Infect* 2019).

**Comparators receive the same information.**

| Method | Forecast of level *k* |
|---|---|
| Engine | Bayesian forecast from earlier levels (the tool) |
| Population prior | The engine with no levels: what a covariate-only dosing calculator would say |
| Previous level, unchanged | "It will be what it was": what a clinician implicitly assumes on the same dose |
| Previous level × dose ratio | Linear scaling by the change in average daily dose: the usual mental arithmetic |

**Metrics**, overall and by stratum (95% CIs from a patient-level bootstrap):
- median prediction error (MDPE, bias) and median absolute prediction error (MDAPE, imprecision): the metrics of the published external evaluations [E39]. Medians are used because a few very low levels, for example after a dose missed the night before the test, make means of relative errors meaningless;
- mean absolute error;
- share within ±20% and within ±0.1 mmol/L;
- coverage of the model's 90% prediction interval;
- whether levels above 1.0 were anticipated.

**Strata:**
- hours since the last dose;
- number of earlier levels;
- age;
- eGFR;
- formulation;
- whether the dose changed;
- days since the previous level;
- how the sample time was recorded.

### Data format

A directory of CSV files (`patients`, `doses`, `levels`, `creatinine`, optional `weights`), specified in `validation/schema.py`. Converters write this format; the evaluator reads it.
- Doses carry a source: `administered`, `prescribed`, `reported` or `assumed`.
- Levels carry how well their time is known: `exact`, `approximate` or `unknown`.

### Running it

```bash
# MIMIC-IV (hosp module), inside your approved environment
python -m validation.mimic_iv /path/to/mimic-iv/3.1 /secure/lithium_canonical
python -m validation.evaluate /secure/lithium_canonical --out report.md --json metrics.json

# Any other dataset: write the canonical CSVs, then run the same evaluator
python -m validation.evaluate /path/to/canonical --procs 8

# Synthetic data, to see the report format
python -m validation.synthetic /tmp/li_synth --n 100
python -m validation.evaluate /tmp/li_synth
```

**What the MIMIC-IV converter does:**
- takes documented administration times and doses from the eMAR, and infers the release type from the product;
- selects serum lithium and creatinine by label (no hard-coded item ids);
- converts creatinine mg/dL to µmol/L;
- takes age from `anchor_age`, and weight and height from `omr`;
- writes one record per admission.

Doses taken at home before admission are not recorded:
- If the first administration is within 36 h of admission, the first documented day's regimen is assumed for the 14 days before it.
- Levels influenced by assumed doses are reported separately, not in the primary analysis.
- Admissions with any level above 2.5 mmol/L (probable poisoning) are excluded.

All thresholds are command-line options. The converter prints only aggregate counts.

### Demonstration on synthetic data (not validation)

`python -m validation.synthetic` makes 100 virtual outpatients followed for a year. Each is titrated the usual way, then has three-monthly levels. The dataset records only what a clinic would:
- the prescribed regimen (patients miss 3–15% of doses, sometimes the night before a test);
- noisy levels, with approximate times for 30% of samples;
- creatinine.

Some patients start a thiazide that is never recorded. The patients come from the same deliberately different "truth" model as the simulations. The run (700 forecast levels) shows what the report looks like and how the comparators behave:

| Method | Bias (MDPE) | Imprecision (MDAPE) | Mean abs. error, mmol/L (95% CI) | Within ±20% | 90% PI coverage |
|---|---|---|---|---|---|
| Engine, Bayesian forecast from earlier levels | −1.0% | 15.2% | 0.124 (0.115–0.135) | 63% | 77% |
| Engine's population prior only | +29.1% | 41.4% | 0.310 (0.273–0.349) | 28% | 66% |
| Previous level, unchanged | −14.3% | 28.6% | 0.204 (0.189–0.219) | 37% | – |
| Previous level × change in daily dose | −2.2% | 17.3% | 0.150 (0.136–0.163) | 55% | – |

**What this shows about the harness:**
- **It separates the methods clearly.** The engine vs "previous level" difference in mean absolute error was −0.080 mmol/L (95% CI −0.095 to −0.068).
- **The strata are informative:**
  - when the dose had changed, the engine's error was about half the previous level's (0.114 vs 0.243);
  - when the dose had not changed, it was barely better (0.137 vs 0.153).
- **It exposes miscalibration.** The engine's 90% intervals covered only 77% of measured levels, dropping to 69% after gaps of more than 90 days. The causes are missed doses, unrecorded interactions and drift larger than the engine assumes.
  - Real data will have the same problems. Interval calibration is therefore a secondary endpoint (section 4), and the error model may need a term for adherence.
- **It catches weak alerts.** Only 10% of the levels above 1.0 were anticipated (P > 1.0 ≥ 20%). Most were caused by the unrecorded thiazide or by drift. That is a reminder that forecasts need the events to be recorded.

None of this is evidence about real patients: both the data and the misspecification are invented.

## 3. Where real individual-level data could come from

| Source | What it has | Access | What it would test |
|---|---|---|---|
| **MIMIC-IV** (Beth Israel Deaconess, Boston) | Inpatients; eMAR administration times; lithium and creatinine | PhysioNet credentialing: training plus a data use agreement. The converter here is ready | Forecasting with exact dose times, during acute illness (a hard test). Home doses before admission are unknown; eMAR covers only part of the database period |
| **NIMH Data Archive**: Bipolar CHOICE, LiTMUS (lithium arms) | Trial data with lithium levels; CHOICE recorded sample timing [E22] | Institutional data access request (check the specific studies' availability) | Outpatient titration and the "any-time" level, in adults on protocolised doses |
| **Central Denmark** (Jacobsen, Köhler-Forsberg, Østergaard; the eLi12 group) | 52,837 lithium tests with timing [E22]; eLi12 development data [E47] | Collaboration | Standardising mistimed levels head-to-head with eLi12, in routine care |
| **Marseille TDM service** (Lereclus, Guilhaumou) | The external evaluation of 10 published models [E39] | Collaboration | A like-for-like benchmark against published models |
| **Sweden** (Millischer et al., n = 2,357) | Population PK with renal and co-medication covariates [E40] | Collaboration | Covariate effects (diuretics, RAAS agents, eGFR) |
| **UK**: Norfolk lithium register (Kirkham) [E20, E29]; CPRD | Long-term levels and eGFR; primary-care prescribing (dose times unknown) | Collaboration; CPRD research approval | Long-term drift, kidney-function decline, recall |
| **Australia**: CADE Clinic (Malhi, Royal North Shore) [E70]; pathology providers | Specialist-clinic cohorts; community levels | Collaboration and HREC approval | The target setting for the product |

There are, in practice, no openly downloadable individual-level lithium PK datasets. Individual concentration–time data can sometimes be digitised from figures in older PK papers. That would test the structural model on a handful of subjects, but not the forecasting task.

## 4. Pre-specified plan for the first real dataset

Write this down before looking at the data. That is what makes a result believable.

1. **Freeze the engine.** Record the git commit. No tuning on the validation data. If recalibration is needed, split the data by patient: develop on one part, then report on held-out patients.
2. **Primary endpoint.** Mean absolute error of the engine's forecast vs "previous level, unchanged", for levels with at least one earlier level. The engine passes only if the patient-bootstrap 95% CI of the difference excludes zero in its favour.
3. **Secondary endpoints:**
   - relative bias within ±10%;
   - 90% prediction-interval coverage of 85–95%;
   - share within ±20%;
   - anticipation of levels above 1.0 mmol/L.
4. **Report every stratum,** including those where the engine does worse. The expected hard cases are:
   - samples within 6 h of a dose;
   - older adults;
   - acute kidney injury;
   - long gaps between levels.
5. **Publish the aggregate report** with the commit hash, whatever it shows.

## 5. Data governance

- **Never commit patient-level data.** Credentialed data (MIMIC-IV, NDA) are governed by their data use agreements.
  - The harness writes only aggregate results by default.
  - The optional per-level output (`--rows`) must stay with the data.
- **Do not paste patient-level rows into online AI tools,** including the assistant that built this prototype. PhysioNet's credentialed-data agreement restricts sharing with third parties, and PhysioNet has published guidance on online AI services; read it first.
  - Run the harness locally, or in the approved environment.
  - Share only the aggregate report.
- **Ethics.** A human research ethics committee (Australia) or IRB approval, or a documented exemption, is needed for any non-public dataset. Local institutions may require notification even for de-identified public data.
