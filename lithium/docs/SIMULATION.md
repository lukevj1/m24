# In-silico trials: methods, results and caveats

Reproduce with `cd lithium && pip install numpy scipy matplotlib && python -m sim.run_all --n 500`
(about 10–15 minutes on 4 cores). Raw per-patient output goes to `sim/results.json` and the summary tables to
[`SIMULATION_RESULTS.md`](SIMULATION_RESULTS.md).

## Why simulate first

Before asking patients or clinicians to trust a model, we want to know whether the *mechanisms* behave as intended.
Three claims need testing:

- Any-time levels can be standardised.
- Bayesian titration converges in fewer tests.
- Interaction forecasts prevent high levels.

The test must be adversarial: the engine must not be allowed to "know" the truth.

## Design

**Virtual population.** Adults resembling a bipolar-disorder clinic:
- 80% of patients are aged 42 ± 12 and 20% are older adults, 74 ± 5.
- 55% are women.
- BMI is 29 ± 5.5.
- eGFR declines with age (mean 112 − 0.85 × (age − 20), SD 14, range 35–135).

**Truth ≠ engine (deliberate misspecification).** Patients are generated from a `truth-perturbed` model that differs from the engine's prior in every way that matters:

| | Truth (generates patients) | Engine prior |
|---|---|---|
| Structure | Two-compartment (distribution phase) | One-compartment (development prior) or published model |
| Renal covariate | Clearance ∝ (de-indexed eGFR)^0.75, 10% lower in women | Clearance ∝ Cockcroft–Gault CrCL, no sex effect |
| Between-subject variability | CL 32%, V 30% | CL 25%, V 25% |
| Within-person drift | Slow SD 0.12 (τ 540 d) + fast SD 0.15 (τ 14 d) | Slow SD 0.10 (τ 730 d) + fast SD 0.12 (τ 21 d) |
| IR absorption | ka 1.5 /h | ka 1.2 /h |

**What the engine sees is only what a clinic would see:**
- **Prescribed doses.** The engine assumes all prescribed doses were taken. In truth, 80% of patients miss 3% of doses and 20% miss 15%.
- **Measured levels.** Assay CV is 4%, plus an additive SD of 0.02 mmol/L, rounded to 0.01.
- **Sample times as recorded.** Model-informed arms: recorded accurately to ±15 min, since the design captures dose and sample times. Usual care: the clinician assumes a 12-hour level.

**Realistic sampling times** after a 21:00 dose:
- 60% in the morning (11–13.5 h);
- 30% in the late morning or afternoon (13.5–19 h);
- 10% after work (19–22 h).

## Experiments

1. **Titration to 0.6–0.8 mmol/L** (once-nightly IR lithium carbonate, 250 mg tablets, 8-week horizon). Three arms:
   - *Usual care, stepwise.* Start 500 mg (250 mg if 65 or over). First level at 1 week, read at face value. Adjust by ±250 mg (±500 mg if < 0.4 or > 1.0). Re-check a week after each change. Stop after two consecutive in-range levels.
   - *Usual care, proportional* (a strong comparator). As above, but new dose = old × 0.7 / level, rounded to 250 mg and capped at ±500 mg per step.
   - *Model-informed.* Model-chosen lead-in (about half the maintenance dose for 3 days), then the model-predicted maintenance dose. Levels at convenient times with recorded timing. Bayesian update and recommendation after every level. Stop after two consecutive "no change" decisions.

   Outcomes: blood tests, days to confirmed stability, whether the *true* 12-hour level on the final dose is in range, days on a dose whose true level exceeds 1.0, and number of dose changes.
2. **Any-time sampling.** Stable once-nightly patients (true 12-hour level spread across 0.4–1.1). One level drawn uniformly 2–22 h after the evening dose. Truth is what a correctly timed 12-hour sample would have shown that morning. Compared: the raw level read at face value, the engine with that level only, and the engine with one earlier level on record.
3. **Thiazide started in a stable patient.** The true clearance effect is drawn with median 0.70 (90% range 0.50–0.90); the engine's catalogue assumes 0.72 (0.55–0.90). Compared: no lithium change vs the engine's pre-emptive dose adjustment. Outcome: the true 12-hour level after starting.
4. **Early sampling after a dose change.** One earlier level on the old dose (or none). The dose is increased by about 40%. A level is drawn on day 1, 2, 3, 4, 5 or 7. Outcome: error in the estimated new steady-state 12-hour level.

## Caveats (read before quoting any number)

- These are simulations. They test whether the *mechanisms* work under stated assumptions; they are not evidence of clinical effectiveness.
- The truth model is invented to be *different*, not to be *right*. Real-world misspecification could be larger or smaller.
- Usual care is stylised. Real practice is more variable, both better (experienced lithium clinics) and worse (lost follow-up, no recheck).
- Adherence, sampling behaviour and drift magnitudes are assumptions that must be estimated from real data.
- The engine's development prior is physiology-based until published models are fitted and validated. Results are reported for whichever prior `run_all` used (see the header of `SIMULATION_RESULTS.md`).

## Results

Full tables: [`SIMULATION_RESULTS.md`](SIMULATION_RESULTS.md). Figures: `docs/figures/`. Headline numbers are summarised in
[`DESIGN.md` §11](DESIGN.md#11-what-the-prototype-shows).

### 2. A level drawn at a convenient time

![Any-time sampling](figures/anytime.png)

Across 500 stable patients sampled 2–22 h after a night-time dose:

- **Error.** Median absolute error against a correctly timed 12-hour level fell from **0.098 mmol/L** (reading the level at face value) to **0.037** (engine).
- **Classification.** The share of samples put in the wrong category (below / in / above 0.6–0.8) fell from **32% to 16%**.
- **Calibration.** The engine's 90% intervals contained the truth **95%** of the time.
- **Time-of-day effect.**
  - Late samples (18–22 h): the gain was largest, error 0.135 → 0.039 and wrong category 43% → 14%.
  - Samples at 10–14 h: the two methods were equal, as they should be.
  - Samples within 6 h of a dose: face-value reading was badly wrong (error 0.32; wrong category 59%). The engine reduced this to 0.074, but still misclassified 36%, and its interval coverage fell to 87%. This is why the product asks for samples at least 6 h after a dose and returns "repeat level first" when a sample can't support a decision.

**An honest negative finding.** Adding a second, older level (from 45 days earlier) made the estimate slightly *worse* (0.048 vs 0.037). The virtual patients' clearance drifts more between visits than the engine assumes, so the engine over-trusts old data. The sensitivity analysis below tests this explanation. The practical lesson: drift magnitudes must be estimated from real longitudinal data before history is trusted.

### 3. A thiazide is started in a stable patient

![Thiazide](figures/thiazide.png)

| True 12-hour level afterwards | No lithium change | Engine-adjusted dose |
|---|---|---|
| > 1.0 mmol/L | 53% | 22% |
| > 1.2 mmol/L | 29% | 6% |
| 0.5–0.9 | 33% | 61% |
| < 0.5 | 1% | 5% |

The pre-emptive adjustment halved the share above 1.0 and cut the share above 1.2 almost five-fold. It is not a substitute for the check level, though. The truth's interaction effect was deliberately drawn stronger and more variable than the engine's catalogue (median clearance ×0.70 vs ×0.76), and the engine's 90% forecast interval for the no-change level contained the truth in only **68%** of patients. The design therefore always pairs a forecast with a level 5–7 days later, and interaction effect sizes should be learned from outcomes as data accrue.
