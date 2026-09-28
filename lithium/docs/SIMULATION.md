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

(Inserted after the full run.)
