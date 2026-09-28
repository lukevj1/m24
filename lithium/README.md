# Lithos: making lithium easy to prescribe

A design proposal and a working research prototype for next-generation lithium monitoring and dosing: something like DoseMeRx for vancomycin, rebuilt around what makes lithium hard.

> **Research prototype. Not a medical device and not for clinical use.** Model parameters are provisional until externally validated.

## Start here

| Read | What's in it |
|---|---|
| [`docs/DESIGN.md`](docs/DESIGN.md) | The answer to "what would you build?": the problem, why earlier fixes failed, design principles, the patient journey, the engine, the standardised 12-hour level, architecture, regulation, evidence plan, adoption and roadmap |
| [`docs/EVIDENCE.md`](docs/EVIDENCE.md) | The evidence brief behind it: 70 tagged findings, each with a verification label |
| [`docs/SIMULATION.md`](docs/SIMULATION.md) | In-silico trials: methods, results and caveats |

## The prototype (`lithos/`)

A Bayesian lithium engine in plain Python (numpy + scipy):

| Module | What it does |
|---|---|
| `pk.py` | Exact linear 1-/2-compartment PK with IR/SR/liquid absorption and time-varying parameters (checked against an ODE solver) |
| `models.py` | Population priors with provenance: a kidney-function two-compartment model, Methaneethorn & Sringam 2019, a development placeholder and a simulation-only "truth" |
| `bayes.py` | MAP + Laplace individualisation; two-timescale clearance drift (Gaussian process); timing-aware error model; model-averaging ensemble |
| `forecast.py` | Standardised 12-hour level from any sample, half-life, concentration–time bands |
| `dosing.py` | Tablet-aware dose finder by probability of target attainment, "repeat level first" logic, initiation plan, age-based default targets |
| `interactions.py` | Interaction catalogue with sourced effect sizes; patient-specific "what if" forecasts and compensating dose |
| `schedule.py` | NICE CG185 calendar as a floor, plus uncertainty-driven tightening |
| `renal.py` | CKD-EPI 2021, Cockcroft–Gault, robust eGFR trends with NICE NG203/KDIGO flags |
| `pregnancy.py` | Gestation-dependent clearance (Westin/Wesseloo), dose plan, postpartum hazard |
| `safety.py` | Toxicity triage, elimination forecast, EXTRIP criteria |
| `report.py`, `demo.py` | Clinician card and five worked scenarios |

```bash
cd lithium
pip install numpy scipy matplotlib pytest
python -m lithos.demo          # five clinical scenarios, simulated patients with hidden truth
python -m pytest -q            # 28 tests
python -m sim.run_all --n 500  # in-silico trials (~10 min on 4 cores)
```
