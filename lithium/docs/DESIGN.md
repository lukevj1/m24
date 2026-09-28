# Lithos: designing lithium care that makes lithium easy to prescribe

> **Status:** design proposal plus a working research prototype. Not a medical device and not for clinical use.
> Model parameters are provisional until externally validated. Evidence tags such as [E21] refer to
> [`EVIDENCE.md`](EVIDENCE.md), where every entry carries a verification label. Simulation methods and results are in
> [`SIMULATION.md`](SIMULATION.md).

*Lithos* is Greek for "stone". Lithium was named after it because the element was first found in a mineral. It is a working name.

---

## The short answer

If I were building the lithium equivalent of DoseMeRx, I would **not build a better dosing calculator**. Bayesian lithium forecasting has been published since the late 1980s [E43], and a priori dosing equations and web calculators exist. Yet none of the established precision-dosing products lists lithium [E50], and prescribing has kept falling [E5].

The obstacle is not the maths. It is the *system*:

- Lithium is taken outpatient, for decades.
- Care is split between a psychiatrist, a GP, a pathology lab, a pharmacy and the patient.
- Nobody owns the recall.
- About half of real-world levels are not drawn at the 12 hours the rules require [E22].
- Only about one patient in six gets fully guideline-compliant monitoring [E17].
- The dangerous moments happen far from whoever understands lithium: a new blood-pressure tablet, gastroenteritis, a heatwave, childbirth, slowly failing kidneys.

So I would build **a lithium care system with a Bayesian engine at its core**. It would do five things that, together, no previous tool does:

1. **Every blood test counts, whenever it is taken.** The engine converts a level drawn at any convenient time, including before steady state, into the standardised 12-hour level every guideline is written in. It gives an honest uncertainty interval and says when a sample is too uninformative to act on.

   This builds directly on Amdisen's standardised 12-hour level [E11] and on the recent **eLi12** work [E47], which showed that estimating the 12-hour level from mistimed samples cuts the error from 25% to 10%. What the engine adds: estimates are individualised, uncertainty-quantified and informed by the patient's own history, and they work when the patient isn't yet at steady state.
2. **The right dose in fewer tests.** The first dose is chosen from kidney function, age and size. Each level updates a personal model, so titration converges in fewer steps and doesn't wait for steady state after every change.
3. **Forecasts, not alerts, at the dangerous moments.** Take a GP choosing ramipril for a patient on lithium. The software forecasts *this* patient's new level and proposes the compensating dose. It books the check level inside the GP's own prescribing screen and flags the first month, when ACE inhibitors and loop diuretics carry a 5–8-fold risk of toxicity admission in older adults [E23].

   Generic alerts have been tried. In one US cohort an interaction alert had fired for every interacting drug that preceded a toxicity episode, and starting one still carried about 30-fold odds of toxicity [E21]. The same forecasting covers sick days, heatwaves, surgery and delivery.
4. **Monitoring follows risk, not the calendar.** Guideline intervals are the floor. The model brings tests forward for people near the ceiling, on interacting drugs or with falling eGFR. It tracks the kidney trajectory for life, because higher levels and toxicity episodes predict renal decline [E28, E29].
5. **Nobody has to remember.** Recall, pathology requests, reminders, results, interpretation, sign-off and patient instructions form one closed loop. Patient, GP, psychiatrist, lab and pharmacist share one timeline.

How it gets adopted matters as much as what it does. Vancomycin Bayesian dosing had a four-society guideline naming it the preferred method, and still reached under half of surveyed US hospitals a year later [E51]. Lithium has no such mandate and no commercial sponsor. It is a cheap generic. So:

- **Start where adoption is frictionless:** pathology labs (a standardised 12-hour level on every lithium report) and specialist lithium clinics.
- **Generate the evidence:** external validation, then a pragmatic trial.
- **Make the standardised level a reporting standard** with guideline and laboratory bodies: lithium's "INR moment".
- **Keep the engine an open, public-good core.**

This is also the practical answer to Professor Gin Malhi's long-running argument that lithium deserves wider use [E65–E70]. He has written on interactions, target levels (the "lithiumeter") and the balance of renal risk against benefit. Most of the fear he argues against is about *managing* lithium, and management is what software can take off people's hands.

---

## 1. Why monitoring is what stops people prescribing lithium

Lithium is among the most effective maintenance treatments for bipolar disorder, with the longest time to treatment failure among common alternatives in UK records [E1] and a distinctive anti-suicide signal [E2]. Its use has declined worldwide over two decades [E5]. The reasons offered are overwhelmingly about *managing* lithium rather than whether it works: toxicity, kidneys, the burden of monitoring and confidence [E6, E70]. (No quantified survey of prescriber barriers could be found. That is a gap, and an easy first study.)

| Stage | What makes it hard today | What happens as a result |
|---|---|---|
| Starting | Baseline bloods; a rule-of-thumb dose; a level a week after starting and after every change, weekly until stable [E7] | Weeks of weekly tests; slow titration; side effects before benefit |
| Every level | Supposed to be drawn 12 h after the evening dose [E8, E10], but about half are not: 44.9–49.7% fall outside 10–14 h, and it's worst in general practice [E22] | Mistimed samples are misread, repeated, or push the dose the wrong way |
| Maintenance | 3-monthly levels in year 1; then 6-monthly, or 3-monthly for a long list of higher-risk groups; renal, thyroid, calcium and weight every 6 months, for life [E7] | Fully compliant monitoring ~16% [E17]; annual level standard met by 30% (48% after a national programme) [E14]; 69% of tests in stable patients requested at the "wrong" interval [E15]; 42% of levels below range in one service [E16] |
| Other prescribers | NSAIDs, ACE inhibitors, ARBs and diuretics come from people who don't manage the lithium | ~30-fold odds of toxicity after starting an interacting drug, with alerts firing [E21]; ACE inhibitor RR 7.6 and loop diuretic RR 5.5 for admission in the first month in older adults [E23] |
| Illness and life | Vomiting, diarrhoea, fever, heat, fasting, surgery; pregnancy and delivery | Acute-on-chronic toxicity between appointments; clearance up ~50% in pregnancy then falling abruptly at delivery [E33–E35] |
| Decades | Slow eGFR decline, hypothyroidism, hypercalcaemia | 29.5% of long-term users have an eGFR < 60 at some point (end-stage failure is rare) [E28]; kidney fear drives avoidance and discontinuation [E6] |

## 2. Why previous fixes did not fix it

**Registers, record books and reminders** treat the problem as compliance. England's 2009 safety alert [E13] and a national audit-and-feedback programme moved annual level monitoring from 30% to 48% over years [E14]. The only randomised trial of EHR reminders missed its primary endpoint [E18]. A county register with active recall reported no known toxicity from inadequate monitoring over a decade [E20], which shows what reliable recall can do but also how local such successes remain. These tools address *whether* a test happens, not *what it means* or *what dose follows*.

**Dosing equations and calculators** (Cooper, Pepin, Zetin, Terao, Abou-Auda) are biased and imprecise in head-to-head comparisons [E43]. Newer work is promising but early:
- a one-compartment benchmark predicted the dose within one 200 mg tablet in a small sample [E45];
- machine learning reached RMSE 0.17–0.20 mmol/L on inpatients [E44];
- one app has a trial protocol [E46];
- eLi12 standardises mistimed samples [E47].

None is embedded in routine workflow, and none handles the event-driven risks.

**Published population models predict poorly on new patients.** In an external evaluation, median absolute errors across ten models ranged from 24% to 138%, with most models under-predicting, largely because they lack kidney-function covariates [E39]. Any engine must therefore combine models and learn from each patient, not trust one equation (§5).

**Precision-dosing platforms** succeeded where there is an inpatient drug, a pharmacist who owns dosing, dense sampling and fast measurable harm. None lists lithium [E50].

### Lessons from three analogues

| Analogue | What made it work | Present for lithium? | What Lithos borrows |
|---|---|---|---|
| **Vancomycin AUC dosing** [E51] | Guideline named Bayesian software as preferred; pharmacist owner; EHR-embedded apps; any-time levels (no true trough needed); nephrotoxicity fell (8% → 0–2% in one trial) | No mandate, no single owner, outpatient, fragmented data | Any-time sampling; EHR-embedded forecasts; a deliberate guideline strategy, because even a mandate reached < 50% adoption in a year |
| **Warfarin** [E52] | The INR standardised the measurement; anticoagulation clinics; time-in-therapeutic-range (TTR) as the quality metric; self-testing and self-management for selected, trained patients reduced clots and deaths | The 12-hour rule is a human workaround, not a standard; no TTR equivalent; validated finger-prick devices exist but aren't integrated [E55] | A standardised 12-hour level on every report; **lithium TTR** as the metric; home/POC testing with model interpretation for selected patients |
| **Clozapine registries** [E53] | "No blood, no drug" enforced monitoring | Hard gating burdens access; the US removed its REMS | *Soft* gating: pharmacy-visible monitoring status and refill-time checks, never a block on a mood stabiliser |

**Conclusion:** the gap is a system, not an algorithm. The algorithm still matters, because it is what makes the system trustworthy and light-touch, but a better calculator alone would change little.

---

## 3. Design principles

1. **Every level counts.** No sample is "wrong". A sample with recorded timing is information, weighted by how informative it is.
2. **Know when not to answer.** If the evidence is too uncertain to justify a change, say so and ask for a better-timed sample instead of changing the dose.
3. **Forecast, don't alert.** Replace "interaction: monitor lithium" with "this patient's level is forecast to rise from 0.68 to 0.90 (0.69–1.20); reduce to 500 mg; level on 14 March". One tap to accept.
4. **The guideline is the floor; uncertainty sets the pace.** Never monitor less often than the configured guideline. Monitor *more* when the forecast warrants. Extending intervals for stable patients is a trial question (§9), though the observational case for it is strong [E19].
5. **Nobody has to remember.** Every obligation (test due, result pending, sign-off needed, patient not reached) is a task with an owner and an escalation path.
6. **Put the answer where the decision is made:** in the pathology report, the GP's prescribing screen, the pharmacy dispensing check and the patient's phone.
7. **Kidneys are a first-class outcome:** the eGFR trajectory, lithium-clearance trajectory and exposure history are tracked for life.
8. **Clinician in the loop; patient as partner.** The engine recommends and a clinician signs. The patient records the one thing only they know (when they took the dose), sees their numbers and gets the right advice at the right moment.
9. **Built to be regulated and trusted:** transparent models and model cards, validation per release, audit trails, post-market monitoring.
10. **Public-good economics.** Open engine and models; funded integration and service. The drug is too cheap for a sponsor and too important to leave without one.

---

## 4. What it does, moment by moment

Each subsection says what the system does and what the prototype in this repository demonstrates (`python -m lithos.demo`).

### 4.1 Deciding to start
- **A pre-lithium work-up** is generated from the record and ordered in one tap: weight/BMI, U&E with calcium, eGFR, TSH, FBC, and an ECG where there is cardiovascular risk [E7]. For anyone who could become pregnant, it adds a pregnancy test and a contraception and pregnancy-planning discussion.
- **Counselling and consent content:**
  - sick-day rules;
  - over-the-counter NSAIDs;
  - fluids and salt;
  - what monitoring involves;
  - why levels matter: relapse is 2.6 times more likely at a median 0.54 than at 0.83 [E3].
- **A baseline kidney-risk profile** (eGFR, albuminuria, diabetes, hypertension, age, nephrotoxic co-medication) sets the target and schedule. The default is once-daily night-time dosing [E32]. The profile also frames an honest kidney conversation from day one.

### 4.2 Starting: the first dose comes from the model
Instead of "start 250–500 mg and creep up weekly", the model ensemble predicts the maintenance dose for the target from eGFR, age, sex and body composition. The system proposes:
- a short **lead-in** at about half that dose, for tolerability;
- then the **predicted maintenance dose**, capped so the forecast probability of exceeding 1.0 mmol/L stays under 5%;
- then the **first level about a week in, at any convenient time** at least 6 h after a dose, with the dose time recorded.

*(Prototype: `dosing.initiation_plan`; demo scenario 2.)*

### 4.3 Getting to target: Bayesian titration
- **Each level updates the individual model.** The recommendation is the real-tablet regimen with the highest posterior probability of a 12-hour level in the target window.
- **The next level is scheduled when it will be informative**, not necessarily at steady state.
- **Uncertainty alone never drives a dose change.** If the only information is, say, a sample drawn 2 h after a dose, the recommendation is "repeat a well-timed level", not "cut the dose". *(Prototype: `dosing.recommend`; demo scenario 1; simulation experiment 1.)*

### 4.4 Staying well: the maintenance loop
- **Every lithium result is reported with its standardised 12-hour equivalent** and an interval, next to the raw value. After a night-time dose, blood can be taken any time the next day.
- **Recall is automatic:**
  - the pathology request is generated;
  - the patient is reminded by app or SMS (carer option);
  - non-attendance escalates;
  - results route to the responsible clinician with a pre-drafted decision.
- **The interval is individual.** NICE's higher-risk list [E7] is derived automatically from the record: age 65+, interacting drugs, renal or thyroid risk, raised calcium, poor control, poor adherence, and a last level ≥ 0.8. The uncertainty forecast can bring the next test forward. *(Prototype: `schedule.plan`.)*
- **Clearance-drift detection.** Each person's lithium clearance is tracked over time, and a sudden change prompts a specific question:
  - clearance down about 30% since the last check: "dehydration? new NSAID or ACE inhibitor? renal decline?";
  - clearance up: "missed doses? a formulation change? pregnancy? an SGLT2 inhibitor [E26]?"

  *(Prototype: two-timescale drift model and diagnostics in `bayes.py`.)*

### 4.5 When something changes
- **Interacting drugs at the point of prescribing and dispensing.** CDS Hooks services on `order-select`/`order-sign` and on `medication-refill` [E63] forecast this patient's new level and propose a dose and a level date. The two kinds of evidence are treated differently:
  - **thiazides** have the largest average pharmacokinetic effect [E24], so they get *dose arithmetic*;
  - **ACE inhibitors and loop diuretics** carry the largest first-month admission risk in older adults [E23], so they get a *time-anchored first-month monitoring plan*, whatever the forecast says.

  *(Prototype: `interactions.what_if`; demo scenario 3; simulation experiment 3.)*
- **Sick days.** A one-tap "I'm unwell" in the app runs a short symptom check (vomiting, diarrhoea, fever, reduced intake, worse tremor, unsteadiness, confusion). It then:
  - returns the local sick-day rule [E10] and notifies the team;
  - offers a level and U&E;
  - resumes the plan when the person is well.
- **Heat.** Local heatwave warnings (in Australia, from the Bureau of Meteorology) trigger fluid and salt messages to patients in the affected region, and flag higher-risk patients to clinicians.
- **Surgery and fasting.** A peri-operative plan (withhold and restart criteria) goes to the surgical team.

### 4.6 The long game: kidneys, thyroid, calcium
- **Renal trend analysis follows NICE NG203 and KDIGO** [E31]:
  - a robust eGFR slope needs ≥ 3 values over ≥ 90 days;
  - a new low is repeated within 2 weeks to exclude acute injury;
  - accelerated-progression and rapid-decline definitions are applied;
  - referral criteria follow NG203.
- **Lithium-specific actions** [E28, E29, E32]: lowest effective level, once-daily dosing, avoiding toxicity episodes, and more frequent monitoring. A single level > 1.0 was followed by a measurable eGFR fall [E29].
- **The continue-or-stop decision is shared, not reflexive.** It is framed as NICE recommends: renal and bipolar-specialist advice, weighed with the patient [E31]. Stopping carries relapse and suicide risk [E4].
- **TSH and calcium trends**, with prompts for PTH; weight; screening questions for polyuria and polydipsia.
- **Exposure metrics are computed for every patient:** cumulative dose, time above 1.0 and AUC. Which exposure metric best predicts renal harm is an open research question (§14). *(Prototype: `renal.egfr_trend`; demo scenario 4.)*

### 4.7 Pregnancy and after delivery
Clearance rises through pregnancy and falls within days of delivery, when relapse risk also peaks. The engine carries both published descriptions of the curve, with between-woman uncertainty:
- Westin 2017 (full text): about 1.2% per week, −34% at week 34 [E34];
- Wesseloo 2017: trimester nadir −36% [E33].

It turns them into a week-by-week forecast dose plan for the perinatal team and schedules:
- levels every 4 weeks, weekly from 36 weeks, and within 24 h of birth;
- twice weekly for 2 weeks postpartum;
- a first-month postpartum target of 0.8–1.0 [E35].

It surfaces rather than hides the unresolved question of what to do before delivery. US labelling advises reducing or stopping 2–3 days beforehand; a 233-level study found no intrapartum rise and advises against [E35].

It also quantifies the classic hazard: in the prototype, simply continuing a 1,500 mg late-pregnancy dose after delivery forecasts a 12-hour level of about 1.05 (90% CrI 0.77–1.40). *(Prototype: `pregnancy.plan`, `pregnancy.postpartum_hazard`.)*

### 4.8 When things go wrong: toxicity
Triage combines the level, symptoms, acute vs chronic context and kidney function, and follows standard practice:
- symptoms override the number;
- don't wait for a 12-hour sample;
- activated charcoal doesn't bind lithium [E10, E36].

The distinctive addition is an **individual elimination forecast**: hours until the level falls below 1.0 mmol/L with lithium withheld. That is literally one of the EXTRIP criteria for extracorporeal treatment [E27]. It goes to the treating team with its assumptions stated. Thresholds are aligned with the local toxicology service (in Australia, the Poisons Information Centre, 13 11 26). *(Prototype: `safety.triage`; demo scenario 5.)*

### 4.9 Stopping, switching, or continuing with CKD
A structured shared-decision view shows the kidney trajectory, psychiatric history and relapse risk, and the options and consequences, including slow tapering over 4 weeks to 3 months [E4].

---

## 5. The engine

The prototype is in `lithos/`: about 2,000 lines of Python and 28 tests.

| Component | Choice | Why |
|---|---|---|
| Structure | Linear 1- or 2-compartment model, first-order absorption, one depot per formulation (IR/SR/liquid). **Piecewise-constant parameters**; exact superposition via eigen-decomposition | Lithium is unmetabolised, renally cleared and linear [E41]. Time-varying parameters represent interactions, pregnancy, drift and renal change exactly. Verified against an ODE solver to 1e-11. |
| Priors | **An ensemble, weighted per patient by Bayesian evidence:** (A) a kidney-function-based two-compartment model assembled from verified physiology [E41]; (B) Methaneethorn & Sringam 2019 [E38]. More published models can be added. | No single published model predicts well externally [E39]; kidney function and body composition matter most [E37, E40]. Model averaging is more robust than picking one [E54]. |
| Individualisation | MAP estimate + Laplace posterior (optional importance-sampling correction) | Fast (tens of ms per model), well understood, auditable |
| **Clearance drift** | Each blood-test occasion has its own clearance deviation. Deviations are the sum of a slow (renal ageing) and a fast (hydration, illness) Ornstein–Uhlenbeck process, i.e. a Gaussian process in time | Decades of therapy: old levels inform the present less than recent ones; forecast uncertainty grows with time since the last level; sudden departures are detectable |
| **Timing-aware error** | Level variance includes (dC/dt × timing SD)² plus a distribution-phase term that decays with time since the dose | Samples on steep parts of the curve are automatically trusted less. This is what makes any-time sampling safe [E41: distribution can take ~10 h in older adults]. |
| Decisions | Probability of target attainment over real tablet combinations; a cap on P(level > 1.0); step limits; hysteresis to avoid churn; **a "repeat level first" state** | Recommendations are probabilistic, practical and honest about uncertainty |
| Recheck timing | From the individual half-life posterior | Older adults and people with CKD need longer; young adults less |
| Diagnostics | Samples in the absorption or distribution phase; implausibly long since the last dose; large standardised residuals; clearance jumps | Data-quality problems become specific questions a clinician can act on |

**Honest status of the priors:**
- Model A's covariate structure comes from published physiology, not from a fitted dataset.
- Model B's between-subject variability and residual error could not be retrieved and are assumed.
- Both must be refitted and externally validated on real outpatient data with recorded dose times before any clinical use (§9).

**What it deliberately does not do (yet):**
- patient-facing autonomous dose changes;
- act without clinician sign-off;
- use sensor inputs;
- answer when its own uncertainty is too wide.

**Learning system.** The system should learn as it is used:
- priors improve with use: refit population parameters, drift magnitudes and interaction effect sizes from pooled, consented, de-identified data (federated where needed);
- prediction error is monitored per release and per subgroup (age, sex, kidney function, ethnicity, formulation);
- model cards are published.

---

## 6. A standardised 12-hour level on every report: lithium's "INR moment"

The 12-hour rule is a human workaround. Without a model, the only way to compare levels is to sample everyone at the same point on the same curve. Amdisen defined that standard in 1977 [E11]. With a model, the *report* can be standardised instead of the *sample*.

Others have started down this road:
- eLi12 estimates the 12-hour level from mistimed samples (error 10% vs 25%) and has a feasibility trial registered [E47];
- SimpLi shows that purely empirical time corrections explain little [E48];
- a 2025 proposal calls for re-introducing Amdisen-standardised levels [E11].

The natural move is to collaborate with these groups, not to compete.

**Proposal.** Pathology reports carry, next to the measured value:
- the **standardised 12-hour level** on the current regimen, with a 90% interval and a confidence grade;
- the inputs used: dose, regimen, last-dose time and sample time;
- quality flags such as "absorption/distribution phase" or "not at steady state; projection".

Only two extra facts are needed at collection: the regimen and when the last dose was taken. The collection-centre form, the patient's app or dispensing data can supply them.

Why it matters:
- It ends "was that a 12-hour level?".
- It makes levels comparable across time, labs and settings, and supports the age-specific reporting ranges experts have asked for [E9].
- It makes a **lithium time-in-therapeutic-range** computable for benchmarking and trials, as TTR is for warfarin [E52].
- It needs a standards fix. LOINC has lithium codes (14334-7) and a trough variant, but no 12-hour post-dose term, and real datasets map codes and units inconsistently [E62].

Two practical traps the design handles:
1. **Lithium-heparin tubes** falsely raise lithium results, and some home-collection kits ship with them, so kits must be lithium-free [E60].
2. **Once- vs twice-daily dosing** give different 12-hour levels for the same total dose [E49]. The standardised level is always defined *for the regimen actually taken*, and the engine translates between regimens from the actual dose times rather than a fixed conversion factor.

---

## 7. Where it lives: architecture and integration

```
                 ┌───────────────────────── Patient ─────────────────────────┐
                 │  App / SMS / IVR / carer: dose-time taps, sick-day button, │
                 │  reminders, results in plain language, mood & side effects │
                 └───────────────▲───────────────────────────────▲───────────┘
                                 │                               │
 Pathology labs ── HL7 v2 ORU ──►│        LITHOS PLATFORM        │◄── eRx / dispensing (refills, new drugs)
 (Li, U&E, TSH, Ca;              │  Engine (PK/Bayes, models)     │
  12-h equivalent back on report)│  Rules (guidelines as data)    │◄── Hospital EMR (SMART on FHIR app)
                                 │  Orchestrator (tasks, recall)  │
 GP systems ◄── CDS Hooks cards ─│  Register & quality dashboard  │──► My Health Record (care plan, results)
 (order-select / order-sign /    │  Audit, model registry,        │
  medication-refill; recall)     │  post-market monitoring        │──► Psychiatrist / lithium clinic workspace
                                 └────────────────────────────────┘
```

- **Data model:**
  - FHIR `MedicationRequest`/`MedicationStatement` (regimen and patient-reported dose times);
  - `Observation` (LOINC 14334-7 plus timing extensions; eGFR; TSH; calcium);
  - `CarePlan` (target, regimen, schedule, sick-day rules);
  - `Task` (every monitoring obligation, with owner and due date);
  - `Communication`.
- **Integration points:**
  - Pathology HL7 v2 feeds.
  - GP software. In Australia that means Best Practice and MedicalDirector, via partner APIs and HL7 import, and CDS Hooks where supported.
  - Hospital EMRs via SMART on FHIR. DoseMeRx's Oracle Health integration at an Australian hospital shows the path [E50].
  - Electronic prescribing and dispensing data (Active Script List) for refills and new interacting drugs, where access is permitted.
  - My Health Record.
- **Deployment:** a stateless engine service (the prototype is essentially this), an orchestrator, clinician and patient front-ends, and a register. The platform is the system of record for *lithium-specific* state only.

---

## 8. Safety, regulation and governance

| Configuration | Likely position (to be confirmed with regulatory counsel) [E64] |
|---|---|
| **A. Clinician-facing** recommendations from lab results and EHR data; transparent basis; clinician signs; not time-critical | Australia: may meet the TGA clinical-decision-support exemption. US: plausibly non-device CDS; the January 2026 guidance extends enforcement discretion to a single clinically appropriate recommendation. EU: MDR Rule 11 means IIa at least, realistically IIb. |
| **B. Patient-facing** dose instructions or self-titration | A device everywhere. US: most plausibly Class II "drug dose calculator", with insulin-dosing software as the predicate. |
| **C. Ingesting POC, wearable or saliva-sensor signals** | Fails CDS exemptions; the sensor is its own IVD submission |

Strategy:
- Launch as configuration A.
- Build to IEC 62304, ISO 14971 and ISO 13485 from day one regardless.
- Add patient-facing dosing and sensor integration later, as separately cleared modules.

Guardrails:
- **Hard stops:** triage pathways above set levels or with neurological symptoms.
- **Caps:** P(level > 1.0) ≤ 5% for a recommended regimen; dose-step limits.
- **Explicit "repeat level first" and "cannot recommend" states.**
- **Mandatory clinician sign-off** and a full audit trail (inputs, model versions, posterior weights, recommendation).
- **Per-release validation, with field monitoring** of prediction error, drift and subgroup performance.
- **Privacy:**
  - Australian Privacy Principles and the My Health Records Act;
  - explicit consent for the app and for research reuse;
  - data minimisation;
  - Australian hosting.

---

## 9. Evidence: how to prove it, in order

| Stage | Question | Design | Primary measures |
|---|---|---|---|
| 0. In silico (this repo) | Do the mechanisms work under stated assumptions, including a misspecified model? | Virtual trials (§11) | Tests to target; standardisation accuracy; interaction outcomes |
| 1. Retrospective validation | How accurate is the engine on real outpatient data? | External validation on lithium-clinic and TDM-service data with *recorded* dose and sample times; several populations; comparison with eLi12 and published models | Bias, precision, % within ±0.1 mmol/L, interval calibration, category agreement |
| 2. Prospective shadow mode | Would its recommendations have been right? | Runs silently alongside care in 2–3 services | Agreement with later levels; clinician-rated appropriateness; alert burden |
| 3. Pragmatic trial | Does it improve care? | Stepped-wedge cluster RCT across services and GP practices | **Primary:** lithium TTR over 12 months. **Secondary:** tests per patient-year; time to target at initiation; levels > 1.2; toxicity admissions; guideline-concordant monitoring; patient burden; clinician confidence; eGFR slope (long-term); **new lithium starts among eligible patients** |
| 4. Adaptive-monitoring trial | Can stable patients safely test less often? | Randomised non-inferiority, risk-adapted vs fixed intervals (building on [E19]) | Levels > 1.0; toxicity; relapse; tests |
| 5. Standards | Make it normal | Work with RANZCP, ISBD/IGSLi, NICE, RCPA and LOINC | Standardised 12-hour reporting; lithium TTR as a quality indicator |

The trial outcome that matters most is the last secondary one: **does making lithium easy make clinicians prescribe it when it's indicated?**

---

## 10. Getting it used: adoption and money

**Why nobody has done this:**
- the drug is generic and cheap, so there is no sponsor;
- care is fragmented, so there is no single buyer;
- regulatory cost is real;
- vancomycin shows that even a mandate plus mature software yields slow uptake, with licence cost the top barrier [E51].

**Wedges, in order:**
1. **Pathology interpretive reporting.** Labs add the standardised 12-hour level and a plain-language interpretation to every lithium result. Every requesting GP benefits without installing anything. Australian pathology is concentrated in a few large providers, so a handful of partnerships reach most practices.
2. **Specialist services and lithium clinics:** the full workspace, register and orchestration. This is where evidence is generated and where the champions are, for example academic mood-disorder clinics such as Professor Malhi's CADE Clinic at Royal North Shore Hospital [E70].
3. **GP software:** interaction forecasts at prescribing and refill, and the recall register. It spreads through practices already receiving the reports.
4. **Patients:**
   - the app: dose-time capture, sick-day button, results, mood and side-effect tracking;
   - later, home and point-of-care testing for selected, trained patients, as with warfarin self-management [E52]. Validated finger-prick devices exist [E55], and postal dried blood spots are an option [E59], which matters most for rural and remote Australia.

**Money:**
- **An open, public-good core:** engine, models, validation code and model cards.
- **A funded service layer:** integration, hosting, support and quality reporting, paid by health services, primary health networks and pathology providers.
- **Research funding** for the trials.
- **Cost offsets:**
  - fewer wasted and repeated tests: a stability-based interval scheme cut tests by about 15% [E19];
  - fewer toxicity admissions;
  - and, if uptake rises, fewer relapses. That is the big one.

---

## 11. What the prototype shows

<!-- SIMULATION-SUMMARY -->

---

## 12. Risks, failure modes and mitigations

| Risk | How it would hurt | Mitigation in the design |
|---|---|---|
| Wrong or missing dose times | Mis-standardised levels | Patient taps and pharmacy data; timing uncertainty in the error model; distribution-phase flags; "repeat level first" state |
| Model misspecification (older adults, CKD, pregnancy, other ancestries, SR products) | Biased forecasts | Model ensemble weighted by evidence; subgroup validation before release; drift terms absorb individual deviation; field monitoring; conservative caps |
| Over-trust in history when clearance changes | Stale forecasts | Two-timescale drift model; clearance-jump diagnostics. In simulation, adding an older level sometimes *hurt* when the engine underestimated drift, so drift must be estimated from real longitudinal data (§11, §14) |
| Automation bias | Clinicians rubber-stamp | Show inputs, curve and uncertainty; require sign-off; human-factors testing; audit |
| Alert fatigue | Ignored warnings | Forecast only when a decision changes; cards carry a pre-filled action |
| Monitoring less often by accident | Missed toxicity | Guideline floor enforced; interval *extension* only inside a trial |
| Digital exclusion | Inequity | SMS/IVR/carer channels; collection-centre capture of dose time; no smartphone required |
| Liability and regulation | Slow launch | Configuration A first; quality system from day one; clear intended-use statement |
| Business model | Nobody pays | Lab and service wedges; public-good core; trial evidence of cost offsets |
| Data breach | Harm and loss of trust | Minimisation, encryption, access audit, Australian hosting |

---

## 13. Roadmap

| Horizon | Deliverables |
|---|---|
| 0–6 months | Confirm guideline rules against RANZCP 2020 and Australian sources; retrieve and implement further published priors; retrospective validation on 1–2 datasets with recorded times (seek collaboration with the eLi12 group); clinician workspace for one specialist service (register, recall, standardised level, recommendations); dose-time capture by SMS/app |
| 6–18 months | Shadow-mode study; pathology partner pilot of standardised reporting; CDS Hooks interaction forecasts in one GP software environment; pregnancy module with a perinatal service; human-factors work and safety case; quality-management certification |
| 18–36 months | Stepped-wedge trial; patient app general release; POC/DBS pilot for rural and remote sites; standards work (reporting, LOINC, lithium TTR); guideline engagement; start the adaptive-monitoring trial |

---

## 14. Open research questions this platform would answer

1. How large and how fast is within-person drift in lithium clearance? It sets how much history to trust.
2. How accurate is a Bayesian standardised 12-hour level across formulations, twice-daily regimens and older adults, compared with eLi12?
3. Which exposure metric (12-hour level, AUC, peak, time above 1.0, cumulative dose) best predicts renal decline?
4. Are there *personal* therapeutic windows (n-of-1 exposure–response using standardised levels and mood data)?
5. Can person-specific calibration make saliva lithium usable [E57]?
6. Is risk-adaptive monitoring non-inferior to fixed intervals?
7. What do prescribers, GPs and patients actually say stops them? No quantified survey exists [E6].
8. Does frictionless monitoring increase appropriate lithium initiation?
