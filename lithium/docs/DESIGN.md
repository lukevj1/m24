# Lithos: designing lithium care that makes lithium easy to prescribe

> **Status:** design proposal plus a working research prototype. Not a medical device and not for clinical use.
> The engine's parameters are provisional until they are externally validated.
> Evidence tags such as [E12] refer to [`EVIDENCE.md`](EVIDENCE.md). Simulation results are in [`SIMULATION.md`](SIMULATION.md).

*Lithos* is the Greek for "stone". Berzelius named lithium after it because the element was first found in a mineral. It's a working name.

---

## The short answer

If I were building the lithium equivalent of DoseMeRx, I would **not build a dosing calculator**.

Bayesian lithium dosing has been published since the 1980s [E40]. Web calculators exist. Yet none of the roughly ten established precision-dosing products offers lithium [E41], and prescribing has not moved. The obstacle is not the maths. It is the *system*: lithium is taken outpatient, for decades, with care split between a psychiatrist, a GP, a pathology lab, a pharmacy and the patient. No one owns the recall. Every level must be drawn at an awkward time. The dangerous moments (a new blood-pressure tablet, gastroenteritis, a heatwave, childbirth, slowly failing kidneys) happen far from the person who understands lithium.

So I would build **a lithium care system with a Bayesian engine at its core**. It would do five things that no previous tool combines:

1. **Every blood test counts, whenever it is taken.** The engine converts a level drawn at any convenient time, even before steady state, into the standardised 12-hour level that every guideline is written in, with an honest uncertainty interval. The 12-hour rule stops being the patient's problem.
2. **The right dose in fewer tests.** The starting dose is chosen from the patient's kidney function, age and size. Each level then updates a personal model, so titration converges in fewer steps and does not wait for steady state.
3. **Forecasts, not alerts, at the dangerous moments.** When a GP selects ramipril, the software forecasts *this* patient's new level, proposes the compensating dose and books the check level, inside the GP's own prescribing screen. The same applies to sick days, heatwaves, surgery and delivery. (Generic alerts have been tried: in one US health system an interaction alert fired for every interacting drug that preceded a toxicity admission [E22].)
4. **Monitoring follows risk, not the calendar.** Guideline intervals are the floor. The model brings tests forward for patients near the ceiling, on interacting drugs or with a falling eGFR. The kidney trajectory is tracked as a first-class outcome for life.
5. **Nobody has to remember.** Recall, pathology requests, reminders, results, interpretation, sign-off and patient instructions form one closed loop. The patient, GP, psychiatrist, lab and pharmacist share one timeline.

How to get it adopted matters as much as what it does. Vancomycin Bayesian dosing had a four-society guideline *mandating* it and still reached under half of US hospitals a year later [E43]. Lithium has no mandate and no commercial sponsor: it is a cheap generic. So the plan runs in this order:

- start where adoption is frictionless: pathology labs and specialist lithium clinics;
- generate the evidence (validation, then a pragmatic trial);
- work with guideline bodies to make the *model-standardised 12-hour-equivalent level* a reportable standard, lithium's "INR moment";
- keep the engine an open, public-good core.

---

## 1. Why monitoring is what stops people prescribing lithium

Lithium remains the most effective long-term treatment for bipolar disorder and has a unique anti-suicide signal [E1–E4], yet its use has fallen in many countries [E5–E7]. When clinicians are asked why, the reasons are overwhelmingly about *managing* lithium rather than whether it works: toxicity, kidney damage, the burden and complexity of monitoring, and lack of confidence [E8–E10]. Every one of those is a design problem.

| Stage | What makes it hard today | What happens as a result |
|---|---|---|
| Starting | Baseline bloods; choosing a dose by rule of thumb; weekly correctly timed levels until stable | 4–8 weeks of weekly tests; slow titration; early side effects before any benefit |
| Every level | Must be drawn ~12 h after the evening dose, at steady state, before any morning dose | Mistimed samples are common; they are misread or repeated, and misreading can push the dose the wrong way |
| Maintenance | 3-monthly levels in year 1, then 3- to 6-monthly for life; renal, thyroid and calcium tests 6-monthly | Recall fails. UK audits met the annual level standard for only 30% of patients, and 48% after a national quality programme [E18]; 69% of tests in stable patients were requested at the wrong interval [E19] |
| Interpretation | Is 0.55 low enough to change? Was that a 12-hour level? Is this eGFR fall lithium? | Low confidence, especially in primary care; defensive practice |
| Other prescribers | NSAIDs, ACE inhibitors, ARBs and thiazides are prescribed by people who don't manage the lithium | Starting an interacting drug was associated with ~30-fold odds of toxicity needing acute care [E22] |
| Illness and life | Vomiting, diarrhoea, heat, fasting, surgery, pregnancy and delivery | Acute-on-chronic toxicity, often in the gaps between appointments |
| Decades | Slow eGFR decline, hypothyroidism, hypercalcaemia | Fear of kidney damage dominates decisions to avoid or stop lithium [E8, E30] |

## 2. Why previous fixes did not fix it

**Registers, record books and reminders** treat the problem as compliance. England's 2009 patient-safety alert [E17] and a national quality programme lifted annual level monitoring from 30% to 48% over years [E18]. The one randomised trial of EHR reminders missed its primary endpoint [E20]. These tools address *whether* a test happens, not *what the result means* or *what dose should follow*.

**Dosing calculators** mostly implement a priori equations from the 1970s–90s (Cooper, Zetin and others) [E42]. Some newer academic models are promising: a one-compartment benchmark predicted the required dose within one 200 mg tablet in small samples [E44], and one app has a randomised-trial protocol but no results yet [E45]. None is validated in routine care, none is embedded in the workflow, and none handles the event-driven risks.

**Precision-dosing platforms** (DoseMeRx, InsightRx, PrecisePK, MwPharm++, BestDose and others) succeeded where there is an inpatient drug, a pharmacist who owns dosing, dense sampling and fast measurable harm. None of these vendors lists lithium [E41].

**Interaction alerts** fire generically and are overridden. In the Kaiser Permanente Colorado cohort, an alert had fired for every interacting drug that preceded a toxicity episode [E22].

### Lessons from three analogues

| Analogue | What made it work | Present for lithium? | What Lithos borrows |
|---|---|---|---|
| **Vancomycin AUC dosing** [E43] | Guideline named Bayesian software as preferred; pharmacist owner; EHR integration; any-time levels (no need for a true trough) | No mandate; no single owner; outpatient; fragmented data | Any-time sampling; EHR-embedded apps; a guideline strategy (because even a mandate got <50% uptake in a year) |
| **Warfarin** [E46, E47] | INR standardised the measurement; anticoagulation clinics; time-in-therapeutic-range (TTR) as a quality metric; self-testing and self-management for selected patients reduced clots and deaths | No standardised "interpretable level" (the 12-hour rule is a human workaround); no TTR equivalent; POC devices exist but are not integrated [E48] | **Li12e** as an INR-like standard; **lithium TTR** as the metric; home/POC testing with model interpretation for selected patients |
| **Clozapine registries** [E49] | "No blood, no drug" enforced monitoring | Hard gating adds access barriers; the US removed its REMS | *Soft* gating: pharmacy-visible status and refill-time checks, never a block on a mood stabiliser |

**Conclusion:** the gap is a system, not an algorithm. The algorithm matters, because it is what makes the system trustworthy and light-touch, but a better calculator alone would change little.

---

## 3. Design principles

1. **Every level counts.** No sample is "wrong". A sample with a recorded time is information; the engine weighs it by how informative it is.
2. **Forecast, don't alert.** Replace "interaction: monitor lithium" with "this patient's level is forecast to rise from 0.72 to 0.98 (0.80–1.21); reduce to 750 mg; check a level on 14 March". One tap to accept.
3. **The guideline is the floor; uncertainty sets the pace.** The software never monitors less often than the configured guideline. It monitors *more* when the forecast warrants. Extending intervals for stable patients is a research question, not a default (§9).
4. **Nobody has to remember.** Every obligation (test due, result pending, sign-off needed, patient not reached) is a task with an owner and an escalation path.
5. **Put the answer where the decision is made:** in the GP's prescribing screen, the pathology report, the pharmacy dispensing check and the patient's phone. Not in a separate psychiatric app nobody opens.
6. **Kidneys are a first-class outcome.** Track the eGFR trajectory and lithium clearance for life. Make the renal conversation explicit, early and shared.
7. **Clinician in the loop; patient as partner.** The engine recommends and a clinician signs. The patient records the one thing only they know (when they took the dose), sees their numbers and gets the right advice at the right moment.
8. **Show the uncertainty; know when not to answer.** Every number carries an interval. The engine says "I can't interpret this sample reliably" rather than guess.
9. **Built to be regulated and trusted.** Transparent models, validation per release, audit trails, post-market monitoring.
10. **Public-good economics.** An open engine and open models, with a funded service layer. The drug is too cheap for a pharma sponsor and too important to leave without one.

---

## 4. What it does, moment by moment

Each subsection says what the system does and what the prototype in this repository already demonstrates (`python -m lithos.demo`).

### 4.1 Deciding to start
- A pre-lithium checklist is generated from the record and ordered in one tap: U&E and eGFR, TSH, calcium, weight/BMI, ECG where indicated, and a pregnancy test and contraception discussion for anyone who could become pregnant.
- Structured counselling and consent: sick-day rules, interacting over-the-counter drugs, fluid intake and what monitoring will involve. Patient materials are available in multiple languages.
- **Baseline renal risk profile:** eGFR, albuminuria, diabetes, hypertension, age and nephrotoxic co-medications. It is used to set the target range and dosing schedule (once-daily at night by default) and to frame the kidney conversation honestly from day one.

### 4.2 Starting: the first dose comes from the model
The population model predicts the maintenance dose for the target from eGFR, age, sex and size. The system proposes:
- a short **lead-in** at about half that dose, for tolerability;
- then the **predicted maintenance dose**, capped so the forecast probability of exceeding 1.0 mmol/L stays under 5%;
- then the **first level about a week in, at any convenient time** of day with the last dose time recorded.

This replaces "start 250–500 mg and creep up weekly". *(Prototype: `dosing.initiation_plan`; demo scenario 2.)*

### 4.3 Getting to target: Bayesian titration
Each level updates the individual model. The recommendation is the regimen, built from real tablet strengths, with the highest posterior probability of a 12-hour level in the target window. The next level is scheduled when it will be informative, which is not necessarily at steady state. *(Prototype: `dosing.recommend`; simulation experiment 1.)*

### 4.4 Staying well: the maintenance loop
- Every result is reported as **Li12e**, the model-standardised 12-hour-equivalent level (§6), alongside the raw value. The patient can have blood taken any time the next day after a night-time dose.
- **Recall is automatic.** The pathology request is generated, the patient is reminded by app or SMS (with a carer option), non-attendance escalates, and results route to the responsible clinician with a pre-drafted decision.
- **Clearance drift detection.** The engine tracks each person's lithium clearance over time. A sudden change prompts a specific question: "Clearance about 30% lower than at the last check: dehydration? new NSAID/ACE inhibitor? renal decline?" Or in the other direction: "missed doses? formulation change? pregnancy?" *(Prototype: two-timescale drift model and diagnostics in `bayes.py`; demo scenario 5.)*

### 4.5 When something changes
- **Interacting drugs at the point of prescribing.** A CDS Hooks service on `order-select` and `order-sign` (and `medication-refill` at dispensing) forecasts this patient's new level and proposes a dose and a level date [E50]. *(Prototype: `interactions.what_if`; demo scenario 3; simulation experiment 3.)*
- **Sick days.** A one-tap "I'm unwell" in the app runs a short symptom check (vomiting, diarrhoea, fever, reduced intake, tremor, unsteadiness, confusion). It returns the local sick-day rule (typically: withhold lithium and seek a level and U&E), notifies the team, and resumes the plan when well.
- **Heatwaves.** Local weather alerts (in Australia, Bureau of Meteorology heatwave warnings) trigger hydration messages to patients in the affected region, with higher-risk patients flagged to their clinicians.
- **Surgery and fasting:** a peri-operative plan (withholding and restart criteria) is pushed to the surgical team.

### 4.6 The long game: kidneys, thyroid, calcium
- A **robust eGFR slope** with a confidence range and lithium-specific decision points [E30–E32]: lower target, once-daily dosing, more frequent monitoring, nephrology input, and when to discuss alternatives.
- The system shows what *not* acting costs as well as what acting costs. Stopping lithium carries relapse and suicide risk, and the decision is shared, not reflexive.
- TSH and calcium trends, with prompts for PTH. Weight trend. Screening questions for polyuria and polydipsia (nephrogenic diabetes insipidus).
- Exposure metrics (cumulative lithium, time above 1.0, AUC) are computed for every patient. This is a research asset, because which exposure metric best predicts renal harm is an open question. *(Prototype: `renal.egfr_trend`; demo scenario 4.)*

### 4.7 Pregnancy and after delivery
Clearance rises through pregnancy and falls abruptly at delivery, when relapse risk is also at its peak [E33–E35]. The engine carries a gestation-dependent clearance curve with between-woman uncertainty and produces a week-by-week forecast dose plan for the perinatal team. It schedules levels every 4 weeks and then weekly near term. It generates an explicit **delivery plan**:
- what to do at the onset of labour;
- the postpartum dose;
- a level within 24 hours of delivery.

It also quantifies the classic hazard: continuing the late-pregnancy dose after delivery. *(Prototype: `pregnancy.plan`, `pregnancy.postpartum_hazard`.)*

### 4.8 When things go wrong: toxicity
Triage combines the level, symptoms, acute vs chronic context and kidney function. The distinctive addition is an **individual elimination forecast**: hours until the level falls below 1.0 mmol/L with lithium withheld. That is literally one of the EXTRIP criteria for extracorporeal treatment [E36]. The forecast goes to the treating team with its assumptions stated. Thresholds are aligned with the local toxicology service (in Australia, the Poisons Information Centre, 13 11 26). *(Prototype: `safety.triage`; demo scenario 5.)*

### 4.9 Stopping, switching, or continuing with CKD
A structured shared-decision view shows, together:
- the kidney trajectory;
- the psychiatric history and relapse risk;
- the options and their consequences, including slow tapering, which carries less rebound risk than abrupt cessation [E4].

---

## 5. The engine

The prototype is in `lithos/` (about 1,500 lines of Python, 25 tests). Its choices are deliberate.

| Component | Choice | Why |
|---|---|---|
| Structure | Linear 1- or 2-compartment model with first-order absorption; one depot per formulation (IR/SR); **piecewise-constant parameters**; exact superposition via eigen-decomposition | Lithium is unmetabolised, renally cleared and linear. Time-varying parameters represent interactions, pregnancy, drift and renal change exactly. Verified against an ODE solver to 1e-11. |
| Priors | Published population PK models (covariates: kidney function, size, age) [E37–E39]; **model averaging** across candidate models by Laplace evidence | No single published model fits every population; averaging is more robust than picking one [E51] |
| Individualisation | MAP estimate + Laplace posterior (optional importance-sampling correction) | Fast (tens of ms), well understood, auditable |
| **Clearance drift** | Each blood-test occasion has its own clearance deviation, modelled as the sum of a slow (renal ageing) and a fast (hydration, illness) Ornstein–Uhlenbeck process, i.e. a Gaussian process in time | Decades of therapy: old levels inform the present less than recent ones; forecast uncertainty grows with time since the last level; sudden departures are detectable |
| **Timing-aware error** | Level variance includes (dC/dt × timing SD)² plus a distribution-phase term that decays with time since the dose | Samples on steep parts of the curve are automatically trusted less; this is what makes any-time sampling safe |
| Decisions | Probability of target attainment over real tablet combinations; cap on P(level > 1.0); step limits; hysteresis to avoid churn | Recommendations are probabilistic and practical |
| Recheck timing | From the individual half-life posterior, not a fixed 5–7 days | Older adults and people with CKD need longer; young adults shorter |
| Diagnostics | Samples in the absorption phase; implausibly long since the last dose; large standardised residuals; clearance jumps | Converts data-quality problems into specific questions a clinician can act on |

**What it deliberately does not do (yet):**
- patient-facing autonomous dose changes;
- decisions without clinician sign-off;
- interpretation of sensor signals;
- anything when its own uncertainty is too wide to be useful. In that case it says so and recommends a better-timed sample.

**Learning system.** Priors should improve with use:
- refit population parameters, drift magnitudes and interaction effect sizes from pooled, consented, de-identified data (federated where needed);
- monitor prediction error per release and per subgroup (age, sex, kidney function, ethnicity, formulation);
- publish the model cards.

---

## 6. Li12e: an "INR moment" for lithium

The 12-hour rule is a human workaround. Without a model, the only way to compare levels is to sample everyone at the same point on the same curve. With a model, the *report* can be standardised instead of the *sample*.

**Proposal.** Pathology reports carry, next to the measured value:
- the **Li12e**: the model-standardised steady-state 12-hour-equivalent level on the current regimen, with a 90% interval;
- the inputs it used: dose, regimen, last-dose time and sample time;
- quality flags such as "absorption phase" or "not at steady state; projection".

The lab needs two extra facts at collection: the regimen, and when the last dose was taken. These can come from the collection-centre form, the patient's app or dispensing data.

Why it matters:
- It ends "was that a 12-hour level?".
- It makes levels comparable across time, labs and settings.
- It makes a **lithium time-in-therapeutic-range** (TTR) metric computable for quality benchmarking and trials, the way TTR works for warfarin [E47].
- It has a standards gap to fill: LOINC has lithium codes (14334-7) and a "trough" variant, but no 12-hour post-dose term. Real datasets also map lithium codes and units inconsistently [E52]. Sample-to-dose timing must be carried explicitly.

Two practical traps the design handles:
1. **Lithium-heparin tubes** (used in some home-collection kits) falsely raise lithium results. Kits must be lithium-free [E48].
2. **Once-daily vs twice-daily dosing** produce different 12-hour levels for the same total dose. Li12e is always defined *for the regimen actually taken*, and the engine can translate between regimens.

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
  Li12e back on report)          │  Rules (guidelines as data)    │◄── Hospital EMR (SMART on FHIR app)
                                 │  Orchestrator (tasks, recall)  │
 GP systems ◄── CDS Hooks cards ─│  Register & quality dashboard  │──► My Health Record (care plan, results)
 (order-select / order-sign /    │  Audit, model registry,        │
  medication-refill; recall)     │  post-market monitoring        │──► Psychiatrist / lithium clinic workspace
                                 └────────────────────────────────┘
```

- **Data model:**
  - FHIR `MedicationRequest`/`MedicationStatement` (regimen and patient-reported dose times)
  - `Observation` (LOINC 14334-7 plus timing extensions; eGFR; TSH; calcium)
  - `CarePlan` (the lithium plan: target, regimen, schedule, sick-day rules)
  - `Task` (every monitoring obligation, with owner and due date)
  - `Communication`
- **Integration points:**
  - Pathology HL7 v2 feeds.
  - GP software. In Australia that means Best Practice and MedicalDirector, via partner APIs and HL7 import, and CDS Hooks where supported.
  - Hospital EMRs via SMART on FHIR. DoseMeRx's FHIR integration with Oracle Health at an Australian hospital shows the path [E41].
  - Electronic prescribing and dispensing (Active Script List) for refills and new interacting drugs, where access is permitted.
  - My Health Record.
- **Deployment:** a stateless engine service (the prototype is effectively this), an orchestration service, and clinician and patient front-ends. Built as a system of record for *lithium-specific* state only; everything else stays in the source systems.

---

## 8. Safety, regulation and governance

| Configuration | Likely position (to be confirmed with regulatory counsel) [E53] |
|---|---|
| **A. Clinician-facing** recommendations from lab results and EHR data; transparent basis; clinician signs; not time-critical | Australia: may meet the TGA clinical-decision-support exemption (advises a health professional; processes no device signal; doesn't replace judgement). US: plausibly non-device CDS. EU: MDR Rule 11 means Class IIa at least, realistically IIb. |
| **B. Patient-facing** dose instructions or self-titration | A regulated device everywhere (US: most plausibly Class II "drug dose calculator", with insulin-dosing software as predicate; Australia: likely Class IIb) |
| **C. Ingesting POC, wearable or saliva-sensor signals** | Fails CDS exemptions; the sensor is its own IVD submission |

Strategy:
- Launch as configuration A.
- Build to IEC 62304, ISO 14971 and ISO 13485 from day one regardless.
- Treat patient-facing dosing and sensor integration as later, separately cleared modules.

Guardrails in the engine and workflow:
- **Hard stops:** triage pathways above set levels or with neurological symptoms. **Caps:** on P(level > 1.0) and on the size of dose steps.
- An explicit **"cannot recommend"** state.
- **Mandatory clinician sign-off**, and a full **audit trail** of inputs, model version and recommendation.
- **Validation per model release**, with **monitoring of prediction error, drift and subgroup performance** in the field, and incident reporting.
- **Privacy:**
  - Australian Privacy Principles and the My Health Records Act;
  - explicit consent for the patient app and for research reuse;
  - data minimisation (the platform stores lithium-relevant data only).

---

## 9. Evidence: how to prove it, in order

| Stage | Question | Design | Primary measures |
|---|---|---|---|
| 0. In silico (this repo) | Do the mechanisms work under stated assumptions, including a misspecified model? | Virtual trials (§11) | Tests to target; accuracy of Li12e; interaction outcomes |
| 1. Retrospective validation | How accurate is the engine on real outpatient data? | External validation on lithium-clinic data with *recorded* dose and sample times; several populations | Bias, precision, % within ±0.1 mmol/L, calibration of intervals; Li12e vs properly timed levels |
| 2. Prospective shadow mode | Would its recommendations have been right? | Runs silently alongside care in 2–3 services | Agreement with subsequent levels; clinician-rated appropriateness; alert burden |
| 3. Pragmatic trial | Does it improve care? | Stepped-wedge cluster RCT across services and GP practices | **Primary:** lithium TTR over 12 months. **Secondary:** tests per patient-year; time to target at initiation; levels > 1.2; toxicity admissions; guideline-concordant monitoring; patient burden and satisfaction; clinician confidence; eGFR slope (long-term); **new lithium starts among eligible patients** |
| 4. Adaptive monitoring trial | Can stable patients safely test less often? | Randomised non-inferiority (risk-adapted vs fixed intervals) | Levels > 1.0; toxicity; relapse; tests |
| 5. Standards | Make it normal | Work with guideline and lab bodies (RANZCP, ISBD/IGSLi, NICE, RCPA, LOINC) | Li12e reporting standard; lithium TTR as a quality indicator |

The trial outcome that matters most is the last secondary one: **does making lithium easy make clinicians prescribe it when it's indicated?**

---

## 10. Getting it used: adoption and money

**Why nobody has done this:**
- the drug is generic and cheap, so there is no sponsor;
- care is fragmented, so there is no single buyer;
- regulatory cost is real;
- the vancomycin precedent shows that even a mandate plus mature software yields slow uptake, with licence cost the top barrier in 2026 [E43].

**Wedges (in order):**
1. **Pathology interpretive reporting.** Labs add Li12e and a plain-language interpretation to every lithium result. Every requesting GP benefits without installing anything. Pathology in Australia is concentrated in a few large providers, so a handful of partnerships reach most practices.
2. **Specialist services and lithium clinics:** the full clinician workspace, register and orchestration. This is where evidence is generated and where the champions are.
3. **GP software:** interaction forecasts at prescribing and refill, and the recall register. It spreads through the practices that already receive Li12e reports.
4. **Patients:** the app (dose-time capture, sick-day button, results, mood and side-effect tracking), then home and POC testing for selected, trained patients, as with warfarin self-management [E46].

**Money:**
- An **open, public-good core**: engine, models, validation code and model cards.
- A **funded service layer** (integration, hosting, support, quality reporting) paid for by health services, primary health networks and pathology providers.
- Research funding for the trials, from medical research funds and philanthropy.
- Cost offsets: fewer wasted or repeated tests (a stability-based interval scheme cut tests by about 15% [E21]), fewer toxicity admissions and, if uptake rises, fewer relapses. The last is the big one.

---

## 11. What the prototype shows

See [`SIMULATION.md`](SIMULATION.md) for methods, full results and caveats. Virtual patients were simulated from a deliberately *different* "truth" model than the engine's prior. They missed doses, had blood drawn at realistic times and had noisy assays.

*(Headline results are inserted from `SIMULATION_RESULTS.md` after the full run.)*

---

## 12. Risks, failure modes and mitigations

| Risk | How it would hurt | Mitigation in the design |
|---|---|---|
| Wrong or missing dose times | Li12e mis-standardised | Patient taps and pharmacy data; timing-uncertainty in the error model; absorption-phase flags; "cannot interpret" state |
| Model misspecification (e.g. older adults, CKD, pregnancy, other ethnicities, SR formulations) | Biased forecasts | Model averaging; subgroup validation before release; drift terms absorb individual deviation; field monitoring; conservative caps |
| Over-trust in history when clearance changes | Stale forecasts | Two-timescale drift model; clearance-jump diagnostics; in simulation, history *hurt* when drift was underestimated, which is why drift must be estimated from real data (§11) |
| Automation bias | Clinicians rubber-stamp | Show inputs, curve and uncertainty; require sign-off; audit; design reviews with human-factors testing |
| Alert fatigue | Ignored warnings | Forecasts only when they change a decision; cards carry a pre-filled action |
| Monitoring less often by accident | Missed toxicity | Guideline floor is enforced; interval *extension* only inside a trial |
| Digital exclusion | Inequity | SMS/IVR/carer channels; collection-centre capture of dose time; no smartphone required |
| Liability and regulation | Slow launch | Configuration A first; quality system from day one; clear intended-use statement |
| Business model | Nobody pays | Lab and service wedges; public-good core; trial evidence of cost offsets |
| Data breach | Harm and trust | Minimisation; encryption; access audit; Australian hosting |

---

## 13. Roadmap

| Horizon | Deliverables |
|---|---|
| 0–6 months | Engine hardened; published priors implemented and model-averaged; retrospective validation on 1–2 lithium-clinic datasets with recorded times; clinician workspace for one specialist service (register, recall, Li12e, recommendations); dose-time capture by SMS/app |
| 6–18 months | Shadow-mode study; pathology partner pilot of Li12e reporting; CDS Hooks interaction forecasts in one GP software environment; pregnancy module with a perinatal service; human-factors and safety case; quality management system certification |
| 18–36 months | Stepped-wedge trial; patient app general release; POC/DBS integration pilot for rural and remote sites; standards work on Li12e and lithium TTR; guideline engagement; start the adaptive-monitoring trial |

---

## 14. Open research questions this platform would answer

1. How large and how fast is within-person drift in lithium clearance? (It sets how much history to trust.)
2. How accurate is Li12e across formulations (SR vs IR), twice-daily regimens and older adults?
3. Which exposure metric (12-hour level, AUC, peak, time above 1.0, cumulative dose) best predicts renal decline?
4. Are there *personal* therapeutic windows (n-of-1 exposure–response using standardised levels and mood data)?
5. Can a person-specific calibration make saliva or sweat lithium usable [E48]?
6. Is risk-adaptive monitoring non-inferior to fixed intervals?
7. Does frictionless monitoring increase appropriate lithium initiation?
