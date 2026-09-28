# Evidence brief

Compiled 27–28 September 2026 by three parallel research passes (pharmacokinetics/dosing; clinical guidelines and safety; digital landscape, sensing and regulation), then consolidated. The full working notes, about 30,000 words with 300+ references, are summarised here. Tags **[E1]**–**[E70]** are cited from [`DESIGN.md`](DESIGN.md) and the code.

## How far to trust this brief

The research environment's network policy blocked most primary sources: PubMed/PMC, publishers, NICE, RANZCP, TGA and FDA sites. The session's web-search quota was also exhausted part-way through. Verification therefore relied on verbatim public copies of primary texts hosted on GitHub (archived NICE guidance, FDA label XML, open-access full texts, abstract corpora, LOINC files, regulatory texts) and on search excerpts.

Every entry carries one of these labels:

| Label | Meaning |
|---|---|
| **V** | Read in the primary source: full text, label or verbatim abstract |
| **X** | Read in a search-engine excerpt of the named source |
| **S** | Read only in a secondary source citing it |
| **U** | Recollection, not confirmed this session; treat as a lead |

**RANZCP 2020, CANMAT/ISBD, BAP 2016, AGNP, Maudsley and the Australian Therapeutic Guidelines/AMH could not be read.** Guideline rules in the prototype follow NICE CG185 (verbatim) and ISBD/IGSLi 2019 (abstract). Confirming against RANZCP 2020 is a priority before any Australian use.

---

## A. Why lithium, and why it is under-used

| # | Finding | Source | Label |
|---|---|---|---|
| E1 | Maintenance efficacy vs placebo: relapse RR 0.66 (any episode), 0.52 (mania), 0.78 (depression). Longest time to treatment failure among common maintenance drugs in UK records: lithium 2.05 y vs quetiapine 0.76, valproate 0.98, olanzapine 1.13 | Severus et al. 2014; Hayes et al. 2016 | S |
| E2 | Anti-suicide signal: suicide OR 0.13 (0.03–0.66) in an RCT meta-analysis. For balance, a VA RCT in people with recent suicidal behaviour found no added benefit (HR 1.10, 0.77–1.55), with low achieved levels (0.54/0.46) | Cipriani et al. 2013 (S); Katz et al. 2022 (V, abstract) | S / V |
| E3 | Levels matter: standard range (median 0.83) vs low range (0.54) gave relapse 13% vs 38% (RR 2.6, 1.3–5.2), with more side effects at the standard range. **Confound:** many low-range patients had their level lowered abruptly at randomisation, and a reanalysis attributed much of the excess relapse to that abrupt change | Gelenberg et al., *NEJM* 1989 (V, abstract); Perlis et al., *Am J Psychiatry* 2002 (S) | V / S |
| E4 | Stopping carries risk: rebound mania; elevated suicidal behaviour after discontinuation. NICE advises tapering over at least 4 weeks, preferably up to 3 months | Suppes 1991; Tondo 1998; Baldessarini 1999; Yerevanian 2007 (S); NICE CG185 (V) | S / V |
| E5 | Declining use "worldwide in the last 20 years". US lithium share of bipolar visits fell from 30.4% to 17.6% | PMID 38057894 (V, abstract); Rhee et al. 2020 (U) | V / U |
| E6 | Discontinuation is common: in a Swedish cohort 54% stopped, 62% of them because of adverse effects. Prescribers across 24 countries monitor but use varying targets. **No quantified surveys of prescriber barriers were found.** This is a gap worth filling with primary research | Swedish cohort (S); Nederlof et al. 2018 (X) | S / X |

## B. Guidelines, targets and the 12-hour convention

| # | Finding | Source | Label |
|---|---|---|---|
| E7 | **NICE CG185 schedule.**<br>• Baseline: weight/BMI, U&E incl. calcium, eGFR, TFT, FBC; ECG if cardiovascular risk.<br>• Levels 1 week after starting and after every dose change, then weekly until stable; 3-monthly in year 1; then 6-monthly.<br>• Stay 3-monthly for: older people; interacting drugs; risk of impaired renal/thyroid function, raised calcium or other complications; poor symptom control; poor adherence; last level ≥ 0.8.<br>• Weight/BMI, U&E with calcium, eGFR and TSH every 6 months.<br>• Avoid NSAIDs. If prescribed, use them regularly (not as needed), with levels monthly until stable, then 3-monthly.<br>• Check for neurotoxicity at every review. | NICE CG185 (archived verbatim copy) | V |
| E8 | NICE targets: 0.6–0.8 mmol/L for people starting lithium; 0.8–1.0 for a trial of at least 6 months after relapse on lithium or with subthreshold symptoms. NG222 (depression): "12 hours post dose"; consider 0.4–0.6 at age 65+ for augmentation | NICE CG185; NG222 | V |
| E9 | ISBD/IGSLi task force:<br>• Standard maintenance 0.60–0.80.<br>• 0.40–0.60 for good response with poor tolerability; 0.80–1.00 for poor response with good tolerability.<br>• Older adults (majority view): 0.40–0.60, maximum 0.70–0.80 at 65–79 and 0.70 above 80. | Nolen et al., *Bipolar Disord* 2019 | V (abstract) |
| E10 | US labelling:<br>• Draw samples 12 h after the previous dose (one ER label says 8–12 h).<br>• Targets: acute 0.8–1.2, maintenance 0.8–1.0; toxic at ≥ 1.5.<br>• Not recommended at CrCl < 30 mL/min.<br>• Keep normal salt and 2.5–3 L of fluid during stabilisation; reduce or suspend lithium with protracted sweating, diarrhoea or febrile infection.<br>• Conversion: 300 mg carbonate = 8 mEq = 5 mL oral solution. | FDA SPL labels | V |
| E11 | Amdisen (1977) defined the standardised 12-hour serum lithium: drawn 12 h (±30 min) after the last dose, noting later samples are less valid. A 2025 proposal calls for re-introducing Amdisen standardisation | Amdisen, *Clin Pharmacokinet* 1977 (S); Mendelsohn, Aftab, Muir, *Bipolar Disord* 2025 (S) | S |
| E12 | RANZCP 2020 mood-disorders guideline (Malhi et al.): content **not accessed** | *ANZJP* 2021;55:7–117 | — |

## C. How monitoring actually goes, and what has helped

| # | Finding | Source | Label |
|---|---|---|---|
| E13 | England/Wales patient-safety alert "Safer lithium therapy" (2009): information booklet, alert card and record book; check tests before prescribing and dispensing; systems to catch interacting drugs | NPSA/2009/PSA005 | X |
| E14 | POMH-UK national audit.<br>• Baseline (patients on lithium ≥ 1 year): NICE standards met for serum level 30%, renal 55%, thyroid 50%.<br>• After the quality programme: four levels a year 30% → 48%; renal 55% → 70%; thyroid 49% → 66%.<br>• Change took years and was "generally modest". | Collins et al. 2010; Paton et al. 2013 | X |
| E15 | Laboratory data from three UK regions (46,555 requests): 19.2% of levels below range, 6.1% above; **69.4% of tests in in-range patients requested outside expected intervals** | *BMC Psychiatry* 2021 | X |
| E16 | Secondary care, south London: 50.7% of levels in range, 42.4% below, 6.9% above | Nikolova et al. 2018 | X |
| E17 | Fully guideline-compliant monitoring was about 16% in Dutch (n = 1,583) and Swedish (n = 4,428) cohorts; 21% of Swedish starts had no baseline creatinine. In Japan 15% had an annual level even after a 2012 regulatory warning. US VA: 19% had a repeat level within 6 months | Various | S |
| E18 | EHR-reminder RCT (KONOTORI, n = 111): levels 0.4–1.0 in 69.1% vs 60.0%, OR 2.14 (0.82–5.58), p = 0.12 (not significant) | Seki et al., *JMIR* 2023 | X |
| E19 | **Risk-adapted intervals:** after 12 months stable at 0.40–0.79, 6-monthly testing is reasonable; at 0.80–0.99 stay 3-monthly (10% vs 2% later reached ≥ 1.0); about 15% fewer tests | "Can we check serum lithium levels less often…", *BJPsych Open* 2022 | X |
| E20 | Norfolk county register with active recall: no known toxicity from inadequate monitoring over about 10 years. UK GP incentives reached 89–93% recorded levels, sustained after withdrawal, but higher achievement was not linked to fewer admissions | Kirkham et al. 2013 (V statement); QOF analysis (V) | V |
| E21 | **Interacting-drug initiation and toxicity (Kaiser Permanente Colorado):**<br>• 3,115 users; 70 episodes needing acute care (2.2%).<br>• Starting an ACE inhibitor, ARB, diuretic or NSAID in the prior month: about 30× odds of toxicity.<br>• **An interaction alert had fired for every interacting drug.** | Heath et al., *Psychiatr Serv* 2018 | X |
| E22 | **Timing reality:** about half of real lithium levels are not 12-hour levels.<br>• Bipolar CHOICE: 44.9% outside 10–14 h.<br>• Central Denmark, 52,837 tests: 49.7% outside, median 13.7 h, worst in general practice. | Jacobsen, Köhler-Forsberg et al., *Bipolar Disord* 2025 | S |

## D. Safety: interactions, kidneys, pregnancy, toxicity

| # | Finding | Source | Label |
|---|---|---|---|
| E23 | People aged 66+ (n = 10,615; 3.9% admitted for toxicity). RR of admission within 1 month of starting: **ACE inhibitor 7.6 (2.6–22.0)**; **loop diuretic 5.5 (1.9–16.1)**. Thiazides and NSAIDs were not independently associated | Juurlink et al., *J Am Geriatr Soc* 2004 | S |
| E24 | ACE inhibitors: steady-state levels +36.1%, clearance −25.5% (n = 20; 4 with toxicity symptoms). Thiazides: +25–40% (review). HCTZ raised levels more than furosemide or placebo in a crossover study | Finley et al. 1996 (S); Finley 2016 review (U); Crabtree et al. 1991 (V, abstract) | S / U / V |
| E25 | NSAID class labelling: mean trough +15%, renal lithium clearance −20%. Celecoxib +17%; meloxicam about +20%; indometacin and piroxicam "significantly" higher; sulindac sparing | FDA labels (V); meloxicam study (V, abstract) | V |
| E26 | SGLT2 inhibitors may **decrease** lithium levels | US lithium label 2026 (V); NZ Medsafe 2023, Malaysia NPRA 2024 (S) | V / S |
| E27 | **EXTRIP** (verbatim).<br>• ECTR **recommended** if kidney function is impaired and [Li⁺] > 4.0 mEq/L, or with decreased consciousness, seizures or life-threatening dysrhythmias irrespective of [Li⁺].<br>• **Suggested** if [Li⁺] > 5.0, significant confusion, or **expected time to [Li⁺] < 1.0 with optimal management > 36 h**.<br>• Continue until clinical improvement or [Li⁺] < 1.0 (minimum 6 h if not measurable); haemodialysis preferred.<br>• Very-low-quality evidence. | Decker et al., *CJASN* 2015 | V (abstract) |
| E28 | Long-term kidneys (IGSLi, 312 patients, 8–48 y of treatment):<br>• 29.5% had ≥ 1 eGFR < 60; no end-stage kidney failure.<br>• eGFR fell about 30% faster than with ageing alone.<br>• Risk factors include **higher serum levels**, longer treatment, older age and lower starting eGFR. | Tondo et al. 2017 | V |
| E29 | A single level of 1.01–1.2 was followed by an eGFR fall of 3.9%, and 1.21–2.0 by 5.4%, over ≤ 3 months (not detectable by 6 months) | Kirkham et al. 2014 (Norfolk, n = 699) | V |
| E30 | UK laboratory cohort, lithium vs controls: CKD stage 3 HR 1.93; hypothyroidism HR 2.31; hypercalcaemia HR 1.43 | Shine et al., *Lancet* 2015 | S |
| E31 | **NICE NG203 (CKD) rules:**<br>• A trajectory needs ≥ 3 eGFRs over ≥ 90 days; repeat a new low within 2 weeks to exclude AKI.<br>• Accelerated progression: sustained fall ≥ 25% with a category change within 12 months, or ≥ 15/yr.<br>• Refer if KFRE 5-year risk > 5% or ACR ≥ 70, among other criteria.<br>• GFR at least annually on lithium.<br>KDIGO rapid progression: > 5/yr. NICE CG185: when deciding whether to continue, seek renal and bipolar-specialist advice. | NICE NG203; KDIGO 2012; NICE CG185 | V |
| E32 | Renal-protective practice: once-daily dosing, lowest effective level, prevent intoxications, psychiatry–nephrology collaboration, amiloride for nephrogenic DI. Renal decline was minimal in older adults when levels were kept < 0.7 | Schoot et al. 2020 (S); Rej et al. 2020 (S) | S |
| E33 | **Pregnancy, Wesseloo et al. 2017** (1,101 levels, 113 pregnancies): level vs preconception −24% (T1), **−36% (T2)**, −21% (T3); +9% postpartum | *Br J Psychiatry* 2017 | S (corroborated via Poels, V) |
| E34 | **Pregnancy, Westin et al. 2017:** dose-adjusted levels fall 1.2% per gestational week (−7% at week 6, −22% at 20, −34% at 34); "doses generally need to be increased by 50%" in T3; steep postpartum rise by day 4 | *BMJ Open* 2017 | V (full text) |
| E35 | **Perinatal management.**<br>• GFR rises up to 50% by the end of T1.<br>• Weekly levels in T3; levels before and 24 h after delivery; restart on the first postpartum evening at 0.8–1.0; twice-weekly levels for 2 weeks (Poels 2018, V).<br>• NICE CG192: 4-weekly, weekly from 36 weeks, within 24 h of birth (X).<br>• US label: reduce or stop 2–3 days before delivery (V). Molenaar 2021 (233 levels): no intrapartum rise in dose-corrected levels, so they advise **against** pre-delivery reduction (S).<br>• Imaz 2025 (1,260 levels): −23.9/−27.6/−16.9% by trimester, +11% postpartum (S). | as listed | V / X / S |
| E36 | Toxicity kinetics: chronic poisoning half-life 36–79 h vs 19–29 h acute-on-chronic. Haemodialysis half-life 3.5–5.7 h, with rebound. Activated charcoal does not bind lithium | Case series (V, abstracts); US label (V) | V |

## E. Pharmacokinetics and models

| # | Finding | Source | Label |
|---|---|---|---|
| E37 | Systematic review of 8 lithium popPK models: typical CL 0.41–9.39 L/h; IIV 12.7–25.1%. **Renal function and body size** were the key covariates | Methaneethorn, *Eur J Drug Metab Pharmacokinet* 2018 | S |
| E38 | Methaneethorn & Sringam 2019 (222 Thai adults, acute mania):<br>• CL/F = 1.43 × (WT/65)^0.425 × (age/38)^−0.242 L/h.<br>• V/F 54 L and ka 0.426 /h (both fixed).<br>• IIV and residual error not retrieved. | *Hum Psychopharmacol* 2019 | V (abstract) / S |
| E39 | **External evaluation of 10 published models:** median absolute prediction error 24–138%; 7/10 under-predicted; only one model (with fat-free mass) acceptable on both datasets. The authors blame missing renal covariates. A **meta-model** (FFM model + GFR model) met acceptability in all subgroups | Lereclus et al., *Pharmaceuticals* 2023; *J Psychopharmacol* 2024 | S |
| E40 | Swedish popPK + GWAS (n = 2,357): age, sex, eGFR, diuretics and RAAS agents explained 61.4% / 49.8% of lithium-clearance variance | Millischer et al., *Lancet Psychiatry* 2022 | S |
| E41 | **Physiology.**<br>• Label: bioavailability ~100%; V 0.7–1 L/kg; no metabolism; 80% of filtered lithium reabsorbed proximally; half-life 18–36 h.<br>• Lithium clearance ≈ 20% of CrCL.<br>• Obese Vss 0.42 vs 0.66 L/kg.<br>• Older adults: distribution half-life 2.7 h, **distribution complete 10.6 ± 3.0 h**.<br>• IR→CR: peak 1.44 → 4.25 h, Cmax −22%, recovery 94.5% vs 90.2%.<br>• Renal lithium clearance is lower at night; half-life lengthens with years of therapy. | US label; HSDB/AMA; primary PK studies | V |
| E42 | Other models:<br>• Jermain 1991 (CL 1.36 L/h, V 32.8 L; LBW + CrCL).<br>• Yukawa 1993 (IIV 25.1%, residual 14.3%).<br>• Taright 1994 (nonparametric Bayesian).<br>• ElDesoky 2008.<br>• Pérez-Castelló 2016 (2-cmt; models non-adherence).<br>• Jin 2022 (CrCL + total daily dose).<br>• Landersdorfer 2017 (paediatric; adult-range after size scaling).<br>• Yoshida 2018 (V 0.79 L/kg, ka 1.5 /h).<br>• Couffignal 2019 (serum–erythrocyte). | as cited | mixed |
| E43 | **Bayesian lithium forecasting dates to the late 1980s.** A priori equations (Cooper, Pepin, Zetin, Terao, Abou-Auda) are biased and imprecise in comparisons: all three a priori models under-predicted in one study; Pepin's mean error was 37% | *Clin Pharmacokinet* 1989; PMID 3411433; Radhakrishnan 2012; Nichols 2014 (read with its corrigendum) | X / S |
| E44 | Machine learning (Taiwan, inpatients): RMSE 0.17–0.20 mmol/L; accuracy (±0.2) 0.68–0.75 | Hsu et al., *Biomedicines* 2021 | V (full text) |
| E45 | One-compartment benchmark calculator predicted the required dose within one 200 mg tablet (n = 58, retrospective) | Kavanagh et al., *J Psychopharmacol* 2025 | X |
| E46 | **PharmOracle:** app with model-informed dosing, Bayesian refinement and reminders. RCT protocol only (2026) | PMID 42529656 | X |

## F. Prior art on standardising mistimed levels

| # | Finding | Source | Label |
|---|---|---|---|
| E47 | **eLi12** (estimated 12-hour serum lithium): mean deviation from the true 12-h level 10% vs 25% for the raw level; closer in 97% of samples drawn 3–24 h post dose. Depends on accurate recall of dose time. A feasibility trial is registered (NCT07306039) | Köhler-Forsberg, Wiuff, Devantier, Østergaard, Faraone, Nierenberg et al., *J Clin Psychiatry* 2025 | S |
| E48 | **SimpLi:** empirical multiplicative standardisation explained little variance (R² 0.108 with time only; 0.177 with dose), which supports mechanistic, dose-history-aware methods | medRxiv 2026 | S |
| E49 | Once- vs twice-daily dosing: the 12-h level is reported 10–26% higher on once-daily for the same dose, but this has never been confirmed prospectively (NCT03811860). On once-daily SR, the 12-h level was 1.3× the 24-h level | Perry 1981 (U); Reddy 2014 (S) | U / S |

## G. Digital landscape and analogues

| # | Finding | Source | Label |
|---|---|---|---|
| E50 | **No mainstream precision-dosing product lists lithium.**<br>• Explicitly absent from NextDose and from Tucuxi's public model repository.<br>• Not found for DoseMeRx, InsightRx, PrecisePK, BestDose, TDMx, ID-ODS or AutoKinetics; MwPharm++ unknown.<br>DoseMe history:<br>• Brisbane; founder/CTO Robert McLeay.<br>• Acquired by Tabula Rasa HealthCare in 2019, then sold to a Fairlong Capital affiliate in 2023.<br>• CE and TGA listed; Oracle Health FHIR integration at Bayside Health (the Alfred), July 2026. | Vendor pages and filings via excerpts; Tucuxi repo (V) | X / V |
| E51 | **Vancomycin.**<br>• The 2020 consensus named Bayesian software as preferred and dropped trough-only monitoring.<br>• Nephrotoxicity 8% → 0–2% with Bayesian AUC dosing (Neely 2018).<br>• A single-sample Bayesian approach saved about US$2,065 per patient.<br>• Only 42.8% of surveyed US hospitals had adopted it a year later; licence cost is the top barrier (65%) in 2026. | Rybak et al. 2020; Neely 2018; Lee 2021; Bradley 2021; *Front Pharmacol* 2026 survey | X |
| E52 | **Warfarin.** Self-testing and self-management (28 RCTs, 8,950 people) gave fewer thromboembolic events and lower mortality in selected, trained patients. The Rosendaal TTR metric | Heneghan et al., *Lancet* 2006 / Cochrane 2016; Rosendaal 1993 | X |
| E53 | Clozapine "no blood, no drug" registries; the US REMS was discontinued in 2025 | Secondary sources | S |
| E54 | Model averaging/selection improves precision-dosing predictions (vancomycin case study) | Uster et al., *Clin Pharmacol Ther* 2021 | U |

## H. Near-patient and novel sensing

| # | Finding | Source | Label |
|---|---|---|---|
| E55 | **Finger-prick capillary electrophoresis (Medimate/FISIC)**:<br>• CE-certified for self-testing.<br>• r = 0.96 vs laboratory, with bias and limits of agreement within ±0.2 mmol/L (*BJPsych* 2026); r = 0.917 in perinatal women.<br>• The company has been bankrupt twice.<br>• A UK randomised feasibility trial (LiPOC, n = 80) is under way. | as cited | X |
| E56 | InstaRead finger-stick: r 0.93–0.98; FDA-cleared and CLIA-waived about 2004; current availability unverified | Glazer et al. 2004 | X |
| E57 | Saliva: r = 0.77 after adjustment; **individual saliva/serum ratios are stable**, so person-specific calibration is plausible | Parkin et al., *Bipolar Disord* 2021; Vlahović 2025 | X |
| E58 | Sweat patches: pilot data only (3 patients); no paired sweat–serum validation. No clinical microneedle lithium sensor | USC 2025; EPFL 2018/2021 | X |
| E59 | Volumetric dried blood spots agree well with serum and allow postal home sampling | *J Pharm Biomed Anal* 2023 and others | X |
| E60 | **Lithium-heparin tubes falsely elevate lithium**, and some home-collection kits ship with them, so kits must be lithium-free | Laboratory practice; Tasso+ kit configuration | U (practice) |
| E61 | Digital phenotyping: MONARCA I (smartphone self-monitoring) showed no difference in primary outcomes. **No validated link from passive signals (tremor, sleep, typing) to lithium levels or toxicity** | Faurholt-Jepsen et al. 2015, 2016 | V |

## I. Interoperability and regulation

| # | Finding | Source | Label |
|---|---|---|---|
| E62 | LOINC:<br>• 14334-7 Lithium [Moles/volume] in Serum or Plasma; 3719-2 (mass).<br>• Trough variants 74806-1 and 3723-4.<br>• **No 12-hour post-dose term.**<br>• Real datasets map lithium codes and units inconsistently. | LOINC files; OMOP/MIMIC mappings | V |
| E63 | CDS Hooks: `order-select`, `order-sign` and `medication-refill` (the latter may fire with no user in context); `medication-prescribe` is deprecated | CDS Hooks specification | V |
| E64 | **Regulation.**<br>• US statute excludes some clinician-facing decision support from the device definition.<br>• FDA CDS guidance (29 Jan 2026) extends enforcement discretion to a single clinically appropriate recommendation; a patient-facing insulin calculator is a device (Example 28).<br>• Drug-dose calculators are FDA product code NDC, Class II.<br>• EU MDR Rule 11: IIa minimum, IIb/III where decisions can cause serious or irreversible harm.<br>• TGA exemption for CDS that advises a health professional, processes no device signal and does not replace judgement. | Statute, guidance and regulation texts (V); TGA criteria (S) | V / S |

## J. Professor Gin Malhi's lithium work (bibliographic; no full texts read)

| # | Citation | Label |
|---|---|---|
| E65 | Malhi GS, Bell E, Outhred T, Berk M. Lithium therapy and its interactions. *Aust Prescr* 2020;43:91–3. doi:10.18773/austprescr.2020.024. A correction (12 June 2020) clarified that non-thiazide diuretics may *alter*, not necessarily raise, levels | S |
| E66 | Malhi GS, Bell E, Porter RJ, et al. Lithium should be borne in mind: five key reasons. *ANZJP* 2020;54(7):659–63 | S |
| E67 | Malhi GS, Masson M, Bellivier F (eds). *The Science and Practice of Lithium Therapy.* Springer, 2017. Not to be confused with Malhi GS, Tanious M, Das P, Berk M, *ANZJP* 2012;46:192–211, which has the same title | S |
| E68 | Malhi GS, Tanious M, Gershon S. The lithiumeter: a measured approach. *Bipolar Disord* 2011;13:219–26: a framework for level targets | S |
| E69 | Malhi GS, Tanious M, Bargh D, Das P, Berk M. Safe and effective use of lithium. *Aust Prescr* 2013;36:18–21. Also Malhi GS, *Lancet* 2012;379:690–2 ("Is the safety of lithium no longer in the balance?") | S |
| E70 | Malhi GS, Bell E, Bassett D, et al. The 2020 RANZCP clinical practice guidelines for mood disorders. *ANZJP* 2021;55(1):7–117; bipolar summary *Bipolar Disord* 2020;22:805–21. Secondary sources cite him arguing that lithium's decline in favour of antipsychotics is ill-advised, and that end-stage renal failure after long-term therapy is "possible, but largely unlikely" | S |
