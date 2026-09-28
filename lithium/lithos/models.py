"""Population pharmacokinetic models used as Bayesian priors.

A population model says what a *typical* patient with given covariates looks
like (fixed effects), how much individuals differ from that typical patient
(between-subject variability, ``omega``), how much one person's clearance
wanders over months and years (``drift``) and how noisy a measured level is
(``sigma``). The Bayesian layer combines this prior with the patient's own
levels.

Why an ensemble. An external evaluation of ten published lithium models
found poor predictive performance on new data (median absolute prediction
errors of 24-138%), attributed largely to missing kidney-function covariates;
a meta-model combining a fat-free-mass model with a GFR model did better
(Lereclus et al. 2023, 2024; docs/EVIDENCE.md). So the default engine weighs
several priors by how well each explains the patient's own levels rather
than trusting any single model.

Every model carries its provenance and a status. Nothing here is validated
for clinical use.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from .pk import Params
from .renal import bsa_du_bois, crcl_cockcroft_gault, egfr_ckd_epi_2021, fat_free_mass


@dataclass(frozen=True)
class Covariates:
    age: float                  # years
    sex: str                    # "F" or "M"
    weight: float               # kg (total body weight)
    scr: float                  # serum creatinine, umol/L
    height: float = 170.0       # cm

    @property
    def female(self) -> bool:
        return self.sex.upper().startswith("F")

    @property
    def bmi(self) -> float:
        return self.weight / (self.height / 100.0) ** 2

    @property
    def ibw(self) -> float:
        """Devine ideal body weight, kg."""
        return (45.5 if self.female else 50.0) + 0.9 * (self.height - 152.4)

    @property
    def dosing_weight(self) -> float:
        """Total body weight, or adjusted body weight when BMI > 30 (Cockcroft-Gault
        with total weight overestimates creatinine clearance in obesity)."""
        if self.bmi <= 30:
            return self.weight
        return self.ibw + 0.4 * (self.weight - self.ibw)

    @property
    def crcl(self) -> float:
        """Cockcroft-Gault creatinine clearance, mL/min (adjusted weight if obese)."""
        return crcl_cockcroft_gault(self.scr, self.age, self.sex, self.dosing_weight)

    @property
    def egfr(self) -> float:
        """CKD-EPI 2021 eGFR, mL/min/1.73 m^2 (what laboratories report)."""
        return egfr_ckd_epi_2021(self.scr, self.age, self.sex)

    @property
    def egfr_abs(self) -> float:
        """eGFR de-indexed to the patient's body surface area, mL/min."""
        return self.egfr * bsa_du_bois(self.weight, self.height) / 1.73

    @property
    def ffm(self) -> float:
        return fat_free_mass(self.weight, self.height, self.sex)


@dataclass(frozen=True)
class Formulation:
    ka: float          # 1/h
    f: float = 1.0     # bioavailability relative to oral solution
    tlag: float = 0.0  # h


@dataclass
class PopModel:
    name: str
    reference: str
    typical: Callable[[Covariates], dict[str, float]]   # -> cl, v1 (and q, v2 for 2-cmt)
    omega: dict[str, float]                             # SD of log-normal BSV; "v" scales V1 and V2 together
    formulations: dict[str, Formulation]
    sigma_prop: float                                   # proportional residual SD
    sigma_add: float                                    # additive residual SD, mmol/L
    # Within-person clearance change = slow drift (ageing kidneys, CKD) + fast
    # drift (hydration, sodium intake, intercurrent illness). Each is an
    # Ornstein-Uhlenbeck process on log-CL: (stationary SD, correlation time in days).
    drift_slow: tuple[float, float] = (0.10, 730.0)
    drift_fast: tuple[float, float] = (0.08, 21.0)
    timing_sd_h: float = 0.5                            # default uncertainty in reported dose/sample times
    # Structural uncertainty soon after a dose (absorption/distribution phase),
    # as a proportional SD that decays with time since the dose. Large for
    # one-compartment models, which cannot represent the distribution phase.
    dist_phase_sd: float = 0.5
    dist_phase_tau_h: float = 2.5
    status: str = "provisional"     # provisional | published | published-adapted | simulation-truth
    notes: str = ""
    omega_corr: dict[tuple[str, str], float] = field(default_factory=dict)
    # Who the model was built for; outside it the model is left out of the ensemble.
    applicable: Callable[[Covariates], bool] | None = None

    def applies_to(self, cov: Covariates) -> bool:
        return self.applicable is None or bool(self.applicable(cov))

    def drift_cov(self, dt_days):
        """Covariance of log-CL drift between two times ``dt_days`` apart."""
        dt = np.abs(np.asarray(dt_days, dtype=float))
        (s1, t1), (s2, t2) = self.drift_slow, self.drift_fast
        return s1 * s1 * np.exp(-dt / t1) + s2 * s2 * np.exp(-dt / t2)

    @property
    def eta_names(self) -> list[str]:
        return list(self.omega)

    def omega_matrix(self) -> np.ndarray:
        names = self.eta_names
        sd = np.array([self.omega[k] for k in names])
        corr = np.eye(len(names))
        for (a, b), r in self.omega_corr.items():
            i, j = names.index(a), names.index(b)
            corr[i, j] = corr[j, i] = r
        return corr * np.outer(sd, sd)

    def individual(
        self,
        cov: Covariates,
        eta: dict[str, float] | None = None,
        drift: float = 0.0,
        cl_mult: float = 1.0,
        v_mult: float = 1.0,
        typical: dict[str, float] | None = None,
    ) -> Params:
        """Individual parameters: typical value x exp(eta) x drift x external effects."""
        eta = eta or {}
        tv = typical if typical is not None else self.typical(cov)
        ev = np.exp(eta.get("v", 0.0))
        return Params(
            cl=tv["cl"] * np.exp(eta.get("cl", 0.0) + drift) * cl_mult,
            v1=tv["v1"] * np.exp(eta.get("v1", 0.0)) * ev * v_mult,
            q=tv.get("q", 0.0) * np.exp(eta.get("q", 0.0)),
            v2=tv.get("v2", 0.0) * np.exp(eta.get("v2", 0.0)) * ev * v_mult,
            ka={k: f.ka * np.exp(eta.get("ka", 0.0)) for k, f in self.formulations.items()},
            f={k: f.f for k, f in self.formulations.items()},
            tlag={k: f.tlag for k, f in self.formulations.items()},
        )


# Absorption: IR tablets peak at 0.25-3 h, SR/CR at 2-6 h (US label); in a
# head-to-head study peak time moved from 1.44 h (IR) to 4.25 h (CR) with 22%
# lower Cmax, and urinary recovery was 94.5% vs 90.2% (docs/EVIDENCE.md).
ABSORPTION = {
    "IR": Formulation(ka=1.2),
    "SR": Formulation(ka=0.40, f=0.95),
    "LIQ": Formulation(ka=2.5),
}


# ---------------------------------------------------------------------------
# Model A: renal-physiological two-compartment prior (published-adapted).
# Assembled from verified relationships rather than fitted to one dataset:
#   * lithium clearance ~20-23.5% of creatinine clearance (AMA/HSDB; Pepin);
#     23% chosen, giving 1.38 L/h at CrCL 100 mL/min, consistent with
#     published typical values of 1.36-1.50 L/h in adults with normal kidneys;
#   * lithium distributes in total body water: V 0.7-1.0 L/kg in lean adults
#     but ~0.42 L/kg in obesity, so total volume scales with fat-free mass
#     (0.85 L/kg FFM, calibrated against published obese/lean and older-adult
#     cohorts - see validation/literature.py);
#   * distribution half-life ~1.4 h in adults and ~2.7 h in older adults, with
#     distribution complete ~10.6 h post dose in the elderly: V1 55% / V2 45%
#     of total (a larger central volume reproduces the published CR/IR peak
#     ratio and older-adult peaks), Q 5.5 L/h scaled by FFM^0.75 and falling
#     with age;
#   * between-subject variability on CL ~24-25% in published models.
# Residual error excludes timing error and drift, which are modelled
# separately, so it is smaller than published "residual" terms (14-16%).
# ---------------------------------------------------------------------------

def _renal_typical(c: Covariates) -> dict[str, float]:
    vt = 0.85 * c.ffm
    return {
        "cl": 0.23 * 0.06 * c.crcl,
        "v1": 0.55 * vt,
        "v2": 0.45 * vt,
        "q": 5.5 * (c.ffm / 55.0) ** 0.75 * (max(c.age, 18.0) / 40.0) ** -1.0,
    }


RENAL_2CMT = PopModel(
    name="lithos-renal-2cmt",
    reference=("Assembled from published physiology: CL = 23% of Cockcroft-Gault CrCL; V = 0.85 L/kg "
               "fat-free mass, 55% central / 45% peripheral; distribution half-life 1.4-2.7 h. "
               "See docs/EVIDENCE.md and validation/literature.py."),
    typical=_renal_typical,
    omega={"cl": 0.25, "v": 0.20},
    formulations=ABSORPTION,
    sigma_prop=0.10,
    sigma_add=0.02,
    dist_phase_sd=0.15,
    dist_phase_tau_h=3.0,
    status="published-adapted",
)


# ---------------------------------------------------------------------------
# Model B: Methaneethorn & Sringam 2019 (Hum Psychopharmacol 34:e2697).
# 222 Thai adults with acute mania; 1-compartment; CL/F = 1.43 x (WT/65)^0.425
# x (age/38)^-0.242 L/h; V/F 54 L and ka 0.426 /h fixed. The published IIV and
# residual error could not be retrieved, so they are ASSUMED here (IIV within
# the 12.7-25.1% range across published models; residual near Yukawa's 14.3%).
# No kidney-function covariate, which is why it gets a lower prior weight.
# ---------------------------------------------------------------------------

def _methaneethorn_typical(c: Covariates) -> dict[str, float]:
    return {"cl": 1.43 * (c.weight / 65.0) ** 0.425 * (max(c.age, 18.0) / 38.0) ** -0.242, "v1": 54.0}


METHANEETHORN_2019 = PopModel(
    name="methaneethorn-2019",
    reference="Methaneethorn J, Sringam S. Hum Psychopharmacol 2019;34:e2697 (PMID 31025773); IIV/residual assumed.",
    typical=_methaneethorn_typical,
    omega={"cl": 0.20, "v1": 0.20},
    formulations={k: Formulation(ka=0.426, f=f.f) for k, f in ABSORPTION.items()},
    sigma_prop=0.12,
    sigma_add=0.02,
    status="published (IIV and residual assumed)",
    # Built in adults with acute mania (mean age ~38) and without a kidney covariate:
    # not used for older adults or anyone with reduced kidney function.
    applicable=lambda c: c.egfr >= 60 and c.age < 65,
)


# ---------------------------------------------------------------------------
# Development placeholder kept for tests and comparison.
# ---------------------------------------------------------------------------

def _provisional_typical(c: Covariates) -> dict[str, float]:
    return {"cl": 0.25 * 0.06 * c.crcl, "v1": 0.7 * c.weight}


PROVISIONAL = PopModel(
    name="provisional-physiological",
    reference="Physiology-based placeholder (CL = 0.25 x CrCL, V = 0.7 L/kg); development only.",
    typical=_provisional_typical,
    omega={"cl": 0.25, "v1": 0.25},
    formulations=ABSORPTION,
    sigma_prop=0.10,
    sigma_add=0.03,
    drift_fast=(0.12, 21.0),
    status="provisional",
)


# ---------------------------------------------------------------------------
# A deliberately different "truth" for in-silico trials: two compartments,
# a non-proportional renal covariate on de-indexed eGFR, a sex effect,
# volume on total (not fat-free) weight, more between-subject variability and
# larger, faster within-person drift than any engine prior assumes.
# Used only to generate virtual patients - never as a dosing prior.
# ---------------------------------------------------------------------------

def _truth_typical(c: Covariates) -> dict[str, float]:
    cl = 1.45 * (c.egfr_abs / 90.0) ** 0.75 * (0.9 if c.female else 1.0)
    v = 0.62 * c.weight ** 0.9 * 70 ** 0.1
    return {"cl": cl, "v1": 0.45 * v, "q": 4.0, "v2": 0.55 * v}


TRUTH_PERTURBED = PopModel(
    name="truth-perturbed",
    reference="Simulation truth only (not a dosing prior).",
    typical=_truth_typical,
    omega={"cl": 0.32, "v1": 0.30},
    formulations={"IR": Formulation(ka=1.5), "SR": Formulation(ka=0.30, f=0.9), "LIQ": Formulation(ka=3.0)},
    sigma_prop=0.0,
    sigma_add=0.0,
    drift_slow=(0.12, 540.0),
    drift_fast=(0.12, 14.0),
    status="simulation-truth",
)

MODELS: dict[str, PopModel] = {m.name: m for m in (RENAL_2CMT, METHANEETHORN_2019, PROVISIONAL, TRUTH_PERTURBED)}

# Default engine: ensemble of priors with prior model probabilities. The
# kidney-aware model gets most of the prior weight; the data decide the rest.
DEFAULT_ENSEMBLE: list[tuple[str, float]] = [(RENAL_2CMT.name, 0.75), (METHANEETHORN_2019.name, 0.25)]
DEFAULT_MODEL = RENAL_2CMT.name


def get_model(name: str | None = None) -> PopModel:
    return MODELS[name or DEFAULT_MODEL]


def get_engine(spec: str | None = None) -> list[tuple[PopModel, float]]:
    """``None`` or "ensemble" -> the default ensemble; a model name -> that model alone."""
    if spec in (None, "ensemble"):
        return [(MODELS[n], w) for n, w in DEFAULT_ENSEMBLE]
    return [(MODELS[spec], 1.0)]
