"""Population pharmacokinetic models used as Bayesian priors.

A population model says what a *typical* patient with given covariates looks
like (fixed effects), how much individuals differ from that typical patient
(between-subject variability, ``omega``), how much one person's clearance
wanders over months and years (``drift``) and how noisy a measured level is
(``sigma``). The Bayesian layer combines this prior with the patient's own
levels.

Every model carries its provenance. Parameters marked *provisional* are
physiologically reasonable placeholders for software development only; they
are not a validated basis for dosing.
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
    weight: float               # kg
    scr: float                  # serum creatinine, umol/L
    height: float = 170.0       # cm

    @property
    def crcl(self) -> float:
        """Cockcroft-Gault creatinine clearance, mL/min."""
        return crcl_cockcroft_gault(self.scr, self.age, self.sex, self.weight)

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
    omega: dict[str, float]                             # SD of log-normal BSV per parameter
    formulations: dict[str, Formulation]
    sigma_prop: float                                   # proportional residual SD
    sigma_add: float                                    # additive residual SD, mmol/L
    # Within-person clearance change = slow drift (ageing kidneys, CKD) + fast
    # drift (hydration, sodium intake, intercurrent illness). Each is an
    # Ornstein-Uhlenbeck process on log-CL: (stationary SD, correlation time).
    drift_slow: tuple[float, float] = (0.10, 730.0)
    drift_fast: tuple[float, float] = (0.12, 21.0)
    timing_sd_h: float = 0.5                            # default uncertainty in reported dose/sample times
    # Structural uncertainty soon after a dose (absorption/distribution phase),
    # as a proportional SD that decays with time since the dose. Large for
    # one-compartment models, which cannot represent the distribution phase.
    dist_phase_sd: float = 0.5
    dist_phase_tau_h: float = 2.5
    status: str = "provisional"                         # provisional | published | published-adapted | validated
    notes: str = ""
    omega_corr: dict[tuple[str, str], float] = field(default_factory=dict)

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
        val = {k: tv.get(k, 0.0) * np.exp(eta.get(k, 0.0)) for k in ("cl", "v1", "q", "v2")}
        ka_mult = np.exp(eta.get("ka", 0.0))
        return Params(
            cl=val["cl"] * np.exp(drift) * cl_mult,
            v1=val["v1"] * v_mult,
            q=val["q"],
            v2=val["v2"] * v_mult,
            ka={k: f.ka * ka_mult for k, f in self.formulations.items()},
            f={k: f.f for k, f in self.formulations.items()},
            tlag={k: f.tlag for k, f in self.formulations.items()},
        )


# ---------------------------------------------------------------------------
# Provisional development model (physiology-based placeholder).
# Lithium is filtered freely and ~80 % is reabsorbed with sodium in the
# proximal tubule, so renal lithium clearance is roughly 20-30 % of creatinine
# clearance; volume approximates total body water (~0.7 L/kg). IR tablets peak
# at ~1-2 h, sustained-release products at ~4-5 h.
# ---------------------------------------------------------------------------

def _provisional_typical(c: Covariates) -> dict[str, float]:
    crcl_lh = c.crcl * 0.06
    return {"cl": 0.25 * crcl_lh, "v1": 0.7 * c.weight}


PROVISIONAL = PopModel(
    name="provisional-physiological",
    reference="Physiology-based placeholder (CL = 0.25 x CrCL, V = 0.7 L/kg); development only.",
    typical=_provisional_typical,
    omega={"cl": 0.25, "v1": 0.25},
    formulations={"IR": Formulation(ka=1.2), "SR": Formulation(ka=0.35, f=0.95)},
    sigma_prop=0.10,
    sigma_add=0.03,
    status="provisional",
)

# ---------------------------------------------------------------------------
# A deliberately different "truth" for in-silico trials: two compartments,
# a non-proportional renal covariate, a sex effect, more between-subject
# variability and faster within-person drift than any engine prior assumes.
# Used only to generate virtual patients - never as a dosing prior.
# ---------------------------------------------------------------------------

def _truth_typical(c: Covariates) -> dict[str, float]:
    cl = 1.45 * (c.egfr_abs / 90.0) ** 0.75 * (0.9 if c.sex.upper().startswith("F") else 1.0)
    v = 0.62 * c.weight ** 0.9 * 70 ** 0.1
    return {"cl": cl, "v1": 0.45 * v, "q": 4.0, "v2": 0.55 * v}


TRUTH_PERTURBED = PopModel(
    name="truth-perturbed",
    reference="Simulation truth only (not a dosing prior).",
    typical=_truth_typical,
    omega={"cl": 0.32, "v1": 0.30},
    formulations={"IR": Formulation(ka=1.5), "SR": Formulation(ka=0.30, f=0.9)},
    sigma_prop=0.0,
    sigma_add=0.0,
    drift_slow=(0.12, 540.0),
    drift_fast=(0.15, 14.0),
    status="simulation-truth",
)

MODELS: dict[str, PopModel] = {m.name: m for m in (PROVISIONAL, TRUTH_PERTURBED)}


def get_model(name: str | None = None) -> PopModel:
    if name is None:
        return MODELS[DEFAULT_MODEL]
    return MODELS[name]


DEFAULT_MODEL = PROVISIONAL.name
