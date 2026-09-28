"""Bayesian individualisation: MAP estimation with a Laplace posterior.

Three choices set this engine apart from a textbook MAP dosing program:

1. **Clearance is allowed to drift.** Lithium is taken for decades while renal
   function, hydration, sodium intake and co-medication change. Each blood
   test occasion gets its own clearance deviation; deviations follow the sum
   of a slow and a fast Ornstein-Uhlenbeck process (a Gaussian process whose
   correlation decays with the time between tests). Old levels inform the
   present, but less than recent ones; forecasts use the GP conditional, so
   uncertainty grows with the time since the last level; and a sudden
   departure from a patient's own history is detectable.

2. **Timing error is part of the error model.** The variance of each level
   includes (dC/dt x timing SD)^2, so a sample drawn on the steep part of the
   curve (soon after a dose) is automatically trusted less than one drawn on
   the flat part. This is what makes "any-time" sampling safe to interpret.

3. **Uncertainty is carried to the decision.** The Laplace posterior (with an
   optional importance-sampling correction) is propagated to every forecast,
   so recommendations come with probabilities, not just point estimates.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import optimize

from .case import Case, Level
from .models import PopModel
from .pk import Params, simulate

OCCASION_WINDOW_H = 72.0   # levels within 3 days of each other share an occasion
LOG2PI = np.log(2 * np.pi)


def _occasions(levels: list[Level]) -> tuple[np.ndarray, np.ndarray]:
    """Group levels into occasions. Returns (occasion time per group, group index per level)."""
    if not levels:
        return np.zeros(0), np.zeros(0, dtype=int)
    times = np.array([lv.time for lv in levels])
    order = np.argsort(times)
    group = np.empty(times.size, dtype=int)
    starts: list[float] = []
    members: list[list[float]] = []
    for i in order:
        if not starts or times[i] - starts[-1] > OCCASION_WINDOW_H:
            starts.append(times[i])
            members.append([])
        members[-1].append(times[i])
        group[i] = len(starts) - 1
    return np.array([np.mean(m) for m in members]), group


@dataclass
class Posterior:
    model: PopModel
    case: Case
    names: list[str]
    mode: np.ndarray
    cov: np.ndarray
    occasion_times: np.ndarray
    ofv: float
    log_evidence: float
    converged: bool
    diagnostics: dict = field(default_factory=dict)

    # -- structure --------------------------------------------------------
    @property
    def n_eta(self) -> int:
        return len(self.model.eta_names)

    @property
    def n_levels(self) -> int:
        return len(self.case.levels)

    def eta(self, u: np.ndarray) -> dict[str, float]:
        return dict(zip(self.model.eta_names, u[: self.n_eta]))

    def _boundaries(self) -> np.ndarray:
        o = self.occasion_times
        return (o[1:] + o[:-1]) / 2.0

    def drift_fitted(self, u: np.ndarray, t: float) -> float:
        if self.occasion_times.size == 0:
            return 0.0
        k = int(np.searchsorted(self._boundaries(), t, side="right"))
        return float(u[self.n_eta + k])

    # -- forecasting ------------------------------------------------------
    def draw(self, n: int, t: float, rng: np.random.Generator | None = None,
             sir: bool = False) -> list[tuple[dict[str, float], float]]:
        """Posterior draws of (eta, clearance drift at time t).

        For t after the last occasion the drift relaxes toward the patient's
        long-run mean, and gains fresh uncertainty, as an OU process would.
        """
        rng = rng or np.random.default_rng(0)
        U = self.sample(n, rng, sir=sir)
        occ = self.occasion_times
        k_tt = float(self.model.drift_cov(0.0))
        if occ.size == 0:
            drifts = rng.normal(0.0, np.sqrt(k_tt), size=n)
        elif t <= occ[-1]:
            drifts = np.array([self.drift_fitted(u, t) for u in U])
        else:
            # Gaussian-process conditional of the drift at t given the occasion drifts.
            K = self.model.drift_cov((occ[:, None] - occ[None, :]) / 24.0)
            k_star = self.model.drift_cov((t - occ) / 24.0)
            w = np.linalg.solve(K, k_star)
            var = max(k_tt - float(k_star @ w), 0.0)
            drifts = U[:, self.n_eta:] @ w + rng.normal(0.0, np.sqrt(var), size=n)
        return [(self.eta(u), float(d)) for u, d in zip(U, drifts)]

    def params_at(self, eta: dict[str, float], drift: float, t: float, cl_mult: float = 1.0) -> Params:
        cov = self.case.covariates(t)
        return self.model.individual(cov, eta, drift, self.case.cl_mult(t) * cl_mult, self.case.v_mult(t))

    def param_draws(self, n: int, t: float, rng: np.random.Generator | None = None,
                    cl_mult: float | np.ndarray = 1.0, sir: bool = False) -> list[Params]:
        """Posterior draws of the individual parameters in force at time t
        (optionally with an extra, possibly per-draw, clearance multiplier)."""
        mult = np.broadcast_to(np.asarray(cl_mult, dtype=float), (n,))
        return [self.params_at(eta, d, t, cl_mult=float(m)) for (eta, d), m in zip(self.draw(n, t, rng, sir), mult)]

    def predict_draws(self, times, n: int, rng: np.random.Generator | None = None, doses=None) -> np.ndarray:
        """Posterior draws of the concentration-time curve (n x len(times))."""
        return np.array([self.predict(times, u, doses=doses) for u in self.sample(n, rng)])

    @property
    def summary_model(self) -> str:
        return f"{self.model.name} [{self.model.status}]"

    def sample(self, n: int, rng: np.random.Generator | None = None, sir: bool = False) -> np.ndarray:
        rng = rng or np.random.default_rng(0)
        if not sir or self.n_levels == 0:
            return rng.multivariate_normal(self.mode, self.cov, size=n, method="cholesky")
        # Sampling-importance-resampling from an inflated Laplace proposal.
        infl = 1.5
        prop = rng.multivariate_normal(self.mode, self.cov * infl, size=4 * n, method="cholesky")
        inv = np.linalg.inv(self.cov * infl)
        dev = prop - self.mode
        logq = -0.5 * np.einsum("ij,jk,ik->i", dev, inv, dev)
        logp = np.array([-0.5 * self.diagnostics["_ofv_fn"](u) for u in prop])
        logw = logp - logq
        w = np.exp(logw - logw.max())
        w /= w.sum()
        self.diagnostics["sir_ess"] = float(1.0 / np.sum(w * w))
        return prop[rng.choice(prop.shape[0], size=n, replace=True, p=w)]

    def predict(self, times, u: np.ndarray | None = None, doses=None) -> np.ndarray:
        """Individual prediction at ``times`` for parameter vector ``u`` (default:
        the MAP estimate), using the case's doses or an alternative dose list.
        Beyond the last occasion the latest clearance state is held constant."""
        u = self.mode if u is None else u
        times = np.atleast_1d(np.asarray(times, float))
        if doses is None:
            return self.diagnostics["_predict_fn"](u, times)[0]
        return simulate(doses, self.diagnostics["_schedule_fn"](u), times)

    def clearance_history(self) -> list[tuple[float, float]]:
        """Relative clearance at each occasion vs the patient's long-run mean."""
        return [(float(t), float(np.exp(self.mode[self.n_eta + k])))
                for k, t in enumerate(self.occasion_times)]

    def clearance_ratio(self, t: float) -> float:
        """Estimated clearance relative to what the model predicts from kidney
        function and size (1.0 = as expected). Much above 1 with low levels
        suggests missed doses as often as fast kidneys."""
        eta = self.eta(self.mode)
        return float(np.exp(eta.get("cl", 0.0) + self.drift_fitted(self.mode, t)))


def fit(model: PopModel, case: Case, *, max_iter: int = 200) -> Posterior:
    """Fit the individual model to the case by MAP estimation."""
    levels = sorted(case.levels, key=lambda lv: lv.time)
    occ_times, _ = _occasions(levels)
    n_eta, n_occ = len(model.eta_names), occ_times.size
    names = [f"eta_{k}" for k in model.eta_names] + [f"drift_{k}" for k in range(n_occ)]

    # Prior covariance: between-subject variability, then OU-correlated drift.
    P = np.zeros((n_eta + n_occ, n_eta + n_occ))
    P[:n_eta, :n_eta] = model.omega_matrix()
    if n_occ:
        P[n_eta:, n_eta:] = model.drift_cov((occ_times[:, None] - occ_times[None, :]) / 24.0)
    P_inv = np.linalg.inv(P)
    _, logdet_2piP = np.linalg.slogdet(2 * np.pi * P)

    # Piecewise-constant segments: occasion boundaries, covariate and effect changes.
    bounds = list((occ_times[1:] + occ_times[:-1]) / 2.0)
    first_t = min([d.time for d in case.doses] + [lv.time for lv in levels] + [0.0])
    breaks = sorted({b for b in bounds + case.change_times if b > first_t})
    starts = [first_t] + breaks
    seg_eval = [s + 1e-6 for s in starts]
    seg_typical = [model.typical(case.covariates(s)) for s in seg_eval]
    seg_cov = [case.covariates(s) for s in seg_eval]
    seg_clm = [case.cl_mult(s) for s in seg_eval]
    seg_vm = [case.v_mult(s) for s in seg_eval]
    seg_occ = [int(np.searchsorted(np.array(bounds), s, side="right")) for s in seg_eval]

    y = np.array([lv.value for lv in levels])
    t_obs = np.array([lv.time for lv in levels])
    tsd = np.array([lv.timing_sd if lv.timing_sd is not None else model.timing_sd_h for lv in levels])
    dose_t = np.sort([d.time for d in case.doses])
    since = np.array([t - dose_t[np.searchsorted(dose_t, t, side="right") - 1]
                      if np.searchsorted(dose_t, t, side="right") > 0 else np.inf for t in t_obs])
    dist_sd = model.dist_phase_sd * np.exp(-since / model.dist_phase_tau_h)

    def schedule(u):
        eta = dict(zip(model.eta_names, u[:n_eta]))
        out = []
        for s, tv, cv, cm, vm, k in zip(starts, seg_typical, seg_cov, seg_clm, seg_vm, seg_occ):
            drift = u[n_eta + k] if n_occ else 0.0
            out.append((s, model.individual(cv, eta, drift, cm, vm, typical=tv)))
        out[0] = (-np.inf, out[0][1])
        return out

    def predict_fn(u, times):
        return simulate(case.doses, schedule(u), times, derivative=True)

    def ofv(u):
        prior = float(u @ P_inv @ u) + logdet_2piP
        if not levels:
            return prior
        f, df = predict_fn(u, t_obs)
        var = model.sigma_add ** 2 + ((model.sigma_prop ** 2 + dist_sd ** 2) * f * f) + (df * tsd) ** 2
        return float(np.sum((y - f) ** 2 / var + np.log(var) + LOG2PI)) + prior

    d = n_eta + n_occ
    u0 = np.zeros(d)
    converged = True
    if levels:
        res = optimize.minimize(ofv, u0, method="BFGS", options={"gtol": 1e-6, "maxiter": max_iter})
        u_hat = res.x if np.isfinite(res.fun) else u0
        if not res.success:
            # BFGS often stops on "precision loss" at a genuine optimum; polish and keep the better point.
            nm = optimize.minimize(ofv, u_hat, method="Nelder-Mead",
                                   options={"maxiter": 4000, "xatol": 1e-7, "fatol": 1e-9})
            if np.isfinite(nm.fun) and nm.fun <= ofv(u_hat):
                u_hat = nm.x
        converged = _max_abs_gradient(ofv, u_hat) < 1e-3
        H = _hessian(ofv, u_hat)
        A = 0.5 * H
        try:
            cov = np.linalg.inv(A)
            np.linalg.cholesky(cov)
        except np.linalg.LinAlgError:
            A = A + np.eye(d) * 1e-6
            cov = np.linalg.inv(A)
            converged = False
        ofv_hat = ofv(u_hat)
        _, logdet_A = np.linalg.slogdet(A)
        log_evidence = -0.5 * ofv_hat + 0.5 * d * LOG2PI - 0.5 * logdet_A
    else:
        u_hat, cov, ofv_hat, log_evidence = u0, P, ofv(u0), 0.0

    post = Posterior(model, case, names, u_hat, cov, occ_times, ofv_hat, float(log_evidence), converged)
    post.diagnostics["_ofv_fn"] = ofv
    post.diagnostics["_predict_fn"] = predict_fn
    post.diagnostics["_schedule_fn"] = schedule
    post.diagnostics["levels"] = _level_diagnostics(post, levels, t_obs, y, tsd, dist_sd)
    return post


def _max_abs_gradient(f, x: np.ndarray, h: float = 1e-5) -> float:
    g = [(f(x + h * e) - f(x - h * e)) / (2 * h) for e in np.eye(x.size)]
    return float(np.max(np.abs(g))) if g else 0.0


def _hessian(f, x: np.ndarray, h: float = 1e-3) -> np.ndarray:
    d = x.size
    H = np.zeros((d, d))
    f0 = f(x)
    for i in range(d):
        ei = np.zeros(d)
        ei[i] = h
        H[i, i] = (f(x + ei) - 2 * f0 + f(x - ei)) / h ** 2
        for j in range(i + 1, d):
            ej = np.zeros(d)
            ej[j] = h
            H[i, j] = H[j, i] = (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej) + f(x - ei - ej)) / (4 * h * h)
    return 0.5 * (H + H.T)


def _level_diagnostics(post: Posterior, levels, t_obs, y, tsd, dist_sd) -> list[dict]:
    """Per-level data-quality checks a clinician can act on."""
    if not levels:
        return []
    model, case = post.model, post.case
    f, df = post.diagnostics["_predict_fn"](post.mode, t_obs)
    var = model.sigma_add ** 2 + ((model.sigma_prop ** 2 + dist_sd ** 2) * f * f) + (df * tsd) ** 2
    out = []
    for lv, fi, vi, dfi in zip(levels, f, var, df):
        last = case.last_dose_before(lv.time)
        since = None if last is None else lv.time - last.time
        flags = []
        if since is not None:
            early = 8.0 if last.formulation == "SR" else 6.0
            if since < early:
                flags.append(f"drawn {since:.1f} h after a dose (absorption/distribution phase): "
                             "informative only with an accurate dose time")
            if since > 30.0:
                flags.append(f"drawn {since:.0f} h after the last recorded dose: check for missed doses")
        z = (lv.value - fi) / np.sqrt(vi)
        if abs(z) > 3.0:
            flags.append(f"level is {abs(z):.1f} SD from the individual fit: check timing, adherence, "
                         "sample handling (e.g. lithium-heparin tube) or a change in clearance")
        out.append({"time": lv.time, "value": lv.value, "fitted": float(fi), "z": float(z),
                    "hours_since_dose": since, "slope_per_h": float(dfi), "flags": flags})
    return out


@dataclass
class Ensemble:
    """Bayesian model averaging over several population priors.

    Each member is fitted separately; its posterior probability is
    proportional to its prior weight times its Laplace marginal likelihood
    (evidence). Forecasts mix the members' posterior draws in proportion.
    """

    members: list[Posterior]
    weights: np.ndarray

    @property
    def case(self) -> Case:
        return self.members[0].case

    @property
    def best(self) -> Posterior:
        return self.members[int(np.argmax(self.weights))]

    @property
    def model(self) -> PopModel:
        return self.best.model

    @property
    def n_levels(self) -> int:
        return self.best.n_levels

    @property
    def diagnostics(self) -> dict:
        return self.best.diagnostics

    @property
    def converged(self) -> bool:
        return all(m.converged for m in self.members)

    @property
    def summary_model(self) -> str:
        parts = [f"{m.model.name} {w:.0%}" for m, w in sorted(zip(self.members, self.weights),
                                                             key=lambda mw: -mw[1])]
        return "ensemble: " + ", ".join(parts)

    def clearance_history(self) -> list[tuple[float, float]]:
        return self.best.clearance_history()

    def clearance_ratio(self, t: float) -> float:
        return float(np.exp(np.sum(self.weights * np.log([m.clearance_ratio(t) for m in self.members]))))

    def _split(self, n: int, rng: np.random.Generator) -> np.ndarray:
        return rng.multinomial(n, self.weights)

    def param_draws(self, n: int, t: float, rng: np.random.Generator | None = None,
                    cl_mult: float | np.ndarray = 1.0, sir: bool = False) -> list[Params]:
        rng = rng or np.random.default_rng(0)
        mult = np.broadcast_to(np.asarray(cl_mult, dtype=float), (n,))
        out: list[Params] = []
        i = 0
        for m, k in zip(self.members, self._split(n, rng)):
            if k:
                out += m.param_draws(int(k), t, rng, cl_mult=mult[i:i + k], sir=sir)
                i += k
        return out

    def predict_draws(self, times, n: int, rng: np.random.Generator | None = None, doses=None) -> np.ndarray:
        rng = rng or np.random.default_rng(0)
        parts = [m.predict_draws(times, int(k), rng, doses=doses)
                 for m, k in zip(self.members, self._split(n, rng)) if k]
        return np.vstack(parts)

    def predict(self, times, u=None, doses=None) -> np.ndarray:
        """Model-averaged point prediction (weighted mean of member MAP predictions)."""
        preds = np.array([m.predict(times, doses=doses) for m in self.members])
        return self.weights @ preds


def individualise(case: Case, engine) -> "Posterior | Ensemble":
    """Fit ``engine`` to ``case``.

    ``engine`` is a PopModel, or a list of (PopModel, prior weight) pairs for
    model averaging (see :func:`lithos.models.get_engine`).
    """
    if isinstance(engine, PopModel):
        return fit(engine, case)
    t_ref = max([lv.time for lv in case.levels] + [d.time for d in case.doses] + [0.0])
    cov = case.covariates(t_ref)
    usable = [(m, w) for m, w in engine if m.applies_to(cov)] or list(engine[:1])
    if len(usable) == 1:
        return fit(usable[0][0], case)
    engine = usable
    members = [fit(m, case) for m, _ in engine]
    prior = np.array([w for _, w in engine], dtype=float)
    logw = np.log(prior / prior.sum()) + np.array([m.log_evidence for m in members])
    w = np.exp(logw - logw.max())
    return Ensemble(members, w / w.sum())
