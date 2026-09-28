"""What the engine knows about one patient, on an hours time axis.

Convention: t = 0 is midnight at the start of day 0, so ``t % 24`` is clock
time. ``Case.t0`` optionally anchors the axis to a calendar date for reports.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timedelta
from typing import Callable, Sequence

import numpy as np

from .models import Covariates
from .pk import Dose
from .units import mg_to_mmol

HOURS_PER_YEAR = 24 * 365.25


@dataclass
class Level:
    time: float                     # h, when the blood was drawn
    value: float                    # mmol/L
    timing_sd: float | None = None  # h, uncertainty in the dose/sample timing (None -> model default)
    note: str = ""


@dataclass
class Effect:
    """A known multiplicative change in clearance and/or volume over a window
    (a co-medication, pregnancy, an intercurrent illness...)."""

    start: float
    end: float | None = None
    cl_mult: float = 1.0
    v_mult: float = 1.0
    label: str = ""

    def active(self, t: float) -> bool:
        return self.start <= t and (self.end is None or t < self.end)


@dataclass
class CovariateTimeline:
    """Demographics plus a creatinine series, carried forward in time."""

    sex: str
    age_at_t0: float
    weight: float
    height: float = 170.0
    creatinine: Sequence[tuple[float, float]] = ((0.0, 80.0),)  # (t_h, umol/L)
    weights: Sequence[tuple[float, float]] = ()                  # optional (t_h, kg)
    # "nearest": use the measurement closest in time (creatinine is usually drawn
    # with the lithium level, so it describes the days *before* the draw as
    # well as after). "locf": last observation carried forward.
    mode: str = "nearest"

    def __call__(self, t: float) -> Covariates:
        pick = _nearest if self.mode == "nearest" else _locf
        scr = pick(self.creatinine, t)
        wt = pick(self.weights, t) if self.weights else self.weight
        return Covariates(age=self.age_at_t0 + t / HOURS_PER_YEAR, sex=self.sex, weight=wt,
                          scr=scr, height=self.height)

    @property
    def change_times(self) -> list[float]:
        out = set()
        for series in (self.creatinine, self.weights):
            ts = sorted(t for t, _ in series)
            if self.mode == "nearest":
                out |= {(a + b) / 2.0 for a, b in zip(ts[:-1], ts[1:])}
            else:
                out |= set(ts)
        return sorted(out)


def _nearest(series: Sequence[tuple[float, float]], t: float) -> float:
    pts = sorted(series)
    ts = [a for a, _ in pts]
    # ties (exact midpoints) resolve to the later measurement, matching change_times
    i = int(np.argmin([abs(a - t) - 1e-9 * (a >= t) for a in ts]))
    return pts[i][1]


def _locf(series: Sequence[tuple[float, float]], t: float) -> float:
    pts = sorted(series)
    val = pts[0][1]
    for ti, vi in pts:
        if ti <= t:
            val = vi
        else:
            break
    return val


@dataclass
class Case:
    covariates: Callable[[float], Covariates] | CovariateTimeline
    doses: list[Dose]
    levels: list[Level] = field(default_factory=list)
    effects: list[Effect] = field(default_factory=list)
    t0: datetime | None = None

    def cl_mult(self, t: float) -> float:
        return float(np.prod([e.cl_mult for e in self.effects if e.active(t)] or [1.0]))

    def v_mult(self, t: float) -> float:
        return float(np.prod([e.v_mult for e in self.effects if e.active(t)] or [1.0]))

    @property
    def change_times(self) -> list[float]:
        ts = set(getattr(self.covariates, "change_times", []))
        for e in self.effects:
            ts.add(e.start)
            if e.end is not None:
                ts.add(e.end)
        return sorted(ts)

    def last_dose_before(self, t: float) -> Dose | None:
        prior = [d for d in self.doses if d.time <= t]
        return max(prior, key=lambda d: d.time) if prior else None

    def when(self, t: float) -> str:
        if self.t0 is None:
            return f"day {int(t // 24)} {int(t % 24):02d}:{int(round((t % 1) * 60)) % 60:02d}"
        return (self.t0 + timedelta(hours=t)).strftime("%a %d %b %Y %H:%M")


# ---------------------------------------------------------------------------
# Regimen helpers
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Administration:
    clock: float          # hour of day, 0-24
    mg: float             # mg of salt
    formulation: str = "IR"
    salt: str = "carbonate"

    @property
    def mmol(self) -> float:
        return mg_to_mmol(self.mg, self.salt)


def expand(start_h: float, end_h: float, admins: Sequence[Administration],
           skip: Sequence[float] = ()) -> list[Dose]:
    """Expand a daily regimen into dose events in [start_h, end_h).
    ``skip`` lists dose times (h) that were missed."""
    doses = []
    day0 = int(np.floor(start_h / 24.0))
    day1 = int(np.ceil(end_h / 24.0))
    skipped = {round(s, 3) for s in skip}
    for day in range(day0, day1 + 1):
        for a in admins:
            t = day * 24.0 + a.clock
            if start_h <= t < end_h and round(t, 3) not in skipped:
                doses.append(Dose(t, a.mmol, a.formulation))
    return doses


def regimen_cycle(admins: Sequence[Administration]) -> list[tuple[float, float, str]]:
    """Convert administrations to the (clock, mmol, formulation) form used by
    :func:`lithos.pk.steady_state`."""
    return [(a.clock, a.mmol, a.formulation) for a in admins]


def li12_clock(admins: Sequence[Administration] | Sequence[tuple[float, float, str]]) -> float:
    """Clock time of the standardised 12-hour level: 12 h after the last dose
    of the day (for twice-daily dosing that is the pre-morning-dose trough)."""
    clocks = [a.clock if isinstance(a, Administration) else a[0] for a in admins]
    evening = max(clocks)
    return (evening + 12.0) % 24.0
