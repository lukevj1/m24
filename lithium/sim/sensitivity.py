"""Sensitivity analysis: does an older level help when the drift model is right?

In the main any-time experiment, adding a level from 45 days earlier made the
engine slightly *less* accurate. The hypothesis is that the virtual patients'
clearance drifts more than the engine assumes, so the engine over-trusts old
data. Here the engine's drift terms are set to the truth's and the experiment
is repeated.

    python -m sim.sensitivity --n 300
"""

from __future__ import annotations

import argparse
import dataclasses
import multiprocessing as mp

import numpy as np

from lithos.models import get_engine, get_model

from . import trials
from .population import make_population


def _engine(matched: bool):
    eng = get_engine()
    if not matched:
        return eng
    truth = get_model("truth-perturbed")
    return [(dataclasses.replace(m, drift_slow=truth.drift_slow, drift_fast=truth.drift_fast), w) for m, w in eng]


_POP: dict = {}


def _one(args):
    i, n, matched, seed = args
    truth = get_model("truth-perturbed")
    if (n, seed) not in _POP:
        _POP[(n, seed)] = make_population(n, truth, seed=seed)
    vp = _POP[(n, seed)][i]
    return trials.any_time(vp, truth, _engine(matched), np.random.default_rng(777 + i))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=300)
    ap.add_argument("--procs", type=int, default=mp.cpu_count())
    ap.add_argument("--seed", type=int, default=4242)
    a = ap.parse_args()
    for matched in (False, True):
        with mp.Pool(a.procs) as pool:
            res = pool.map(_one, [(i, a.n, matched, a.seed) for i in range(a.n)], chunksize=4)
        e1 = np.median([abs(r["engine_med"] - r["truth12"]) for r in res])
        e2 = np.median([abs(r["hist_med"] - r["truth12"]) for r in res])
        label = "drift matched to truth" if matched else "engine's own drift"
        print(f"{label:24s} this level only {e1:.3f} | + earlier level {e2:.3f} mmol/L (n={len(res)})")


if __name__ == "__main__":
    main()
