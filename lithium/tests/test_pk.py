import numpy as np
import pytest
from scipy.integrate import solve_ivp

from lithos.pk import Dose, Params, simulate, steady_state


def one_cmt(cl=1.4, v=50.0, ka=1.0):
    return Params(cl=cl, v1=v, ka={"IR": ka}, f={"IR": 1.0})


def test_single_dose_matches_closed_form():
    p = one_cmt()
    t = np.array([0.5, 1, 2, 5, 10, 24, 48])
    k, ka, D = 1.4 / 50, 1.0, 27.07
    exact = D * ka / (50 * (ka - k)) * (np.exp(-k * t) - np.exp(-ka * t))
    assert np.allclose(simulate([Dose(0, D)], p, t), exact, atol=1e-12)


def test_degenerate_ka_equals_k():
    p = one_cmt(cl=1.0, v=50.0, ka=0.02)
    t = np.array([10.0, 50.0])
    exact = 10 * 0.02 * t * np.exp(-0.02 * t) / 50
    assert np.allclose(simulate([Dose(0, 10.0)], p, t), exact, rtol=1e-3)


def _ode_reference(doses, sched, grid, forms):
    fidx = {f: i for i, f in enumerate(forms)}
    c = len(forms)

    def params_at(t):
        p = sched[0][1]
        for s, pp in sched:
            if t >= s:
                p = pp
        return p

    events = sorted(set([d.time for d in doses] + [s for s, _ in sched if np.isfinite(s)] + [grid[-1] + 1]))
    x, t, out = np.zeros(c + 2), -1.0, np.full(grid.size, np.nan)
    for ev in events:
        p = params_at(t + 1e-9)

        def rhs(_, y, p=p):
            dy = np.zeros(c + 2)
            for f, i in fidx.items():
                dy[i] = -p.ka[f] * y[i]
                dy[c] += p.ka[f] * y[i]
            dy[c] += -(p.cl + p.q) / p.v1 * y[c] + p.q / p.v2 * y[c + 1]
            dy[c + 1] = p.q / p.v1 * y[c] - p.q / p.v2 * y[c + 1]
            return dy

        sol = solve_ivp(rhs, (t, ev), x, dense_output=True, rtol=1e-11, atol=1e-13, method="DOP853")
        m = (grid >= t) & (grid < ev)
        if m.any():
            out[m] = sol.sol(grid[m])[c] / p.v1
        x, t = sol.y[:, -1].copy(), ev
        for d in doses:
            if d.time == ev:
                x[fidx[d.formulation]] += d.amount * params_at(ev + 1e-9).f.get(d.formulation, 1.0)
    return out


def test_two_compartment_mixed_formulations_time_varying_matches_ode():
    p1 = Params(cl=1.5, v1=20, q=3.0, v2=25, ka={"IR": 1.2, "SR": 0.35}, f={"IR": 1.0, "SR": 0.9})
    p2 = Params(cl=1.0, v1=22, q=3.0, v2=25, ka={"IR": 1.2, "SR": 0.35}, f={"IR": 1.0, "SR": 0.9})
    doses = [Dose(24 * d + 21, 13.5, "SR") for d in range(8)] + [Dose(24 * d + 8, 6.77, "IR") for d in range(8)]
    sched = [(0.0, p1), (100.0, p2)]
    grid = np.linspace(0, 200, 801)
    ref = _ode_reference(doses, sched, grid, ["IR", "SR"])
    assert np.nanmax(np.abs(ref - simulate(doses, sched, grid))) < 1e-8


def test_steady_state_equals_long_simulation():
    p = Params(cl=1.2, v1=18, q=2.5, v2=30, ka={"IR": 1.1}, f={"IR": 1.0})
    reg = [(8.0, 13.5, "IR"), (20.0, 13.5, "IR")]
    doses = [Dose(24 * d + h, a, f) for d in range(80) for (h, a, f) in reg]
    clock = np.array([0.0, 4.0, 8.0, 9.5, 12.0, 20.0, 23.9])
    assert np.allclose(simulate(doses, p, 79 * 24 + clock), steady_state(reg, p, clock), atol=1e-10)


def test_steady_state_auc_equals_dose_over_clearance():
    p = Params(cl=1.3, v1=20, q=3.0, v2=30, ka={"SR": 0.4}, f={"SR": 0.9})
    reg = [(21.0, 24.0, "SR")]
    t = np.linspace(0, 24, 20001)
    c = steady_state(reg, p, t)
    auc = np.sum((c[1:] + c[:-1]) / 2 * np.diff(t))
    assert auc == pytest.approx(0.9 * 24.0 / 1.3, rel=1e-4)


def test_terminal_half_life_two_compartment():
    p = Params(cl=1.2, v1=18, q=2.5, v2=30, ka={"IR": 1.5}, f={"IR": 1.0})
    t = np.array([200.0, 260.0])
    c = simulate([Dose(0, 20.0)], p, t)
    slope_hl = np.log(2) * (t[1] - t[0]) / np.log(c[0] / c[1])
    assert slope_hl == pytest.approx(p.terminal_half_life(), rel=1e-6)


def test_derivative_matches_finite_difference():
    p = Params(cl=1.2, v1=18, q=2.5, v2=30, ka={"IR": 1.5}, f={"IR": 1.0})
    doses = [Dose(24 * d + 21, 20.0) for d in range(5)]
    t = np.array([100.0, 110.3])
    c, dc = simulate(doses, p, t, derivative=True)
    h = 1e-5
    fd = (simulate(doses, p, t + h) - simulate(doses, p, t - h)) / (2 * h)
    assert np.allclose(dc, fd, rtol=1e-5, atol=1e-9)
